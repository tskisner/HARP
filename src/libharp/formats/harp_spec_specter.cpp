// @COPYRIGHT@

#include <harp_internal.hpp>


using namespace std;
using namespace harp;



static const char * spec_specter_key_path = "path";
static const char * spec_specter_key_objonly = "objonly";
static const char * spec_specter_key_nspec = "nspec";
static const char * spec_specter_key_nlambda = "nlambda";


harp::spec_specter::spec_specter ( boost::property_tree::ptree const & props ) : spec ( props ) {

  //cerr << "spec specter props = " << endl;
  //ptree_print ( props );

  path_ = props.get < string > ( spec_specter_key_path, "" );

  boost::optional < string > objval = props.get_optional < string > ( spec_specter_key_objonly );
  string obj = boost::get_optional_value_or ( objval, "FALSE" );
  objonly_ = ( obj == "TRUE" );

  if ( path_ == "" ) {

    nspec_ = props.get < size_t > ( spec_specter_key_nspec );

    nlambda_ = props.get < size_t > ( spec_specter_key_nlambda );

    spechdu_ = 1;

    lambdahdu_ = 2;

    targethdu_ = 3;

  } else {

    int np;
    int myp;

    MPI_Comm_size ( MPI_COMM_WORLD, &np );
    MPI_Comm_rank ( MPI_COMM_WORLD, &myp );

    if ( myp == 0 ) {

      fitsfile * fp;

      fits::open_read ( fp, path_ );

      if ( objonly_ ) {
        spechdu_ = fits::img_seek ( fp, "EXTNAME", "OBJPHOT" );
      } else {
        spechdu_ = fits::img_seek ( fp, "EXTNAME", "FLUX" );
      }

      fits::img_dims ( fp, nspec_, nlambda_ );

      lambdahdu_ = fits::img_seek ( fp, "EXTNAME", "WAVELENGTH" );
      targethdu_ = fits::bin_seek ( fp, "EXTNAME", "TARGETINFO" );

      fits::close ( fp );

    }

    // broadcast 

    int ret = MPI_Bcast ( (void*)(&nspec_), 1, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD );
    mpi_check ( MPI_COMM_WORLD, ret );

    ret = MPI_Bcast ( (void*)(&nlambda_), 1, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD );
    mpi_check ( MPI_COMM_WORLD, ret );

    ret = MPI_Bcast ( (void*)(&spechdu_), 1, MPI_INT, 0, MPI_COMM_WORLD );
    mpi_check ( MPI_COMM_WORLD, ret );

    ret = MPI_Bcast ( (void*)(&lambdahdu_), 1, MPI_INT, 0, MPI_COMM_WORLD );
    mpi_check ( MPI_COMM_WORLD, ret );

    ret = MPI_Bcast ( (void*)(&targethdu_), 1, MPI_INT, 0, MPI_COMM_WORLD );
    mpi_check ( MPI_COMM_WORLD, ret );

  }

  nglobal_ = nspec_ * nlambda_;
  
}


harp::spec_specter::~spec_specter ( ) {
  
  cleanup();
  
}


boost::property_tree::ptree harp::spec_specter::serialize ( ) {
  boost::property_tree::ptree ret;


  return ret;
}


void harp::spec_specter::read ( matrix_dist & data, std::vector < double > & lambda, std::vector < bool > & sky ) {

  data.ResizeTo ( nglobal_, 1 );
  dist_matrix_zero ( data );

  lambda.resize ( nlambda_ );
  sky.resize ( nspec_ );

  int np;
  int myp;
  int ret;

  MPI_Comm_size ( MPI_COMM_WORLD, &np );
  MPI_Comm_rank ( MPI_COMM_WORLD, &myp );

  fitsfile * fp;

  // read and broadcast the spectral data

  matrix_local iobuffer ( nglobal_, 1 );

  if ( myp == 0 ) {
    fits::open_read ( fp, path_ );
    fits::img_seek ( fp, spechdu_ );   
    fits::img_read ( fp, iobuffer );
  }

  ret = MPI_Bcast ( (void*)(iobuffer.Buffer()), nglobal_, MPI_DOUBLE, 0, MPI_COMM_WORLD );
  mpi_check ( MPI_COMM_WORLD, ret );

  // Each process sets its local elements of the distributed spectra

  int hlocal = data.LocalHeight();
  int wlocal = data.LocalWidth(); // this should be one...

  int rowoff = data.ColShift();
  int rowstride = data.ColStride();
  int row;

  for ( int i = 0; i < hlocal; ++i ) {
    row = rowoff + i * rowstride;
    data.SetLocal ( i, 0, iobuffer.Get ( row, 0 ) );
  }

  // read and broadcast the wavelength vector

  iobuffer.ResizeTo ( nlambda_, 1 );

  if ( myp == 0 ) {
    fits::img_seek ( fp, lambdahdu_ );

    fits::img_read ( fp, iobuffer );

    for ( size_t w = 0; w < nlambda_; ++w ) {
      lambda[w] = iobuffer.Get ( w, 0 );
    }
  }

  ret = MPI_Bcast ( (void*)(&(lambda[0])), nlambda_, MPI_DOUBLE, 0, MPI_COMM_WORLD );
  mpi_check ( MPI_COMM_WORLD, ret );

  // read and broadcast the sky flag

  vector < uint8_t > skyflags ( nspec_ );

  if ( myp == 0 ) {
    fits::bin_seek ( fp, targethdu_ );
    vector < string > colnames ( 1 );
    vector < int > cols ( 1 );
    colnames[0] = "OBJTYPE";
    cols = fits::bin_columns ( fp, colnames );

    vector < string > objnames;
    fits::bin_read_strings ( fp, 0, nspec_ - 1, cols[0], objnames );

    for ( size_t i = 0; i < nspec_; ++i ) {
      if ( objnames[i] == "SKY" ) {
        skyflags[i] = 1;
      } else {
        skyflags[i] = 0;
      }
    }

    fits::close ( fp );
  }

  ret = MPI_Bcast ( (void*)(&(skyflags[0])), nspec_, MPI_CHAR, 0, MPI_COMM_WORLD );
  mpi_check ( MPI_COMM_WORLD, ret );

  for ( size_t i = 0; i < nspec_; ++i ) {
    if ( skyflags[i] == 0 ) {
      sky[i] = false;
    } else {
      sky[i] = true;
    }
  }

  return;
}


void harp::spec_specter::write ( std::string const & path, matrix_dist & data, std::vector < double > const & lambda, std::vector < bool > const & sky ) {

  fitsfile * fp;

  int ret;
  int status = 0;

  int np;
  int myp;

  MPI_Comm_size ( MPI_COMM_WORLD, &np );
  MPI_Comm_rank ( MPI_COMM_WORLD, &myp );

  // create file

  if ( myp == 0 ) {
    fits::create ( fp, path );
    fits::img_append ( fp, nspec_, nlambda_ );
    fits::write_key ( fp, "EXTNAME", "FLUX", "" );

    fits::img_append ( fp, 1, nlambda_ );
    fits::write_key ( fp, "EXTNAME", "WAVELENGTH", "" );

    matrix_local lambda_loc ( nlambda_, 1 );
    for ( size_t i = 0; i < nlambda_; ++i ) {
      lambda_loc.Set ( i, 0, lambda[i] );
    }
    fits::img_write ( fp, lambda_loc );

    char ** ttype = (char**) malloc ( 3 * sizeof(char*) );
    char ** tform = (char**) malloc ( 3 * sizeof(char*) );
    char ** tunit = (char**) malloc ( 3 * sizeof(char*) );
    for ( size_t i = 0; i < 3; ++i ) {
      ttype[i] = (char*) malloc ( FLEN_VALUE );
      tform[i] = (char*) malloc ( FLEN_VALUE );
      tunit[i] = (char*) malloc ( FLEN_VALUE );
    }
    strcpy ( ttype[0], "OBJTYPE" );
    strcpy ( tform[0], "6A" );
    strcpy ( tunit[0], "" );
    strcpy ( ttype[1], "Z" );
    strcpy ( tform[1], "1E" );
    strcpy ( tunit[1], "" );
    strcpy ( ttype[2], "O2FLUX" );
    strcpy ( tform[2], "1E" );
    strcpy ( tunit[2], "" );

    ret = fits_create_tbl ( fp, BINARY_TBL, nspec_, 3, ttype, tform, tunit, "TARGETINFO", &status );
    fits::check ( status );

    for ( size_t i = 0; i < 3; ++i ) {
      free ( ttype[i] );
      free ( tform[i] );
      free ( tunit[i] );
    }
    free ( ttype );
    free ( tform );
    free ( tunit );

    char ** charray;
    charray = (char**) malloc ( nspec_ * sizeof( char* ) );
    if ( ! charray ) {
      HARP_THROW( "cannot allocate temp char array" );
    }
    for ( size_t i = 0; i < nspec_; ++i ) {
      charray[i] = (char*) malloc ( 50 * sizeof (char) );
      if ( ! charray[i] ) {
        HARP_THROW( "cannot allocate temp char array member" );
      }
      if ( sky[i] ) {
        strcpy ( charray[i], "SKY" );
      } else {
        strcpy ( charray[i], "ELG" );
      }
    }

    ret = fits_write_col_str ( fp, 1, 1, 1, nspec_, charray, &status );
    fits::check ( status );

    for ( size_t i = 0; i < nspec_; ++i ) {
      free ( charray[i] );
    }
    free ( charray );


  }

  // reduce data to root process

  matrix_local data_loc ( data.Height(), 1 );
  local_matrix_zero ( data_loc );

  elem::AxpyInterface < double > globloc;
  globloc.Attach( elem::GLOBAL_TO_LOCAL, data );
  globloc.Axpy ( 1.0, data_loc, 0, 0 );
  globloc.Detach();

  // write data

  if ( myp == 0 ) {
    fits::img_seek ( fp, spechdu_ );
    fits::img_write ( fp, data_loc );
    fits::close ( fp );
  }

  return;
}



