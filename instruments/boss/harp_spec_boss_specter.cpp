// @COPYRIGHT@

#include <harp_internal.hpp>


using namespace std;
using namespace harp;

static const char * format_boss_specter = "boss_specter";

static const char * boss_specter_key_path = "path";
static const char * boss_specter_key_objonly = "objonly";


harp::spec_boss_specter::spec_boss_specter ( boost::property_tree::ptree const & props ) : spec ( format_boss_specter, props ) {

  path_ = props.get < string > ( boss_specter_key_path );

  boost::optional < string > objval = props.get_optional < string > ( boss_specter_key_objonly );
  string obj = boost::get_optional_value_or ( objval, "FALSE" );
  objonly_ = ( obj == "TRUE" );

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

  nglobal_ = nspec_ * nlambda_;
  
}


harp::spec_boss_specter::~spec_boss_specter ( ) {
  
  cleanup();
  
}


boost::property_tree::ptree harp::spec_boss_specter::serialize ( ) {
  boost::property_tree::ptree ret;


  return ret;
}


void harp::spec_boss_specter::read ( matrix_dist & data, std::vector < double > & lambda, std::vector < bool > & sky ) {

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


void harp::spec_boss_specter::write ( std::string const & path, matrix_dist & data, std::vector < double > const & lambda, std::vector < bool > const & sky ) {

  HARP_THROW( "boss specter format spec writing not yet implemented" );

  return;
}



