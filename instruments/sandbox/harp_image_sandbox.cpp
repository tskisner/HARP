// @COPYRIGHT@

#include <harp_internal.hpp>


using namespace std;
using namespace harp;

static const char * format_sandbox = "sandbox";

static const char * sandbox_image_key_path = "path";
static const char * sandbox_image_key_signal = "signal";
static const char * sandbox_image_key_sky = "sky";
static const char * sandbox_image_key_noise = "noise";
static const char * sandbox_image_key_rows = "rows";
static const char * sandbox_image_key_cols = "cols";


harp::image_sandbox::image_sandbox ( boost::property_tree::ptree const & props ) : image ( props ) {

  sighdu_ = props.get ( sandbox_image_key_signal, 1 );

  nsehdu_ = props.get ( sandbox_image_key_noise, 2 );

  skyhdu_ = props.get ( sandbox_image_key_sky, 3 );

  path_ = props.get ( sandbox_image_key_path, "" );

  if ( path_ == "" ) {

    // no path specified- must specify rows / cols

    rows_ = props.get < size_t > ( sandbox_image_key_rows );

    cols_ = props.get < size_t > ( sandbox_image_key_cols );

  } else {
    
    // read rows / cols from the FITS header

    int np;
    int myp;

    MPI_Comm_size ( MPI_COMM_WORLD, &np );
    MPI_Comm_rank ( MPI_COMM_WORLD, &myp );

    fitsfile *fp;
    
    if ( myp == 0 ) {

      fits::open_read ( fp, path_ );

      fits::img_seek ( fp, sighdu_ );
      
      fits::img_dims ( fp, rows_, cols_ );

    }

    int ret = MPI_Bcast ( (void*)(&rows_), 1, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD );
    mpi_check ( MPI_COMM_WORLD, ret );

    ret = MPI_Bcast ( (void*)(&cols_), 1, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD );
    mpi_check ( MPI_COMM_WORLD, ret );

    // read and broadcast the sky flag

    size_t nspec;
    vector < uint8_t > skyflags;

    if ( myp == 0 ) {
      fits::bin_seek ( fp, skyhdu_ );

      long nrows;
      int status = 0;
      ret = fits_get_num_rows ( fp, &nrows, &status );
      fits::check ( status );

      nspec = nrows;
      skyflags.resize( nspec );

      vector < string > colnames ( 1 );
      vector < int > cols ( 1 );
      colnames[0] = "OBJTYPE";
      cols = fits::bin_columns ( fp, colnames );

      vector < string > objnames;
      fits::bin_read_strings ( fp, 0, nspec - 1, cols[0], objnames );

      for ( size_t i = 0; i < nspec; ++i ) {
        if ( objnames[i] == "SKY" ) {
          skyflags[i] = 1;
        } else {
          skyflags[i] = 0;
        }
      }

      fits::close ( fp );
    }

    ret = MPI_Bcast ( (void*)(&nspec), 1, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD );
    mpi_check ( MPI_COMM_WORLD, ret );

    if ( myp != 0 ) {
      skyflags.resize ( nspec );
    }

    ret = MPI_Bcast ( (void*)(&(skyflags[0])), nspec, MPI_CHAR, 0, MPI_COMM_WORLD );
    mpi_check ( MPI_COMM_WORLD, ret );

    sky_.resize ( nspec );

    for ( size_t i = 0; i < nspec; ++i ) {
      if ( skyflags[i] == 0 ) {
        sky_[i] = false;
      } else {
        sky_[i] = true;
      }
    }
    
  }  
  
}


harp::image_sandbox::~image_sandbox ( ) {
  
  cleanup();
  
}


boost::property_tree::ptree harp::image_sandbox::serialize ( ) {
  boost::property_tree::ptree ret;

  ret.put ( "format", image::format() );

  if ( sighdu_ != 1 ) {
    ret.put ( sandbox_image_key_signal, sighdu_ );
  }

  if ( nsehdu_ != 2 ) {
    ret.put ( sandbox_image_key_noise, nsehdu_ );
  }

  if ( path_ == "" ) {
    ret.put ( sandbox_image_key_rows, rows_ );
    ret.put ( sandbox_image_key_cols, cols_ );
  } else {
    ret.put ( sandbox_image_key_path, path_ );
  }

  return ret;
}


void harp::image_sandbox::read ( matrix_local & data ) {

  int np;
  int myp;

  MPI_Comm_size ( MPI_COMM_WORLD, &np );
  MPI_Comm_rank ( MPI_COMM_WORLD, &myp );

  size_t imgsize;
  
  if ( myp == 0 ) {
  
    fitsfile *fp;

    fits::open_read ( fp, path_ );

    fits::img_seek ( fp, sighdu_ );
    
    fits::img_read ( fp, data );
    
    fits::close ( fp );

    imgsize = data.Height();

  }

  int ret = MPI_Bcast ( (void*)(&imgsize), 1, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD );
  mpi_check ( MPI_COMM_WORLD, ret );

  if ( myp != 0 ) {
    data.ResizeTo ( imgsize, 1 );
  }

  ret = MPI_Bcast ( (void*)(data.Buffer()), imgsize, MPI_DOUBLE, 0, MPI_COMM_WORLD );
  mpi_check ( MPI_COMM_WORLD, ret );
  
  return;
}


void harp::image_sandbox::write ( std::string const & path, matrix_local & data ) {
  
  int np;
  int myp;

  MPI_Comm_size ( MPI_COMM_WORLD, &np );
  MPI_Comm_rank ( MPI_COMM_WORLD, &myp );
  
  if ( myp == 0 ) {

    fitsfile *fp;
    
    fits::open_readwrite ( fp, path );
    
    int nh = fits::nhdus ( fp );

    if ( nh < sighdu_ ) {
      while ( nh < sighdu_ ) {
        fits::img_append ( fp, rows_, cols_ );
        ++nh;
      }
    } else {
      fits::img_seek ( fp, sighdu_ );
    }
    
    fits::img_write ( fp, data );
    
    fits::close ( fp );
  }
  
  return;
}


void harp::image_sandbox::read_noise ( matrix_local & data ) {

  int np;
  int myp;

  MPI_Comm_size ( MPI_COMM_WORLD, &np );
  MPI_Comm_rank ( MPI_COMM_WORLD, &myp );

  size_t imgsize;
  
  if ( myp == 0 ) {
  
    fitsfile *fp;

    fits::open_read ( fp, path_ );

    fits::img_seek ( fp, nsehdu_ );
    
    fits::img_read ( fp, data );
    
    fits::close ( fp );

    imgsize = data.Height();
  }

  int ret = MPI_Bcast ( (void*)(&imgsize), 1, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD );
  mpi_check ( MPI_COMM_WORLD, ret );

  if ( myp != 0 ) {
    data.ResizeTo ( imgsize, 1 );
  }

  ret = MPI_Bcast ( (void*)(data.Buffer()), imgsize, MPI_DOUBLE, 0, MPI_COMM_WORLD );
  mpi_check ( MPI_COMM_WORLD, ret );
  
  return;
}


void harp::image_sandbox::write_noise ( std::string const & path, matrix_local & data ) {
  
  int np;
  int myp;

  MPI_Comm_size ( MPI_COMM_WORLD, &np );
  MPI_Comm_rank ( MPI_COMM_WORLD, &myp );
  
  if ( myp == 0 ) {
    fitsfile *fp;
    
    fits::open_readwrite ( fp, path );
    
    int nh = fits::nhdus ( fp );

    if ( nh < nsehdu_ ) {
      while ( nh < nsehdu_ ) {
        fits::img_append ( fp, rows_, cols_ );
        ++nh;
      }
    } else {
      fits::img_seek ( fp, nsehdu_ );
    }
    
    fits::img_write ( fp, data );
    
    fits::close ( fp );
  }
  
  return;
}



