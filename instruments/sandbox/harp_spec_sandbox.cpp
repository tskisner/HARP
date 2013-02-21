// @COPYRIGHT@

#include <harp_internal.hpp>


using namespace std;
using namespace harp;

static const char * format_sandbox = "sandbox";

static const char * sandbox_spec_key_path = "path";
static const char * sandbox_spec_key_hdu = "hdu";
static const char * sandbox_spec_key_nspec = "nspec";
static const char * sandbox_spec_key_specsize = "specsize";


harp::spec_sandbox::spec_sandbox ( boost::property_tree::ptree const & props ) : spec ( format_sandbox, props ) {
  
  hdu_ = props.get ( sandbox_spec_key_hdu, 1 );

  path_ = props.get ( sandbox_spec_key_path, "" );

  if ( path_ == "" ) {

    // no path specified- must specify spectra size

    nspec_ = props.get < size_t > ( sandbox_spec_key_nspec );

    specsize_ = props.get < size_t > ( sandbox_spec_key_specsize );

  } else {
    
    // read size from the FITS header

    int np;
    int myp;

    MPI_Comm_size ( MPI_COMM_WORLD, &np );
    MPI_Comm_rank ( MPI_COMM_WORLD, &myp );
    
    if ( myp == 0 ) {
    
      fitsfile *fp;

      fits::open_read ( fp, path_ );

      fits::img_seek ( fp, hdu_ );
      
      size_t rows, cols;
      
      fits::img_dims ( fp, rows, cols );
      
      nspec_ = rows;
      specsize_ = cols;
      
      fits::close ( fp );
    }

    int ret = MPI_Bcast ( (void*)(&nspec_), 1, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD );
    mpi_check ( MPI_COMM_WORLD, ret );

    ret = MPI_Bcast ( (void*)(&specsize_), 1, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD );
    mpi_check ( MPI_COMM_WORLD, ret );
    
  }
  
  size_ = nspec_ * specsize_;
  
}


harp::spec_sandbox::~spec_sandbox ( ) {
  
  cleanup();
  
}


boost::property_tree::ptree harp::spec_sandbox::serialize ( ) {
  boost::property_tree::ptree ret;

  ret.put ( "format", spec::format() );

  if ( hdu_ != 1 ) {
    ret.put ( sandbox_spec_key_hdu, hdu_ );
  }

  if ( path_ == "" ) {
    ret.put ( sandbox_spec_key_nspec, nspec_ );
    ret.put ( sandbox_spec_key_specsize, specsize_ );
  } else {
    ret.put ( sandbox_spec_key_path, path_ );
  }

  return ret;
}


void harp::spec_sandbox::read ( matrix_dist & data ) {

  
  return;
}


void harp::spec_sandbox::write ( string const & path, matrix_dist & data ) {

  
  return;
}



