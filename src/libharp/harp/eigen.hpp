// @COPYRIGHT@

#ifndef HARP_EIGEN_HPP
#define HARP_EIGEN_HPP

#include <climits>

#include <boost/random.hpp>

#include <ietl/interface/ublas.h>
#include <ietl/vectorspace.h>
#include <ietl/lanczos.h>
#include <ietl/iteration.h>



namespace harp {

  int lapack_ev ( int dim, double * mat, double * eigen );

  typedef ietl::wrapper_vectorspace < vec_dense > vs_dense;
  typedef boost::lagged_fibonacci607 svd_gen;

  template < typename M, typename V, typename R >
  void extract ( M & inout, V & z, V & spectra, V & errors, R & resolution, std::vector < size_t > & Rblocks ) {

    size_t dim = inout.size1();

    if ( inout.size2() != dim ) {
      MOAT_THROW( "input/output matrix must be square" );
    }

    if ( ( resolution.size1() != dim ) || ( resolution.size2() != dim ) ) {
      MOAT_THROW( "resolution matrix must be the same size as the input matrix" );
    }

    if ( z.size() != dim ) {
      MOAT_THROW( "input Z vector size must match number of columns in the input matrix" );
    }

    spectra.clear();
    spectra.resize ( dim );

    errors.clear();
    errors.resize ( dim );

    size_t blocktotal = 0;
    std::vector < size_t > :: iterator blit;
    for ( blit = Rblocks.begin(); blit != Rblocks.end(); ++blit ) {
      blocktotal += (*blit);
    }

    if ( blocktotal != dim ) {
      MOAT_THROW( "sum of all resolution blocksizes must match matrix dimensions" );
    }

    V eigenvals ( dim );


    // Dense LAPACK implementation


    size_t i, j, k;

    // eigendecompose


    
    vec_dense finvcov ( dim * dim );
    mat_densecol W ( dim, dim );

    for ( j = 0; j < dim; ++j ) {
      for ( i = 0; i < dim; ++i ) {
        finvcov[ j * dim + i ] = inout ( i, j );
      }
    }

    int fret = lapack_ev ( (int)dim, &(finvcov(0)), &(eigenvals[0]) );

    for ( j = 0; j < dim; ++j ) {
      for ( i = 0; i < dim; ++i ) {
        W ( j, i ) = finvcov[ j * dim + i ];
      }
    }

    size_t neigen = (size_t)( dim );
 
    std::cerr << "Condition number = " << eigenvals[neigen-1] / eigenvals[0] << std::endl;

    


    // compute norm

    std::cerr << "compute norm" << std::endl;

    vec_dense norm ( dim );

    vec_dense sqrtvals ( dim );

    moat::sf::sqrt ( neigen, &(eigenvals[0]), &(sqrtvals[0]) );

    mat_comprow diag ( dim, dim, dim );
    for ( i = 0; i < neigen; ++i ) {
      diag ( i, i ) = sqrtvals[i];
    }
    for ( i = neigen; i < dim; ++i ) {
      diag ( i, i ) = 0.0;
    }

    vec_dense vtemp1 ( dim );
    vec_dense vtemp2 ( dim );
    
    mat_denserow mtemp1 ( dim, dim );

    vec_dense rc ( dim );
    vec_dense ones ( dim );

    for ( i = 0; i < dim; ++i ) {
      ones[i] = 1.0;
    }

    boost::numeric::ublas::axpy_prod ( W, ones, vtemp1, true );

    boost::numeric::ublas::axpy_prod ( diag, vtemp1, vtemp2, true );

    boost::numeric::ublas::axpy_prod ( vtemp2, W, norm, true );

    // compute R explicitly

    std::cerr << "compute R" << std::endl;

    boost::numeric::ublas::axpy_prod ( diag, W, mtemp1, true );

    boost::numeric::ublas::axpy_prod ( boost::numeric::ublas::trans ( W ), mtemp1, resolution, true );

    for ( i = 0; i < dim; ++i ) {
      for ( j = 0; j < dim; ++j ) {
        resolution ( i, j ) /= norm[i];
      }
    }


    // compute C explicity

    std::cerr << "compute C" << std::endl;

    mat_denserow cov ( dim, dim );

    for ( i = 0; i < neigen; ++i ) {
      diag ( i, i ) = 1.0 / eigenvals[i];
    }
    for ( i = neigen; i < dim; ++i ) {
      diag ( i, i ) = 0.0;
    }

    boost::numeric::ublas::axpy_prod ( diag, W, mtemp1, true );

    boost::numeric::ublas::axpy_prod ( boost::numeric::ublas::trans ( W ), mtemp1, cov, true );

    boost::numeric::ublas::axpy_prod ( resolution, cov, mtemp1, true );

    boost::numeric::ublas::axpy_prod ( mtemp1, z, spectra, true );

    for ( size_t i = 0; i < dim; ++i ) {
      errors[i] = 1.0 / norm[i];
    }


    /*
    // IETL Lanczos specific code

    int idim = (int)dim;
    vs_dense vs ( idim );
    svd_gen gen;

    ietl::lanczos < M, vs_dense > lanczos ( inout, vs );

    // create iterator

    int max_iter = 3 * dim;
    ietl::fixed_lanczos_iteration < double > iter ( max_iter );

    // compute eigenvalues

    std::vector < double > eigen;
    std::vector < double > err;
    std::vector < int > multiplicity;  
    
    lanczos.calculate_eigenvalues ( iter, gen );

    eigen = lanczos.eigenvalues();
    err = lanczos.errors();
    multiplicity = lanczos.multiplicities();

    // determine the truncation threshold

    double largest = eigen[ eigen.size() - 1 ];
    double thresh = largest * std::numeric_limits < double > :: epsilon();

    // set iterators from largest eigenvalue to smallest

    std::vector < double > :: iterator firsteigen = eigen.begin();
    std::vector < double > :: iterator lasteigen = eigen.end();

    size_t nonsing = dim;

    while ( (*firsteigen) < thresh ) {
      ++firsteigen;
      --nonsing;
    }

    // For now, get ALL eigenvectors

    ietl::Info < double > info; // (m1, m2, ma, eigenvalue, residual, status).

    std::vector < vec_dense > eigenvecs;
    
    lanczos.eigenvectors ( firsteigen, lasteigen, std::back_inserter ( eigenvecs ), info, gen ); 
    
    std::cout << "Printing eigen Vectors:" << std::endl << std::endl; 
    
    std::cout << " Information about the eigen vectors computations:\n\n";
    for ( int i = 0; i < info.size(); i++ ) {
      std::cout << " m1(" << i+1 << "): " << info.m1(i) << ", m2(" << i+1 << "): "
                << info.m2(i) << ", ma(" << i+1 << "): " << info.ma(i) << " eigenvalue("
                << i+1 << "): " << info.eigenvalue(i) << " residual(" << i+1 << "): "
                << info.residual(i) << " error_info(" << i+1 << "): "
                << info.error_info(i) << std::endl << std::endl;
    }

    std::vector < double > :: iterator it;

    int cur = 0;
    for ( it = start; it != end; ++it ) {
      std::cout << (*it) << ":  " << (eigenvecs[cur])[0] << " ..." << std::endl;
      ++cur;
    }

    */

    return;
  }


  template < typename M, typename E >
  void svd ( M & mat, vec_dense & eigenvals, E & W ) {

    int N = eigenvals.size();

    std::cerr << "SVD using matrix size " << N << std::endl;

    vs_dense vec ( eigenvals.size() );
    svd_gen gen;

    std::cerr << "constructing lanczos" << std::endl;

    ietl::lanczos < M, vs_dense > lanczos ( mat, vec );

    // Creation of an iteration object: 
       
    int max_iter = 2 * N;
    
    std::cerr << "fixed iterations" << std::endl;

    std::vector < double > eigen;
    std::vector < double > err;
    std::vector < int > multiplicity;
    ietl::fixed_lanczos_iteration < double > iter ( max_iter );
    
    lanczos.calculate_eigenvalues ( iter, gen );
    eigen = lanczos.eigenvalues();
    err = lanczos.errors();
    multiplicity = lanczos.multiplicities();
    
    // Printing eigenvalues with error & multiplicities:  
    std::cout << "#        eigenvalue            error         multiplicity\n";
    std::cout.precision(10);
    
    for ( int i = 0; i < eigen.size(); ++i ) {
      eigenvals[i] = eigen[i];
      std::cout << i << "\t" << eigen[i] << "\t" << err[i] << "\t" << multiplicity[i] << "\n";
    }

    double largest = eigenvals[ eigenvals.size() - 1 ];
    double thresh = largest * std::numeric_limits < double > :: epsilon();

    std::vector < double > :: iterator start = eigen.begin();
    std::vector < double > :: iterator end = eigen.end();

    size_t nvec = eigenvals.size();

    /*
    while ( (*start) < thresh ) {
      ++start;
      --nvec;
    }
    */

    ietl::Info < double > info; // (m1, m2, ma, eigenvalue, residual, status).

    std::vector < vec_dense > eigenvecs;
    
    lanczos.eigenvectors ( start, end, std::back_inserter ( eigenvecs ), info, gen ); 
    
    std::cout << "Printing eigen Vectors:" << std::endl << std::endl;
    
    /*
    std::cout << " Information about the eigen vectors computations:\n\n";
    for ( int i = 0; i < info.size(); i++ ) {
      std::cout << " m1(" << i+1 << "): " << info.m1(i) << ", m2(" << i+1 << "): "
                << info.m2(i) << ", ma(" << i+1 << "): " << info.ma(i) << " eigenvalue("
                << i+1 << "): " << info.eigenvalue(i) << " residual(" << i+1 << "): "
                << info.residual(i) << " error_info(" << i+1 << "): "
                << info.error_info(i) << std::endl;
    }

    std::vector < double > :: iterator it;

    int cur = 0;
    for ( it = start; it != end; ++it ) {

      std::cout << (*it) << ":  " << (eigenvecs[cur])[0] << " " << (eigenvecs[cur])[1] << " " << (eigenvecs[cur])[2] << " ..." << std::endl;
      ++cur;
    }

    */

    int i;
    int j;
    for ( int j = 0; j < N; ++j ) {
      for ( i = 0; i < N; ++i ) {
        W ( i, j ) = (eigenvecs[j])[i];
      }
    }


    return;
  }


}

#endif
