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
      
    }


    vs_dense vs ( dim );
    svd_gen gen;

    ietl::lanczos < M, vs_dense > lanczos ( invcov, vs );

    // create iterator

    int max_iter = 2 * dim;
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

    std::vector < double > :: reverse_iterator firsteigen = eigen.rbegin();
    std::vector < double > :: reverse_iterator lasteigen = firsteigen;

    size_t nonsing = 0;

    while ( ( lasteigen != eigen.rend() ) && ( (*lasteigen) > thresh ) ) {
      ++lasteigen;
      ++nonsing;
    }

    /*
    ietl::Info < double > info; // (m1, m2, ma, eigenvalue, residual, status).
    
    lanczos.eigenvectors ( start, end, std::back_inserter ( eigenvecs ), info, gen ); 
    
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


  template < typename M >
  void svd ( M & mat, vec_dense & eigenvals, std::vector < vec_dense > & eigenvecs ) {

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

    eigenvecs.clear();
    //eigenvecs.resize ( nvec );
 
    ietl::Info < double > info; // (m1, m2, ma, eigenvalue, residual, status).
    
    lanczos.eigenvectors ( start, end, std::back_inserter ( eigenvecs ), info, gen ); 
    
    std::cout << "Printing eigen Vectors:" << std::endl << std::endl; 
    
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

    return;
  }


}

#endif
