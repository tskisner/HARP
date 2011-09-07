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

  template < typename M >
  void svd ( M & mat, vec_dense & eigenvals, std::vector < vec_dense > & eigenvecs ) {

    int N = mat.size1();

    vs_dense vec ( eigenvals.size() );
    svd_gen gen;

    ietl::lanczos < M, vs_dense > lanczos ( mat, vec );

    // Creation of an iteration object: 
       
    int max_iter = 10 * N;
    
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

    while ( (*start) < thresh ) {
      ++start;
      --nvec;
    }

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
                << info.error_info(i) << std::endl << std::endl;
    }

    std::vector < double > :: iterator it;

    int cur = 0;
    for ( it = start; it != end; ++it ) {
      std::cout << (*it) << ":  " << (eigenvecs[cur])[0] << " ..." << std::endl;
      ++cur;
    }

    return;
  }


}

#endif
