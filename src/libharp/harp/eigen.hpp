// @COPYRIGHT@

#ifndef HARP_EIGEN_HPP
#define HARP_EIGEN_HPP

#include <boost/random.hpp>

#include <ietl/interface/ublas.h>
#include <ietl/vectorspace.h>
#include <ietl/lanczos.h>
#include <ietl/iteration.h>



namespace harp {

  typedef ietl::wrapper_vectorspace < vec_dense > vs_dense;
  typedef boost::lagged_fibonacci607 svd_gen;

  template < typename M >
  void svd ( M & mat, vec_dense & eigenvals ) {

    int N = mat.size1();

    vs_dense vec ( eigenvals.size() );
    svd_gen gen;

    ietl::lanczos < M, vs_dense > lanczos ( mat, vec );

    // Creation of an iteration object: 
       
    int max_iter = 10*N;
    
    double rel_tol = 500*std::numeric_limits<double>::epsilon();
    
    double abs_tol = std::pow(std::numeric_limits<double>::epsilon(),2./3.);  
    
    std::cout << "Calculation of 20 lowest converged eigenvalues\n\n";
    std::cout << "-----------------------------------\n\n";
    int n_lowest_eigenval = 20;

    std::vector<double> eigen;
    std::vector<double> err;
    std::vector<int> multiplicity;
    ietl::lanczos_iteration_nlowest<double> iter(max_iter, n_lowest_eigenval, rel_tol, abs_tol);
    
    lanczos.calculate_eigenvalues(iter,gen);
    eigen = lanczos.eigenvalues();
    err = lanczos.errors();
    multiplicity = lanczos.multiplicities();
    
     // Printing eigenvalues with error & multiplicities:  
    std::cout << "#        eigenvalue            error         multiplicity\n";
    std::cout.precision(10);
    for (int i=0;i<eigen.size();++i)
      std::cout << i << "\t" << eigen[i] << "\t" << err[i] << "\t" 
                << multiplicity[i] << "\n";
    



    return;
  }


}

#endif
