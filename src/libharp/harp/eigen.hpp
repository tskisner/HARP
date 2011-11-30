// @COPYRIGHT@

#ifndef HARP_EIGEN_HPP
#define HARP_EIGEN_HPP

#include <iostream>
#include <fstream>
#include <climits>

#include <boost/random.hpp>

#include <moat.hpp>

#include <ietl/interface/ublas.h>
#include <ietl/vectorspace.h>
#include <ietl/lanczos.h>
#include <ietl/iteration.h>


namespace harp {

  int lapack_ev ( int dim, double * mat, double * eigen );

  typedef ietl::wrapper_vectorspace < vec_dense > vs_dense;
  typedef boost::lagged_fibonacci607 svd_gen;

  template < typename M, typename V >
  void extract ( M & inout, V & z, V & spectra, V & errors, std::vector < V > & intestspec, std::vector < V > & outtestspec ) {

    size_t dim = inout.size1();

    if ( inout.size2() != dim ) {
      MOAT_THROW( "input/output matrix must be square" );
    }

    if ( z.size() != dim ) {
      MOAT_THROW( "input Z vector size must match number of columns in the input matrix" );
    }

    spectra.clear();
    spectra.resize ( dim );

    errors.clear();
    errors.resize ( dim );

    size_t ntest = intestspec.size();
    outtestspec.clear();

    size_t i, j, k, m;

    for ( i = 0; i < ntest; ++i ) {
      if ( intestspec[i].size() != dim ) {
        MOAT_THROW( "input test spectra must have same size as dimensions of the inverse covariance ");
      }
      V tempvec ( dim );
      outtestspec.push_back ( tempvec );
    }

    //V eigenvals ( dim );
    size_t neigen = dim;


    /*
    // Dense LAPACK implementation

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
    */
    
    // sparse eigen decomposition


    int idim = (int)dim;
    vs_dense vs ( idim );
    svd_gen gen;

    ietl::lanczos < M, vs_dense > lanczos ( inout, vs );

    // compute eigenvalues

    std::vector < double > eigen;
    std::vector < double > err;
    std::vector < int > multiplicity;

    // create iterator

    int max_iter = 2 * dim;

    double rel_tol = 100 * std::numeric_limits < double > :: epsilon();
    double abs_tol = rel_tol; //std::pow ( std::numeric_limits < double > :: epsilon(), 2.0 / 3.0 );

    ietl::lanczos_iteration_nhighest < double > iter ( max_iter, dim, rel_tol, abs_tol );

    //ietl::fixed_lanczos_iteration < double > iter ( max_iter );

    lanczos.calculate_eigenvalues ( iter, gen );

    std::cout << "relative error = " << iter.relative_tolerance() << std::endl;
    std::cout << "absolute error = " << iter.absolute_tolerance() << std::endl;
    std::cout << "number of iterations: " << iter.iterations() << std::endl;

    eigen = lanczos.eigenvalues();
    err = lanczos.errors();
    multiplicity = lanczos.multiplicities();

    neigen = eigen.size();

    // Printing eigenvalues with error & multiplicities:  
    std::cout << "#        eigenvalue            error         multiplicity\n";
    std::cout.precision(10);
    for (int i=0;i<eigen.size();++i)
      std::cout << i << "  " << eigen[i] << "  " << err[i] << "  " << multiplicity[i] << "\n";

    // determine the truncation threshold

    double largest = eigen[ neigen - 1 ];
    double thresh = largest * std::numeric_limits < double > :: epsilon();
    std::cerr << "threshold = " << thresh << std::endl;

    // set iterator bounds

    std::vector < double > :: reverse_iterator firsteigen = eigen.rbegin();
    std::vector < double > :: reverse_iterator lasteigen = eigen.rend();

    /*
    size_t nonsing = eigen.size();

    while ( (*firsteigen) < thresh ) {
      ++firsteigen;
      --nonsing;
    }
    */

    // scroll through eigenvectors

    ietl::Info < double > info; // (m1, m2, ma, eigenvalue, residual, status).

    std::vector < vec_dense > eigenvecs;

    std::cerr << "fetching eigenvectors..." << std::endl;
    
    while ( firsteigen != lasteigen ) {

      lanczos.eigenvectors ( firsteigen, firsteigen-1, std::back_inserter ( eigenvecs ), info, gen );

      for ( int i = 0; i < info.size(); i++ ) {
        std::cout << " m1(" << i+1 << "): " << info.m1(i) << ", m2(" << i+1 << "): "
                << info.m2(i) << ", ma(" << i+1 << "): " << info.ma(i) << " eigenvalue("
                << i+1 << "): " << info.eigenvalue(i) << " residual(" << i+1 << "): "
                << info.residual(i) << " error_info(" << i+1 << "): "
                << info.error_info(i) << std::endl;
      }

      ++firsteigen;
    }



    vec_dense sqrtvals ( neigen );

    moat::sf::sqrt ( neigen, &(eigen[0]), &(sqrtvals[0]) );

    
 
    std::cerr << "Condition number = " << eigen[neigen-1] / eigen[0] << std::endl;



    // Using one eigenvector at a time, perform computation of normalization vector,
    // application of resolution matrix on test spectra, and solution of final spectra.
    
    V norm ( dim );
    V norm_temp1 ( dim );
    V norm_temp2 ( dim );

    V rc ( dim );
    V rc_temp1 ( dim );
    V rc_temp2 ( dim );

    std::vector < V > test_temp1 ( ntest );
    std::vector < V > test_temp2 ( ntest );

    for ( i = 0; i < ntest; ++i ) {
      test_temp1[i].resize ( dim );
      test_temp2[i].resize ( dim );
    }

    double * evec;

    // eigenvectors row-major

    size_t cur;

    for ( cur = 0; cur < dim; ++cur ) {
      norm_temp1[cur] = 0.0;
      rc_temp1[cur] = 0.0;
      for ( k = 0; k < ntest; ++k ) {
        ( test_temp1[k] )[cur] = 0.0;
      }
    }
    cur = 0;

    for ( i = 0; i < neigen; ++i ) {
      
      //evec = &( finvcov [ i * dim ] );
      evec = &( (eigenvecs[i])[0] );

      for ( m = 0; m < (size_t)multiplicity[i]; ++m ) {

        for ( j = 0; j < dim; ++j ) {

          norm_temp1[cur] += evec[j];

          rc_temp1[cur] += evec[j] * z[j];

          for ( k = 0; k < ntest; ++k ) {
            ( test_temp1[k] )[cur] += evec[j] * ( intestspec[k] )[j];
          }

        }

        ++cur;
      }
    }

    cur = 0;

    for ( i = 0; i < neigen; ++i ) { 

      for ( m = 0; m < (size_t)multiplicity[i]; ++m ) {

        norm_temp2[cur] = norm_temp1[cur] * sqrtvals[i];
      
        rc_temp2[cur] = rc_temp1[cur] / sqrtvals[i];
      
        for ( k = 0; k < ntest; ++k ) {
          ( test_temp2[k] )[cur] = ( test_temp1[k] )[cur] * sqrtvals[i];
          ( outtestspec[k] )[cur] = 0.0;
        }

        norm[cur] = 0.0;
        rc[cur] = 0.0;

        ++cur;

      }

    }

    cur = 0;

    for ( i = 0; i < neigen; ++i ) {
      
      //evec = &( finvcov [ i * dim ] );
      evec = &( (eigenvecs[i])[0] );

      for ( m = 0; m < (size_t)multiplicity[i]; ++m ) {

        for ( j = 0; j < dim; ++j ) {

          norm[cur] += evec[j] * norm_temp2[cur];

          rc[cur] += evec[j] * rc_temp2[cur];

          for ( k = 0; k < ntest; ++k ) {
            ( outtestspec[k] )[cur] += evec[j] * ( test_temp2[k] )[cur];
          }
        }

        ++cur;

      }

    }


    /*
    // compute norm

    std::cerr << "compute norm" << std::endl;

    vec_dense testnorm ( dim );
    

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

    vec_dense ones ( dim );

    vec_dense testspectra ( dim );

    for ( i = 0; i < dim; ++i ) {
      ones[i] = 1.0;
    }

    boost::numeric::ublas::axpy_prod ( W, ones, vtemp1, true );

    boost::numeric::ublas::axpy_prod ( diag, vtemp1, vtemp2, true );

    boost::numeric::ublas::axpy_prod ( vtemp2, W, testnorm, true );



    // compute R explicitly

    std::cerr << "compute R" << std::endl;

    mat_denserow resolution ( dim, dim );

    boost::numeric::ublas::axpy_prod ( diag, W, mtemp1, true );

    boost::numeric::ublas::axpy_prod ( boost::numeric::ublas::trans ( W ), mtemp1, resolution, true );

    for ( i = 0; i < dim; ++i ) {
      for ( j = 0; j < dim; ++j ) {
        resolution ( i, j ) /= testnorm[i];
      }
    }


    // apply resolution to test spectra

    for ( i = 0; i < ntest; ++i ) {
      boost::numeric::ublas::axpy_prod ( resolution, intestspec[i], outtestspec[i], true );
    }

    // compute C explicity

    std::cerr << "compute RC" << std::endl;

    mat_denserow testrc ( dim, dim );

    for ( i = 0; i < neigen; ++i ) {
      diag ( i, i ) = 1.0 / sqrtvals[i];
    }
    for ( i = neigen; i < dim; ++i ) {
      diag ( i, i ) = 0.0;
    }

    boost::numeric::ublas::axpy_prod ( diag, W, mtemp1, true );

    boost::numeric::ublas::axpy_prod ( boost::numeric::ublas::trans ( W ), mtemp1, testrc, true );


    vec_dense unnorm ( dim );
    boost::numeric::ublas::axpy_prod ( testrc, z, unnorm, true );


    std::string outdata = "testdata/test_medium_debug_norm.out";

    std::fstream out;
    out.open ( outdata.c_str(), std::ios::out );
  
    for ( i = 0; i < dim; ++i ) {
      out << i << " " << rc[i] << " " << unnorm[i] << std::endl;
    }
  
    out.close();

    for ( i = 0; i < dim; ++i ) {
      for ( j = 0; j < dim; ++j ) {
        testrc ( i, j ) /= testnorm[i];
      }
    }

    boost::numeric::ublas::axpy_prod ( testrc, z, testspectra, true );
    */

    /*
    for ( size_t i = 0; i < dim; ++i ) {
      norm[i] = 1.0 / norm[i];
      spectra[i] = rc[i] * norm[i];
      errors[i] = norm[i];
    }

    for ( k = 0; k < ntest; ++k ) {
      for ( j = 0; j < dim; ++j ) {
        ( outtestspec[k] )[j] *= norm[j];
      }
    }
    */

    /*
    // IETL Lanczos specific code

    

     
    
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


}

#endif
