// @COPYRIGHT@

#ifndef HARP_EXTRACT_HPP
#define HARP_EXTRACT_HPP

#include <boost/numeric/ublas/matrix_proxy.hpp>


namespace harp {

  template < class P, class N, class V, class W >
  void noise_weighted_spec ( boost::numeric::ublas::matrix_expression < P > const & psf, boost::numeric::ublas::matrix_expression < N > const & invnoise, V & img, W & z ) {

    typedef V img_type;
    typedef W z_type;

//  template < class P, class N, class V >
//  void noise_weighted_spec ( P & psf, N & invnoise, V & img, vec_dense & z ) {
    size_t npix = psf().size1();
    size_t nbins = psf().size2();

    if ( invnoise().size1() != npix ) {
      std::ostringstream o;
      o << "number of rows in inverse noise covariance (" << invnoise().size1() << ") does not match number of pixels in PSF (" << npix << ")";
      MOAT_THROW( o.str().c_str() );
    }

    if ( invnoise().size2() != npix ) {
      std::ostringstream o;
      o << "number of columns in inverse noise covariance (" << invnoise().size2() << ") does not match number of pixels in PSF (" << npix << ")";
      MOAT_THROW( o.str().c_str() );
    }

    if ( img.size() != npix ) {
      std::ostringstream o;
      o << "number of elements in image vector (" << img.size() << ") does not match number of pixels in PSF (" << npix << ")";
      MOAT_THROW( o.str().c_str() );
    }

    z.resize ( nbins );
    z_type temp ( npix );
    
    boost::numeric::ublas::axpy_prod ( invnoise, img, temp, true );

    boost::numeric::ublas::axpy_prod ( boost::numeric::ublas::trans ( psf ), temp, z, true );
    

    //moat::la::multiply_mv < N, V, vec_dense > ( invnoise, img, temp, false, true, true, std::string("") );

    //moat::la::multiply_mv < P, vec_dense, vec_dense > ( psf, temp, z, true, true, true, std::string("") );

    return;
  }


  template < class P, class N, class R >
  void inverse_covariance ( boost::numeric::ublas::matrix_expression < P > const & psf, boost::numeric::ublas::matrix_expression < N > const & invnoise, boost::numeric::ublas::matrix_expression < R > & output ) {

    // check consistent sizes

    if ( invnoise().size1() != invnoise().size2() ) {
      MOAT_THROW( "inverse noise covariance must be square" );
    }

    if ( output().size1() != output().size2() ) {
      MOAT_THROW( "output inverse spectral covariance must be square" );
    }

    if ( invnoise().size1() != psf().size1() ) {
      std::ostringstream o;
      o << "number of rows in inverse noise covariance (" << invnoise().size1() << ") does not match number of pixels in PSF (" << psf().size1() << ")";
      MOAT_THROW( o.str().c_str() );
    }

    if ( output().size1() != psf().size2() ) {
      std::ostringstream o;
      o << "dimension of output spectral covariance (" << output().size2() << ") does not match number of bins in PSF (" << psf().size2() << ")";
      MOAT_THROW( o.str().c_str() );
    }

    size_t nbins = psf().size2();
    size_t npix = psf().size1();

    // construct the inverse spectral covariance matrix

    mat_dynrow builder ( nbins, nbins );

    mat_dynrow temp ( npix, nbins );

    std::cerr << "  multiply N^-1 A..." << std::endl;

    boost::numeric::ublas::axpy_prod ( invnoise, psf, temp, true );

    std::cerr << "  multiply A^T N^-1..." << std::endl;
    
    boost::numeric::ublas::axpy_prod ( boost::numeric::ublas::trans ( psf ), temp, builder, true );

    output().assign ( builder );

    return;
  }



  template < class C, class V, class W >
  void extract_dense ( boost::numeric::ublas::matrix_expression < C > const & invcov, V & z, W & f ) {

    typedef V z_type;
    typedef W f_type;

    size_t nbins = invcov().size1();

    if ( invcov().size2() != nbins ) {
      std::ostringstream o;
      o << "inverse spectral covariance must be square";
      MOAT_THROW( o.str().c_str() );
    }

    if ( z.size() != nbins ) {
      std::ostringstream o;
      o << "number of elements in noise weighted image vector (" << z.size() << ") does not match number of flux bins in spectral covariance (" << nbins << ")";
      MOAT_THROW( o.str().c_str() );
    }
    
    f.resize ( nbins );

    vec_dense sqrteig ( nbins );
    vec_dense temp1 ( nbins );

    mat_denserow tempmat1 ( nbins, nbins );
    mat_denserow resmat ( nbins, nbins );

    // eigendecompose the inverse covariance.  Note that the *rows*
    // of this matrix are the eigenvectors.

    boost::numeric::ublas::EigenvalueDecomposition eig ( invcov );

    mat_denserow eigval = eig.getD();
    mat_denserow eigvec = eig.getV();

    // sqrt of the eigenvalues

    //moat::sf::fast_sqrt ( nbins, &(eigval[0]), &(sqrteig[0]) );

    // find D^1/2

    mat_densecol eigdiag ( nbins, nbins );

    vec_dense :: iterator vit;

    // DEBUG, test eigen-recomposition

    /*
    for ( size_t i = 0; i < nbins; ++i ) {
      eigdiag( i, i ) = eigval(i,i);
    }

    boost::numeric::ublas::axpy_prod ( eigdiag, boost::numeric::ublas::trans ( eigvec ), tempmat1, true );

    boost::numeric::ublas::axpy_prod ( eigvec, tempmat1, resmat, true );
    
    //moat::la::multiply_mm < mat_denserow, mat_diag, mat_densecol > ( eigvec, eigdiag, tempmat1, true, true, true, std::string("") );

    //moat::la::multiply_mm < mat_densecol, mat_denserow, mat_densecol > ( tempmat1, eigvec, resmat, false, true, true, std::string("") );

    mat_denserow :: iterator1 rowit;
    mat_denserow :: iterator2 colit;

    double checkval;
    for ( rowit = resmat.begin1(); rowit != resmat.end1(); ++rowit ) {
      for ( colit = rowit.begin(); colit != rowit.end(); ++colit ) {
        checkval = invcov ( colit.index1(), colit.index2() );
        if ( ( fabs ( checkval ) > moat::EPSILON_SINGLE ) && fabs ( (*colit) - checkval ) / checkval > 1.0e-10 ) {
          std::cerr << "recompose error on (" << colit.index1() << ", " << colit.index2() << ") " << (*colit) << " != " << checkval << std::endl;
        }
      }
    }
    */

    for ( size_t i = 0; i < nbins; ++i ) {
      eigdiag( i, i ) = sqrt ( eigval(i,i) );
    }


    // compute the resolution matrix

    boost::numeric::ublas::axpy_prod ( eigdiag, boost::numeric::ublas::trans ( eigvec ), tempmat1, true );

    boost::numeric::ublas::axpy_prod ( eigvec, tempmat1, resmat, true );

    //moat::la::multiply_mm < mat_densecol, mat_densecol, mat_densecol > ( eigdiag, boost::numeric::ublas::trans ( eigvec ), tempmat1, false, true, true, std::string("") );

    //moat::la::multiply_mm < mat_densecol, mat_densecol, mat_densecol > ( eigvec, tempmat1, resmat, false, true, true, std::string("") );

    mat_denserow :: iterator1 rowit;
    mat_denserow :: iterator2 colit;

    temp1.clear();
    for ( rowit = resmat.begin1(); rowit != resmat.end1(); ++rowit ) {
      for ( colit = rowit.begin(); colit != rowit.end(); ++colit ) {
        temp1 [ colit.index1() ] += (*colit);
      }
    }

    for ( rowit = resmat.begin1(); rowit != resmat.end1(); ++rowit ) {
      for ( colit = rowit.begin(); colit != rowit.end(); ++colit ) {
        (*colit) /= temp1 [ colit.index1() ];
      }
    }

    // compute covariance

    mat_densecol covmat ( nbins, nbins );

    for ( size_t i = 0; i < nbins; ++i ) {
      eigdiag( i, i ) = 1.0 / eigval(i,i);
    }

    boost::numeric::ublas::axpy_prod ( eigdiag, boost::numeric::ublas::trans ( eigvec ), tempmat1, true );

    boost::numeric::ublas::axpy_prod ( eigvec, tempmat1, covmat, true );

    //moat::la::multiply_mm < mat_densecol, mat_densecol, mat_densecol > ( eigdiag, boost::numeric::ublas::trans ( eigvec ), tempmat1, false, true, true, std::string("") );

    //moat::la::multiply_mm < mat_densecol, mat_densecol, mat_densecol > ( eigvec, tempmat1, covmat, false, true, true, std::string("") );

    // compute the resolution-convolved spectra
    // Rf = RCZ

    boost::numeric::ublas::axpy_prod ( covmat, z, temp1, true );

    boost::numeric::ublas::axpy_prod ( resmat, temp1, f, true );

    //moat::la::multiply_mv < mat_densecol, vec_dense, vec_dense > ( covmat, z, temp1, false, true, true, std::string("") );

    //moat::la::multiply_mv < mat_densecol, vec_dense, vec_dense > ( resmat, temp1, f, false, true, true, std::string("") );
    

    return;
  }


  template < class P, class S, class R >
  void append_sky_comp ( boost::numeric::ublas::matrix_expression < P > const & psf, boost::numeric::ublas::matrix_expression < S > const & sky, size_t nspec, boost::numeric::ublas::matrix_expression < R > & output ) {

    typedef boost::numeric::ublas::matrix_range < R > R_view;


    // verify that length of sky components equals number of columns in design matrix

    size_t npix = psf().size1();
    size_t nbins = psf().size2();
    size_t ncomp = sky().size1();
    size_t nspecbins = (size_t) ( nbins / nspec );

    if ( sky().size2() != nspecbins ) {
      std::ostringstream o;
      o << "number of columns in sky component matrix (" << sky().size2() << ") does not match number of bins per spectrum in PSF (" << nspecbins << ")";
      MOAT_THROW( o.str().c_str() );
    }

    if ( ncomp > nspecbins ) {
      std::ostringstream o;
      o << "number of sky components (" << ncomp << ") exceeds the number of flux bins per spectrum (" << nbins << ")";
      MOAT_THROW( o.str().c_str() );
    }

    // resize output matrix and copy the original design matrix into the
    // first column block

    output().resize ( npix, nbins + ncomp * nspecbins );

    R_view original ( output(), mv_range( 0, npix ), mv_range( 0, nbins ) );

    original.assign ( psf() );

    // construct component matrix

    mat_compcol compmat ( nbins, ncomp * nspec, ncomp * nbins );

    for ( size_t i = 0; i < nspec; ++i ) {
      for ( size_t j = 0; j < ncomp; ++j ) {
        for ( size_t k = 0; k < nspecbins; ++k ) {
          compmat ( i * nspecbins + k, i * ncomp + j ) = sky() ( j, k );
        }
      }
    }

    // fill in sky portion of output design matrix

    R_view skyblock ( output(), mv_range ( 0, npix ), mv_range ( nbins, nbins + ncomp * nspecbins ) );

    boost::numeric::ublas::axpy_prod ( psf(), compmat, skyblock, true );

    return;
  }




}

#endif

