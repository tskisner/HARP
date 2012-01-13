// @COPYRIGHT@

#ifndef HARP_EXTRACT_HPP
#define HARP_EXTRACT_HPP


namespace harp {

  template < class P, class N, class V >
  void noise_weighted_spec ( P & psf, N & invnoise, V & img, vec_dense & z ) {
    size_t npix = psf.size1();
    size_t nbins = psf.size2();

    if ( invnoise.size1() != npix ) {
      std::ostringstream o;
      o << "number of rows in inverse noise covariance (" << invnoise.size1() << ") does not match number of pixels in PSF (" << npix << ")";
      MOAT_THROW( o.str().c_str() );
    }

    if ( invnoise.size2() != npix ) {
      std::ostringstream o;
      o << "number of columns in inverse noise covariance (" << invnoise.size2() << ") does not match number of pixels in PSF (" << npix << ")";
      MOAT_THROW( o.str().c_str() );
    }

    if ( img.size() != npix ) {
      std::ostringstream o;
      o << "number of elements in image vector (" << img.size() << ") does not match number of pixels in PSF (" << npix << ")";
      MOAT_THROW( o.str().c_str() );
    }

    z.resize ( nbins );
    vec_dense temp ( npix );

    moat::la::multiply_mv < N, V, vec_dense > ( invnoise, img, temp, false, true, true, std::string("") );

    moat::la::multiply_mv < P, vec_dense, vec_dense > ( psf, temp, z, true, true, true, std::string("") );

    return;
  }


  template < class C >
  void extract_dense ( C & invcov, vec_dense & z, vec_dense & f ) {

    size_t nbins = invcov.size1();

    if ( invcov.size2() != nbins ) {
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

    mat_densecol tempmat1 ( nbins, nbins );
    mat_densecol resmat ( nbins, nbins );

    // eigendecompose the inverse covariance

    boost::numeric::ublas::EigenvalueDecomposition eig ( invcov );

    mat_denserow eigval = eig.getD();
    mat_denserow eigvec = eig.getV();

    // sqrt of the eigenvalues

    //moat::sf::fast_sqrt ( nbins, &(eigval[0]), &(sqrteig[0]) );

    // find D^1/2

    mat_diag eigdiag ( nbins, nbins );

    vec_dense :: iterator vit;

    // DEBUG, test eigen-recomposition

    for ( size_t i = 0; i < nbins; ++i ) {
      eigdiag( i, i ) = eigval(i,i);
    }

    moat::la::multiply_mm < mat_denserow, mat_diag, mat_densecol > ( eigvec, eigdiag, tempmat1, true, true, true, std::string("") );

    moat::la::multiply_mm < mat_densecol, mat_denserow, mat_densecol > ( tempmat1, eigvec, resmat, false, true, true, std::string("") );

    mat_densecol :: iterator2 rowit;
    mat_densecol :: iterator1 colit;

    double checkval;
    for ( colit = resmat.begin1(); colit != resmat.end1(); ++colit ) {
      for ( rowit = colit.begin(); rowit != colit.end(); ++rowit ) {
        checkval = invcov ( rowit.index1(), rowit.index2() );
        if ( fabs ( (*rowit) - checkval ) / checkval > 1.0e-6 ) {
          std::cerr << "recompose error on (" << rowit.index1() << ", " << rowit.index2() << ") " << (*rowit) << " != " << checkval << std::endl;
        }
      }
    }

    for ( size_t i = 0; i < nbins; ++i ) {
      eigdiag( i, i ) = sqrt ( eigval(i,i) );
    }


    // compute the resolution matrix

    moat::la::multiply_mm < mat_densecol, mat_diag, mat_densecol > ( eigvec, eigdiag, tempmat1, true, true, true, std::string("") );

    moat::la::multiply_mm < mat_densecol, mat_densecol, mat_densecol > ( tempmat1, eigvec, resmat, false, true, true, std::string("") );

    temp1.clear();
    for ( colit = resmat.begin1(); colit != resmat.end1(); ++colit ) {
      for ( rowit = colit.begin(); rowit != colit.end(); ++rowit ) {
        temp1 [ rowit.index2() ] += (*rowit);
      }
    }

    for ( colit = resmat.begin1(); colit != resmat.end1(); ++colit ) {
      for ( rowit = colit.begin(); rowit != colit.end(); ++rowit ) {
        (*rowit) /= temp1 [ rowit.index2() ];
      }
    }

    // compute covariance

    mat_densecol covmat ( nbins, nbins );

    for ( size_t i = 0; i < nbins; ++i ) {
      eigdiag( i, i ) = 1.0 / eigdiag(i,i);
    }

    moat::la::multiply_mm < mat_densecol, mat_diag, mat_densecol > ( eigvec, eigdiag, tempmat1, true, true, true, std::string("") );

    moat::la::multiply_mm < mat_densecol, mat_densecol, mat_densecol > ( tempmat1, eigvec, covmat, false, true, true, std::string("") );

    // compute the resolution-convolved spectra
    // Rf = RCZ

    moat::la::multiply_mv < mat_densecol, vec_dense, vec_dense > ( covmat, z, temp1, false, true, true, std::string("") );

    moat::la::multiply_mv < mat_densecol, vec_dense, vec_dense > ( resmat, temp1, f, false, true, true, std::string("") );
    

    return;
  }



}

#endif

