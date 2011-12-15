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

    vec_dense temp1 ( nbins );
    vec_dense temp2 ( nbins );
    vec_dense temp3 ( nbins );

    // eigendecompose the inverse covariance

    boost::numeric::ublas::EigenvalueDecomposition eig ( invcov );

    temp1 = eig.getRealEigenvalues();
    mat_densecol eigvec = eig.getV();

    // sqrt of the eigenvalues

    moat::sf::fast_sqrt ( nbins, &(temp1[0]), &(temp2[0]) );

    // compute the diagonal normalization matrix

    mat_densecol :: iterator2 rowit;
    mat_densecol :: iterator1 colit;

    for ( colit = eigvec.begin1(); colit != eigvec.end1(); ++colit ) {
      for ( rowit = colit.begin(); rowit != colit.end(); ++rowit ) {
        temp3 [ rowit.index1() ] += temp2[ rowit.index1() ] * (*rowit);
      }
    }

    moat::la::multiply_mv < mat_densecol, vec_dense, vec_dense > ( eigvec, temp3, temp1, false, true, true, std::string("") );


    // compute the resolution-convolved spectra

    // find S^-1

    vec_dense :: iterator vit;

    for ( vit = temp1.begin(); vit != temp1.end(); ++vit ) {
      (*vit) = 1.0 / (*vit);
    }

    mat_diag sinv ( nbins, temp1.data() );

    // find D^-1/2

    for ( vit = temp2.begin(); vit != temp2.end(); ++vit ) {
      (*vit) = 1.0 / (*vit);
    }

    mat_diag eigdiag ( nbins, temp2.data() );

    // multiply it all together

    // temp1 = W x z

    moat::la::multiply_mv < mat_densecol, vec_dense, vec_dense > ( eigvec, z, temp1, false, true, true, std::string("") );

    // temp2 = D^-1/2 x ( W x z )

    moat::la::multiply_mv < mat_diag, vec_dense, vec_dense > ( eigdiag, temp1, temp2, false, true, true, std::string("") );

    // temp3 = W^T x ( D^-1/2 x ( W x z ) )

    moat::la::multiply_mv < mat_densecol, vec_dense, vec_dense > ( eigvec, temp2, temp3, true, true, true, std::string("") );

    // f = S^-1 x ( W^T x ( D^-1/2 x ( W x z ) ) )

    moat::la::multiply_mv < mat_diag, vec_dense, vec_dense > ( sinv, temp3, f, false, true, true, std::string("") );

    return;
  }



}

#endif

