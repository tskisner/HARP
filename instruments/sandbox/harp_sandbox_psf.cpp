// @COPYRIGHT@

#include <harp_internal.hpp>


using namespace std;
using namespace harp;

static const char * format_sandbox = "sandbox";

static const char * sandbox_psf_key_path = "path";
static const char * sandbox_psf_key_corr = "corr";

static const char * sandbox_psf_key_name = "PSFPARAM";

static const char * sandbox_psf_hdu_x = "X";
static const char * sandbox_psf_hdu_y = "Y";
static const char * sandbox_psf_hdu_lambda = "Wavelength";
static const char * sandbox_psf_hdu_amp = "Amplitude";
static const char * sandbox_psf_hdu_maj = "MajorAxis";
static const char * sandbox_psf_hdu_min = "MinorAxis";
static const char * sandbox_psf_hdu_ang = "Angle";




harp::psf_sandbox::psf_sandbox ( boost::property_tree::ptree const & props ) : psf ( format_sandbox, props ) {
  
  path_ = props.get < string > ( sandbox_psf_key_path );

  pixcorr_ = props.get < int > ( sandbox_psf_key_corr );

  
  fitsfile *fp;

  fits::open_read ( fp, path_ );

  int hdu = fits::img_seek ( fp, sandbox_psf_key_name, sandbox_psf_hdu_x );
  hdus_[ sandbox_psf_hdu_x ] = hdu;
  fits::img_dims ( fp, nspec_, specsize_ );

  size_t rows, cols;
  
  hdu = fits::img_seek ( fp, sandbox_psf_key_name, sandbox_psf_hdu_y );
  hdus_[ sandbox_psf_hdu_y ] = hdu;
  fits::img_dims ( fp, rows, cols );
  if ( ( rows != nspec_ ) || ( cols != specsize_ ) ) {
    MOAT_THROW( "sandbox_psf: PSF file must have identical dimensions for all HDUs" );
  }

  hdu = fits::img_seek ( fp, sandbox_psf_key_name, sandbox_psf_hdu_lambda );
  hdus_[ sandbox_psf_hdu_lambda ] = hdu;
  fits::img_dims ( fp, rows, cols );
  if ( ( rows != nspec_ ) || ( cols != specsize_ ) ) {
    MOAT_THROW( "sandbox_psf: PSF file must have identical dimensions for all HDUs" );
  }

  hdu = fits::img_seek ( fp, sandbox_psf_key_name, sandbox_psf_hdu_amp );
  hdus_[ sandbox_psf_hdu_amp ] = hdu;
  fits::img_dims ( fp, rows, cols );
  if ( ( rows != nspec_ ) || ( cols != specsize_ ) ) {
    MOAT_THROW( "sandbox_psf: PSF file must have identical dimensions for all HDUs" );
  }
  
  hdu = fits::img_seek ( fp, sandbox_psf_key_name, sandbox_psf_hdu_maj );
  hdus_[ sandbox_psf_hdu_maj ] = hdu;
  fits::img_dims ( fp, rows, cols );
  if ( ( rows != nspec_ ) || ( cols != specsize_ ) ) {
    MOAT_THROW( "sandbox_psf: PSF file must have identical dimensions for all HDUs" );
  }
  
  hdu = fits::img_seek ( fp, sandbox_psf_key_name, sandbox_psf_hdu_min );
  hdus_[ sandbox_psf_hdu_min ] = hdu;
  fits::img_dims ( fp, rows, cols );
  if ( ( rows != nspec_ ) || ( cols != specsize_ ) ) {
    MOAT_THROW( "sandbox_psf: PSF file must have identical dimensions for all HDUs" );
  }
  
  hdu = fits::img_seek ( fp, sandbox_psf_key_name, sandbox_psf_hdu_ang );
  hdus_[ sandbox_psf_hdu_ang ] = hdu;
  fits::img_dims ( fp, rows, cols );
  if ( ( rows != nspec_ ) || ( cols != specsize_ ) ) {
    MOAT_THROW( "sandbox_psf: PSF file must have identical dimensions for all HDUs" );
  }
  
  fits::close ( fp );

  nglobal_ = nspec_ * specsize_;
  
}


harp::psf_sandbox::~psf_sandbox ( ) {
  
  cleanup();
  
}


boost::property_tree::ptree harp::psf_sandbox::serialize ( ) {
  boost::property_tree::ptree ret;

  ret.put ( "format", psf::format() );

  ret.put ( sandbox_psf_key_path, path_ );

  ret.put ( sandbox_psf_key_corr, pixcorr_ );

  return ret;
}


void harp::psf_sandbox::cache_spec ( size_t first, size_t last ) {
  
  fitsfile *fp = NULL;

  for ( size_t spec = first; spec <= last; ++spec ) {
    map < size_t, psf_sandbox_resp > :: iterator check;
    
    check = resp_.find ( spec );

    if ( check == resp_.end() ) {
      //cerr << "cache_spec reading spectrum " << spec << endl;
      if ( ! fp ) {
        //cerr << "cache_spec opening file " << path_ << endl;
        fits::open_read ( fp, path_ );
      }
      
      psf_sandbox_resp temp;
      resp_[ spec ] = temp;
      
      //cerr << "cache_spec resizing vectors to " << specsize_ << " elements" << endl;
      
      resp_[ spec ].x.resize ( specsize_ );
      resp_[ spec ].y.resize ( specsize_ );
      resp_[ spec ].lambda.resize ( specsize_ );
      resp_[ spec ].amp.resize ( specsize_ );
      resp_[ spec ].maj.resize ( specsize_ );
      resp_[ spec ].min.resize ( specsize_ );
      resp_[ spec ].ang.resize ( specsize_ );
      
      fits::img_seek ( fp, hdus_[ sandbox_psf_hdu_x ] );      
      fits::img_read_row ( fp, spec, resp_[ spec ].x );
      
      fits::img_seek ( fp, hdus_[ sandbox_psf_hdu_y ] );      
      fits::img_read_row ( fp, spec, resp_[ spec ].y );
      
      fits::img_seek ( fp, hdus_[ sandbox_psf_hdu_lambda ] );      
      fits::img_read_row ( fp, spec, resp_[ spec ].lambda );
      
      fits::img_seek ( fp, hdus_[ sandbox_psf_hdu_amp ] );      
      fits::img_read_row ( fp, spec, resp_[ spec ].amp );
      
      fits::img_seek ( fp, hdus_[ sandbox_psf_hdu_maj ] );      
      fits::img_read_row ( fp, spec, resp_[ spec ].maj );
      
      fits::img_seek ( fp, hdus_[ sandbox_psf_hdu_min ] );      
      fits::img_read_row ( fp, spec, resp_[ spec ].min );
      
      fits::img_seek ( fp, hdus_[ sandbox_psf_hdu_ang ] );      
      fits::img_read_row ( fp, spec, resp_[ spec ].ang );
      
      //for ( size_t i = 0; i < specsize_; ++i ) {
        //cout << "spec[" << spec << "](" << i << ") = " << resp_[ spec ].x[i] << " " << resp_[ spec ].y[i] << " " << resp_[ spec ].lambda[i] << " " << resp_[ spec ].amp[i] << " " << resp_[ spec ].maj[i] << " " << resp_[ spec ].min[i] << " " << resp_[ spec ].ang[i] << endl;
      //}
      
    }
  }
  
  if ( fp ) {
    fits::close ( fp );
  }
  
  return;
}


void harp::psf_sandbox::lambda ( size_t specnum, vec_dense & data ) {

  cache_spec ( specnum, specnum );
  
  data.resize ( specsize_ );
  
  for ( size_t i = 0; i < specsize_; ++i ) {
    data[ i ] = resp_[ specnum ].lambda[ i ];
  }
  
  return;
}


void harp::psf_sandbox::extent ( size_t firstspec, size_t lastspec, size_t firstbin, size_t lastbin, size_t & firstX, size_t & firstY, size_t & lastX, size_t & lastY ) {
  
  cache_spec ( firstspec, lastspec );
  
  double minX = 1000000000.0;
  double minY = 1000000000.0;
  double maxX = 0.0;
  double maxY = 0.0;
  
  for ( size_t curspec = firstspec; curspec <= lastspec; ++curspec ) {
    for ( size_t curbin = firstbin; curbin <= lastbin; ++curbin ) {
      if ( resp_[ curspec ].x[ curbin ] < minX ) {
        minX = resp_[ curspec ].x[ curbin ];
      }
      if ( resp_[ curspec ].x[ curbin ] > maxX ) {
        maxX = resp_[ curspec ].x[ curbin ];
      }
      if ( resp_[ curspec ].y[ curbin ] < minY ) {
        minY = resp_[ curspec ].y[ curbin ];
      }
      if ( resp_[ curspec ].y[ curbin ] > maxY ) {
        maxY = resp_[ curspec ].y[ curbin ];
      }
    }
  }
  
  if ( (int)minX - (int)pixcorr_ < 0 ) {
    firstX = 0;
  } else {
    firstX = (int)minX - (int)pixcorr_;
  }
  
  if ( (int)minY - (int)pixcorr_ < 0 ) {
    firstY = 0;
  } else {
    firstY = (int)minY - (int)pixcorr_;
  }
  
  lastX = (int)maxX + (int)pixcorr_;
  lastY = (int)maxY + (int)pixcorr_;
  
  return;
}


void harp::psf_sandbox::gauss_sample ( vec_dense & vals, vec_dense & xrel, vec_dense & yrel, double amp, double maj, double min, double ang ) {
  
  amp /= maj * min * moat::TWOPI;
  
  size_t nvals = xrel.size();
  
  double cang = cos ( ang );
  double sang = sin ( ang );
  
  double invmaj = 1.0 / maj;
  double invmin = 1.0 / min;
  
  double * buf = moat::double_alloc ( nvals );
  
  double xt, yt;
  
  size_t i;
  
  #ifdef _OPENMP
  #pragma omp parallel for default(none) private(i, xt, yt) shared(nvals, xrel, yrel, buf, cang, sang, invmaj, invmin) schedule(static)
  #endif
  for ( i = 0; i < nvals; ++i ) {
    xt = xrel[i] * cang + yrel[i] * sang;
    yt = - xrel[i] * sang + yrel[i] * cang;
    buf[i] = - 0.5 * ( xt * xt * invmaj * invmaj + yt * yt * invmin * invmin );
  }
  
  moat::sf::fast_exp ( nvals, buf, &(vals[0]) );
  
  #ifdef _OPENMP
  #pragma omp parallel for default(none) private(i) shared(nvals, vals, amp) schedule(static)
  #endif
  for ( i = 0; i < nvals; ++i ) {
    vals[i] *= amp;
  }
  
  free ( buf );
  
  return;
}


void harp::psf_sandbox::gauss_sample_alt ( vec_dense & vals, vec_dense & xrel, vec_dense & yrel, double amp, double maj, double min, double ang ) {
  
  amp /= maj * min * moat::TWOPI;
  
  size_t nvals = xrel.size();
  
  double * buf = moat::double_alloc ( nvals );
  
  double cang = cos ( ang );
  double sang = sin ( ang );
  double s2ang = sin ( 2.0 * ang );
  
  double mjmj = 0.5 / ( maj * maj );
  double mnmn = 0.5 / ( min * min );
  
  double cxx = mjmj * cang * cang + mnmn * sang * sang;
  double cxy = 0.5 * mjmj * s2ang - 0.5 * mnmn * s2ang;
  double cyy  = mjmj * sang * sang + mnmn * cang * cang;

  size_t i;
  
  #ifdef _OPENMP
  #pragma omp parallel for default(none) private(i) shared(nvals, xrel, yrel, buf, cxx, cxy, cyy) schedule(static)
  #endif
  for ( i = 0; i < nvals; ++i ) {
    buf[i] = - ( cxx * xrel[i] * xrel[i] + 2.0 * cxy * xrel[i] * yrel[i] + cyy * yrel[i] * yrel[i] );
  }
  
  moat::sf::fast_exp ( nvals, buf, &(vals[0]) );
  
  #ifdef _OPENMP
  #pragma omp parallel for default(none) private(i) shared(nvals, vals, amp) schedule(static)
  #endif
  for ( i = 0; i < nvals; ++i ) {
    vals[i] *= amp;
  }
  
  free ( buf );

  return;
}


size_t harp::psf_sandbox::valid_range ( size_t const & firstX, size_t const & lastX, size_t const & firstY, size_t const & lastY, size_t & startX, size_t & stopX, size_t & startY, size_t & stopY, size_t & spec, size_t & bin ) {

  bool valid = true;
  
  extent ( spec, spec, bin, bin, startX, startY, stopX, stopY );
  
  if ( startX < firstX ) {
    startX = firstX;
  } else if ( startX > lastX ) {
    valid = false;
  }
  
  if ( stopX > lastX ) {
    stopX = lastX;
  } else if ( stopX < firstX ) {
    valid = false;
  }
  
  if ( startY < firstY ) {
    startY = firstY;
  } else if ( startY > lastY ) {
    valid = false;
  }
  
  if ( stopY > lastY ) {
    stopY = lastY;
  } else if ( stopY < firstY ) {
    valid = false;
  }
  
  size_t nz = 0;
  if ( valid ) {
    nz = ( stopY - startY + 1 ) * ( stopX - startX + 1);
  }
    
  return nz;
}


void harp::psf_sandbox::projection ( string profcalc, string profremap, size_t firstspec, size_t lastspec, size_t firstbin, size_t lastbin, size_t firstX, size_t lastX, size_t firstY, size_t lastY, mat_compcol & data ) {
  
  size_t nx = lastX - firstX + 1;
  size_t ny = lastY - firstY + 1;
  
  size_t specsize = data.size2();
  
  size_t npix = data.size1();
  
  if ( ( npix != nx * ny ) || ( specsize != ( lastspec - firstspec + 1 ) * ( lastbin - firstbin + 1 ) ) ) {
    std::ostringstream o;
    o << "sandbox_psf: PSF projection ranges must match dimensions of projection data (" << npix << " x " << specsize << ")";
    MOAT_THROW( o.str().c_str() );
  }

  moat::profile * prof = moat::profile::get ( );
  if ( profcalc != "" ) {
    prof->start ( profcalc );
  }
  
  cache_spec ( firstspec, lastspec );
  
  vector < size_t > binlist ( specsize );
  
  size_t b = 0;
  for ( size_t i = firstspec; i <= lastspec; ++i ) {
    for ( size_t j = firstbin; j <= lastbin; ++j ) {
      binlist[b] = i * specsize_ + j;
      ++b;
    }
  }
  
  //cerr << "computing number of non-zeros" << endl;
  
  size_t nonzeros = 0;
  
  #ifdef _OPENMP
  #pragma omp parallel for default(none) private(b) shared(specsize, binlist, firstX, firstY, lastX, lastY, nonzeros) schedule(static)
  #endif
  for ( b = 0; b < specsize; ++b ) {
    
    size_t specbin = binlist[b] % specsize_;
    size_t spec = (size_t)(binlist[b] / specsize_);
    
    size_t startX;
    size_t stopX;
    size_t startY;
    size_t stopY;
    
    size_t nz = valid_range ( firstX, lastX, firstY, lastY, startX, stopX, startY, stopY, spec, specbin );
    
    #ifdef _OPENMP
    #pragma omp critical
    #endif
    {
      nonzeros += nz;
    }
  }
  
  //cerr << "  = " << nonzeros << endl;
  
  
  // We want to fill the matrix using a mapped matrix, but the output will be a 
  // compressed matrix to speed up axpy_prod computations.  We use a temporary 
  // matrix and then assign to the output.
  
  mat_dyncol builder ( npix, specsize );

  int lastfrac = 0;
  size_t complete = 0;

  //fprintf ( stderr, "  Sampling PSF [          ]\r" );
  
  for ( b = 0; b < specsize; ++b ) {
    
    size_t specbin = binlist[b] % specsize_;
    
    size_t spec = (size_t)(binlist[b] / specsize_);
    
    size_t startX;
    size_t stopX;
    size_t startY;
    size_t stopY;
    
    size_t nvalid = valid_range ( firstX, lastX, firstY, lastY, startX, stopX, startY, stopY, spec, specbin );

    if ( nvalid > 0 ) {
      
      double amp = resp_[ spec ].amp[ specbin ];
      double maj = resp_[ spec ].maj[ specbin ];
      double min = resp_[ spec ].min[ specbin ];
      double ang = resp_[ spec ].ang[ specbin ];
      double xcenter = resp_[ spec ].x[ specbin ];
      double ycenter = resp_[ spec ].y[ specbin ];
    
      vec_dense fxdist ( nvalid );
      vec_dense fydist ( nvalid );
    
      size_t pix = 0;
    
      //cerr << "  computing distances" << endl;
    
      for ( size_t imgrow = startY; imgrow <= stopY; ++imgrow ) {
      
        for ( size_t imgcol = startX; imgcol <= stopX; ++imgcol ) {
      
          fxdist[pix] = (double)imgcol - xcenter;
          fydist[pix] = (double)imgrow - ycenter;
        
          ++pix;
        }
      
      }
    
      vec_dense vals ( nvalid );

      gauss_sample ( vals, fxdist, fydist, amp, maj, min, ang );
    
      size_t datarow;
    
      pix = 0;
    
      size_t rowoff;
    
      for ( size_t imgrow = startY; imgrow <= stopY; ++imgrow ) {
        rowoff = ( imgrow - firstY ) * ( lastX - firstX + 1 );
    
        for ( size_t imgcol = startX; imgcol <= stopX; ++imgcol ) {
      
          datarow = rowoff + ( imgcol - firstX );

          builder ( datarow, b ) += vals[pix];
        
          ++pix;
        }
      }
      
    }

    ++complete;

    char msg[256];
    int progfrac = (int) ( 10 * complete / specsize );
    /*
    if ( progfrac != lastfrac ) {
      for ( int p = 0; p < progfrac; ++p ) {
        msg[p] = '*';
      }
      for ( int p = progfrac; p < 10; ++p ) {
        msg[p] = ' ';
      }
      msg[10] = '\0';    
      fprintf ( stderr, "  Sampling PSF [%s]\r", msg );
    }
    */
    lastfrac = progfrac;
    
  }
  
  if ( profcalc != "" ) {
    prof->stop ( profcalc );
  }

  if ( profremap != "" ) {
    prof->start ( profremap );
  }
  
  // copy to output matrix
  
  complete = 0;
  char msg[256];
  int progfrac;
  lastfrac = 0;

  data.reserve ( nonzeros );

  mat_dyncol::iterator2 itcol;
  mat_dyncol::iterator1 itrow;

  for ( itcol = builder.begin2(); itcol != builder.end2(); ++itcol ) {
    for ( itrow = itcol.begin(); itrow != itcol.end(); ++itrow ) {
      data ( itrow.index1(), itrow.index2() ) = (*itrow);
      ++complete;
    }
    progfrac = (int) ( 10 * complete / nonzeros );
    
    /*
    if ( progfrac != lastfrac ) {
      for ( int p = 0; p < progfrac; ++p ) {
        msg[p] = '*';
      }
      for ( int p = progfrac; p < 10; ++p ) {
        msg[p] = ' ';
      }
      msg[10] = '\0';    
      fprintf ( stderr, "  Memory Remap [%s]\r", msg );
    }
    */
    lastfrac = progfrac;
  }
  
  
  if ( profremap != "" ) {
    prof->stop ( profremap );
  }
  
  return;
}




