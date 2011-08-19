// @COPYRIGHT@

#include <harp_internal.hpp>


using namespace std;
using namespace harp;

static const char * format_boss = "boss";

static const char * boss_psf_key_path = "path";
static const char * boss_psf_key_corr = "corr";


static const char * boss_psf_hdu_name = "PSFPARAM";
static const char * boss_psf_type_name = "PSFTYPE";
static const char * boss_psf_xpix_name = "NPIX_X";
static const char * boss_psf_ypix_name = "NPIX_Y";


static const char * boss_gausspsf_hdu_x = "X";
static const char * boss_gausspsf_hdu_y = "Y";
static const char * boss_gausspsf_hdu_lambda = "Wavelength";
static const char * boss_gausspsf_hdu_amp = "Amplitude";
static const char * boss_gausspsf_hdu_maj = "MajorAxis";
static const char * boss_gausspsf_hdu_min = "MinorAxis";
static const char * boss_gausspsf_hdu_ang = "Angle";


static const int boss_pixpsf_hdu_x = 1;
static const int boss_pixpsf_hdu_y = 2;
static const int boss_pixpsf_hdu_lambda = 3;
static const int boss_pixpsf_hdu_nexp = 4;      // icoeff xexp yexp
static const int boss_pixpsf_hdu_xyscale = 5;   // ifiber igroup x0 xscale y0 yscale
static const int boss_pixpsf_hdu_psfimage = 6;  // igroup icoeff iy ix


harp::psf_boss::psf_boss ( std::map < std::string, std::string > const & params ) : psf ( format_boss, params ) {
  
  map < std::string, std::string > :: const_iterator val;
  
  val = params.find( boss_psf_key_path );
  
  if ( val == params.end() ) {
    MOAT_THROW( "boss_psf: must specify path to PSF file" );
  }
    
  path_ = val->second;
  
  fitsfile *fp;
  int hdu;
  size_t rows, cols;
  fits::open_read ( fp, path_ );
  
  // determine PSF type
  
  fits::img_seek ( fp, 1 );
  type_ = fits::key_string ( fp, boss_psf_type_name );
  
  if ( type_ == "GAUSS2D" ) {
    
    val = params.find( boss_psf_key_corr );

    if ( val == params.end() ) {
      MOAT_THROW( "boss_psf: must specify pixel space correlation length for GAUSS2D psf" );
    }

    xpixcorr_ = atoi ( val->second.c_str() );
    ypixcorr_ = xpixcorr_;
    
    xpixwidth_ = 2 * xpixcorr_ + 1;
    ypixwidth_ = 2 * ypixcorr_ + 1;
  
    hdu = fits::img_seek ( fp, boss_psf_hdu_name, boss_gausspsf_hdu_x );
    hdus_[ boss_gausspsf_hdu_x ] = hdu;
    
    // read spectral size
    
    fits::img_dims ( fp, nspec_, nbins_ );
    xpix_ = fits::key_long ( fp, boss_psf_xpix_name );
    ypix_ = fits::key_long ( fp, boss_psf_ypix_name );
  
    hdu = fits::img_seek ( fp, boss_psf_hdu_name, boss_gausspsf_hdu_y );
    hdus_[ boss_gausspsf_hdu_y ] = hdu;
    fits::img_dims ( fp, rows, cols );
    if ( ( rows != nspec_ ) || ( cols != nbins_ ) ) {
      MOAT_THROW( "boss_psf: PSF file must have identical dimensions for all HDUs" );
    }
  
    hdu = fits::img_seek ( fp, boss_psf_hdu_name, boss_gausspsf_hdu_lambda );
    hdus_[ boss_gausspsf_hdu_lambda ] = hdu;
    fits::img_dims ( fp, rows, cols );
    if ( ( rows != nspec_ ) || ( cols != nbins_ ) ) {
      MOAT_THROW( "boss_psf: PSF file must have identical dimensions for all HDUs" );
    }
  
    hdu = fits::img_seek ( fp, boss_psf_hdu_name, boss_gausspsf_hdu_amp );
    hdus_[ boss_gausspsf_hdu_amp ] = hdu;
    fits::img_dims ( fp, rows, cols );
    if ( ( rows != nspec_ ) || ( cols != nbins_ ) ) {
      MOAT_THROW( "boss_psf: PSF file must have identical dimensions for all HDUs" );
    }
  
    hdu = fits::img_seek ( fp, boss_psf_hdu_name, boss_gausspsf_hdu_maj );
    hdus_[ boss_gausspsf_hdu_maj ] = hdu;
    fits::img_dims ( fp, rows, cols );
    if ( ( rows != nspec_ ) || ( cols != nbins_ ) ) {
      MOAT_THROW( "boss_psf: PSF file must have identical dimensions for all HDUs" );
    }
  
    hdu = fits::img_seek ( fp, boss_psf_hdu_name, boss_gausspsf_hdu_min );
    hdus_[ boss_gausspsf_hdu_min ] = hdu;
    fits::img_dims ( fp, rows, cols );
    if ( ( rows != nspec_ ) || ( cols != nbins_ ) ) {
      MOAT_THROW( "boss_psf: PSF file must have identical dimensions for all HDUs" );
    }
  
    hdu = fits::img_seek ( fp, boss_psf_hdu_name, boss_gausspsf_hdu_ang );
    hdus_[ boss_gausspsf_hdu_ang ] = hdu;
    fits::img_dims ( fp, rows, cols );
    if ( ( rows != nspec_ ) || ( cols != nbins_ ) ) {
      MOAT_THROW( "boss_psf: PSF file must have identical dimensions for all HDUs" );
    }
    
  } else if ( type_ == "PCA-PIX" ) {
    
    // read spectral size
    
    fits::img_seek ( fp, boss_pixpsf_hdu_x );
    
    fits::img_dims ( fp, nspec_, nbins_ );
    xpix_ = fits::key_long ( fp, boss_psf_xpix_name );
    ypix_ = fits::key_long ( fp, boss_psf_ypix_name );
    
    // read correlation size, PCA coefficients, and number of groups
    
    fits::img_seek ( fp, boss_pixpsf_hdu_psfimage );
    xpixwidth_ = fits::key_long ( fp, "NAXIS1" );
    ypixwidth_ = fits::key_long ( fp, "NAXIS2" );
    ncoeff_ = fits::key_long ( fp, "NAXIS3" );
    ngroup_ = fits::key_long ( fp, "NAXIS4" );
    
    xpixcorr_ = (size_t)(xpixwidth_ - 1) / 2;
    ypixcorr_ = (size_t)(ypixwidth_ - 1) / 2;
    
    int ret;
    int status = 0;
    long fpixel[4];
    long lpixel[4];
    long inc[4] = {1, 1, 1, 1};
    int anynul;

    fpixel[0] = 1;
    fpixel[1] = 1;
    fpixel[2] = 1;
    fpixel[3] = 1;

    lpixel[0] = (long)xpixwidth_;
    lpixel[1] = (long)ypixwidth_;
    lpixel[2] = (long)ncoeff_;
    lpixel[3] = (long)ngroup_;
    
    size_t bufsize = xpixwidth_ * ypixwidth_ * ncoeff_ * ngroup_;
    
    double * buffer = moat::double_alloc ( bufsize );

    ret = fits_read_subset ( fp, TDOUBLE, fpixel, lpixel, inc, 0, buffer, &anynul, &status );
    fits::check ( status );
    
    size_t global = 0;
    for ( size_t i = 0; i < ngroup_; ++i ) {
      map < size_t, data_vec > curimage;
      for ( size_t j = 0; j < ncoeff_; ++j ) {
        curimage[ j ] = data_vec ( xpixwidth_ * ypixwidth_ );
        for ( size_t k = 0; k < xpixwidth_ * ypixwidth_; ++k ) {
          (curimage[ j ])[ k ] = buffer[ global ];
          ++global;
        }
         
      }
      psfimage_[ i ] = curimage;
    }
    
    free ( buffer );    
    
    // read x/y exponents
    
    fits::bin_seek ( fp, boss_pixpsf_hdu_nexp );
    
    vector < string > cols;
    cols.push_back ( "XEXP" );
    cols.push_back ( "YEXP" );
    vector < int > colids = fits::bin_columns ( fp, cols );
    
    vector < data_vec > tempdata ( 2 );
    fits::bin_read ( fp, 0, ncoeff_-1, colids, tempdata );
    
    xexp_ = tempdata[0];
    yexp_ = tempdata[1];
    
    // read scale factors
    
    fits::bin_seek ( fp, boss_pixpsf_hdu_xyscale );
    
    cols.clear();
    cols.push_back ( "IGROUP" );
    cols.push_back ( "X0" );
    cols.push_back ( "XSCALE" );
    cols.push_back ( "Y0" );
    cols.push_back ( "YSCALE" );
    colids = fits::bin_columns ( fp, cols );
    
    tempdata.clear();
    tempdata.resize ( 5 );
    fits::bin_read ( fp, 0, nspec_-1, colids, tempdata );
    
    psf_boss_pixscale scale;
    
    for ( size_t i = 0; i < nspec_; ++i ) {
      scale.igroup = (size_t) (tempdata[0])[i];
      scale.x0 = (tempdata[1])[i];
      scale.xscale = (tempdata[2])[i];
      scale.y0 = (tempdata[3])[i];
      scale.yscale = (tempdata[4])[i];
      xyscale_[ i ] = scale;
      
      //char msg[256];
      //sprintf(msg,"group = %lu  x0 = %0.16e  y0 = %0.16e  xscale = %0.16e  yscale = %0.16e", scale.igroup, scale.x0, scale.y0, scale.xscale, scale.yscale );
      //cerr << msg << endl;
    }
    
  } else if ( type_ == "GAUSS-HERMITE" ) {

    MOAT_THROW( "GAUSS-HERMITE psf not yet supported" );

  } else {
    MOAT_THROW( "unknown type in boss psf projection" );
  }
  
  fits::close ( fp );
  
}


harp::psf_boss::~psf_boss ( ) {
  
  cleanup();
  
}


void harp::psf_boss::cache_spec ( size_t first, size_t last ) {
  
  fitsfile *fp = NULL;

  for ( size_t spec = first; spec <= last; ++spec ) {
    map < size_t, psf_boss_resp > :: iterator check;
    
    check = resp_.find ( spec );

    if ( check == resp_.end() ) {
      //cerr << "cache_spec reading spectrum " << spec << endl;
      if ( ! fp ) {
        //cerr << "cache_spec opening file " << path_ << endl;
        fits::open_read ( fp, path_ );
      }
      
      psf_boss_resp temp;
      resp_[ spec ] = temp;
      
      if ( type_ == "GAUSS2D" ) {
      
        //cerr << "cache_spec resizing vectors to " << nbins_ << " elements" << endl;
      
        resp_[ spec ].x.resize ( nbins_ );
        resp_[ spec ].y.resize ( nbins_ );
        resp_[ spec ].lambda.resize ( nbins_ );
        resp_[ spec ].amp.resize ( nbins_ );
        resp_[ spec ].maj.resize ( nbins_ );
        resp_[ spec ].min.resize ( nbins_ );
        resp_[ spec ].ang.resize ( nbins_ );
      
        fits::img_seek ( fp, hdus_[ boss_gausspsf_hdu_x ] );      
        fits::img_read_row ( fp, spec, resp_[ spec ].x );
      
        fits::img_seek ( fp, hdus_[ boss_gausspsf_hdu_y ] );      
        fits::img_read_row ( fp, spec, resp_[ spec ].y );
      
        fits::img_seek ( fp, hdus_[ boss_gausspsf_hdu_lambda ] );      
        fits::img_read_row ( fp, spec, resp_[ spec ].lambda );
      
        fits::img_seek ( fp, hdus_[ boss_gausspsf_hdu_amp ] );      
        fits::img_read_row ( fp, spec, resp_[ spec ].amp );
      
        fits::img_seek ( fp, hdus_[ boss_gausspsf_hdu_maj ] );      
        fits::img_read_row ( fp, spec, resp_[ spec ].maj );
      
        fits::img_seek ( fp, hdus_[ boss_gausspsf_hdu_min ] );      
        fits::img_read_row ( fp, spec, resp_[ spec ].min );
      
        fits::img_seek ( fp, hdus_[ boss_gausspsf_hdu_ang ] );      
        fits::img_read_row ( fp, spec, resp_[ spec ].ang );
      
      
        //for ( size_t i = 0; i < nbins_; ++i ) {
          //cout << "spec[" << spec << "](" << i << ") = " << resp_[ spec ].x[i] << " " << resp_[ spec ].y[i] << " " << resp_[ spec ].lambda[i] << " " << resp_[ spec ].amp[i] << " " << resp_[ spec ].maj[i] << " " << resp_[ spec ].min[i] << " " << resp_[ spec ].ang[i] << endl;
        //}
      
      } else if ( type_ == "PCA-PIX" ) { 

        resp_[ spec ].x.resize ( nbins_ );
        resp_[ spec ].y.resize ( nbins_ );
        resp_[ spec ].lambda.resize ( nbins_ );
      
        fits::img_seek ( fp, boss_pixpsf_hdu_x );      
        fits::img_read_row ( fp, spec, resp_[ spec ].x );
      
        fits::img_seek ( fp, boss_pixpsf_hdu_y );      
        fits::img_read_row ( fp, spec, resp_[ spec ].y );
      
        fits::img_seek ( fp, boss_pixpsf_hdu_lambda );      
        fits::img_read_row ( fp, spec, resp_[ spec ].lambda );
        
      } else if ( type_ == "GAUSS-HERMITE" ) {



      } else {
        MOAT_THROW( "unknown type in boss psf projection" );
      }
      
      
    }
  }
  
  if ( fp ) {
    fits::close ( fp );
  }
  
  return;
}


void harp::psf_boss::lambda ( size_t specnum, data_vec & data ) {
  
  return;
}


void harp::psf_boss::extent ( size_t firstspec, size_t lastspec, size_t firstbin, size_t lastbin, size_t & firstX, size_t & firstY, size_t & lastX, size_t & lastY ) {
  
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
  
  if ( (int)minX - (int)xpixcorr_ < 0 ) {
    firstX = 0;
  } else {
    firstX = (int)minX - (int)xpixcorr_;
  }
  
  if ( (int)minY - (int)ypixcorr_ < 0 ) {
    firstY = 0;
  } else {
    firstY = (int)minY - (int)ypixcorr_;
  }
  
  lastX = (int)maxX + (int)xpixcorr_;
  lastY = (int)maxY + (int)ypixcorr_;
  
  return;
}


size_t harp::psf_boss::valid_range ( size_t const & firstX, size_t const & lastX, size_t const & firstY, size_t const & lastY, size_t & startX, size_t & stopX, size_t & startY, size_t & stopY, size_t & spec, size_t & bin ) {

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


void harp::psf_boss::gauss_sample ( data_vec & vals, data_vec & xrel, data_vec & yrel, double amp, double maj, double min, double ang ) {
  
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


void harp::psf_boss::projection ( string profcalc, string profremap, size_t firstspec, size_t lastspec, size_t firstbin, size_t lastbin, size_t firstX, size_t lastX, size_t firstY, size_t lastY, comp_rowmat & data ) {
  
  size_t nx = lastX - firstX + 1;
  size_t ny = lastY - firstY + 1;
  
  size_t nbins = data.size2();
  
  size_t npix = data.size1();
  
  if ( ( npix != nx * ny ) || ( nbins != ( lastspec - firstspec + 1 ) * ( lastbin - firstbin + 1 ) ) ) {
    std::ostringstream o;
    o << "boss_psf: PSF projection ranges must match dimensions of projection data (" << npix << " x " << nbins << ")";
    MOAT_THROW( o.str().c_str() );
  }

  moat::profile * prof = moat::profile::get ( );
  if ( profcalc != "" ) {
    prof->start ( profcalc );
  }
  
  cache_spec ( firstspec, lastspec );
  
  vector < size_t > binlist ( nbins );
  
  size_t b = 0;
  for ( size_t i = firstspec; i <= lastspec; ++i ) {
    for ( size_t j = firstbin; j <= lastbin; ++j ) {
      binlist[b] = i * nbins_ + j;
      ++b;
    }
  }
  
  //cerr << "computing number of non-zeros" << endl;
  
  size_t nonzeros = 0;
  
  #ifdef _OPENMP
  #pragma omp parallel for default(none) private(b) shared(nbins, binlist, firstX, firstY, lastX, lastY, nonzeros) schedule(static)
  #endif
  for ( b = 0; b < nbins; ++b ) {
    
    size_t specbin = binlist[b] % nbins_;
    size_t spec = (size_t)(binlist[b] / nbins_);
    
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
  
  
  // We want to fill the matrix using a mapped matrix, but the output will be a compressed matrix to speed up
  // axpy_prod computations.  We use a temporary matrix and then assign to the output.
  
  sparse_rowmat builder ( npix, nbins );
  
  int lastfrac = 0;
  size_t complete = 0;

  //fprintf ( stderr, "  Sampling PSF [          ]\r" );
  
  for ( b = 0; b < nbins; ++b ) {
    
    size_t specbin = binlist[b] % nbins_;
    
    size_t spec = (size_t)(binlist[b] / nbins_);
    
    size_t startX;
    size_t stopX;
    size_t startY;
    size_t stopY;
    
    size_t nvalid = valid_range ( firstX, lastX, firstY, lastY, startX, stopX, startY, stopY, spec, specbin );

    if ( nvalid > 0 ) {
      
      double xcenter;
      double ycenter;
      size_t pix;
      size_t datarow;
      size_t rowoff;
      data_vec vals ( nvalid );
      
      //cerr << "spec " << spec << ", bin " << specbin << " has " << nvalid << " nonzeroes" << endl;
      
      if ( type_ == "GAUSS2D" ) {
      
        double amp = resp_[ spec ].amp[specbin];
        double maj = resp_[ spec ].maj[specbin];
        double min = resp_[ spec ].min[specbin];
        double ang = resp_[ spec ].ang[specbin];
        xcenter = resp_[ spec ].x[specbin];
        ycenter = resp_[ spec ].y[specbin];
    
        data_vec fxdist ( nvalid );
        data_vec fydist ( nvalid );
    
        pix = 0;
    
        //cerr << "  computing distances" << endl;
    
        for ( size_t imgrow = startY; imgrow <= stopY; ++imgrow ) {
      
          for ( size_t imgcol = startX; imgcol <= stopX; ++imgcol ) {
      
            fxdist[pix] = (double)imgcol - xcenter;
            fydist[pix] = (double)imgrow - ycenter;
        
            ++pix;
          }
      
        }

        gauss_sample ( vals, fxdist, fydist, amp, maj, min, ang );
    
        pix = 0;
    
        for ( size_t imgrow = startY; imgrow <= stopY; ++imgrow ) {
          rowoff = ( imgrow - firstY ) * ( lastX - firstX + 1 );
    
          for ( size_t imgcol = startX; imgcol <= stopX; ++imgcol ) {
      
            datarow = rowoff + ( imgcol - firstX );

            builder ( datarow, b ) += vals[pix];
        
            ++pix;
          }
        }
        
      } else if ( type_ == "PCA-PIX" ) {
        
        xcenter = resp_[ spec ].x[specbin];
        ycenter = resp_[ spec ].y[specbin];
        
        size_t igroup = xyscale_[ spec ].igroup;
        double x0 = xyscale_[ spec ].x0;
        double y0 = xyscale_[ spec ].y0;
        double xscale = xyscale_[ spec ].xscale;
        double yscale = xyscale_[ spec ].yscale;
        
        double xx = xscale * ( xcenter - x0 );
        double yy = yscale * ( ycenter - y0 );
        
        //if ( ( spec == 50 ) && ( specbin == 50 ) ) {
        
          cerr << "PIX_PCA:  spec " << spec << ", bin " << specbin << ":" << endl;
          cerr << "PIX_PCA:    igroup = " << igroup << endl; 
          cerr << "PIX_PCA:    x0 = " << x0 << endl;
          cerr << "PIX_PCA:    xscale = " << xscale << endl;
          cerr << "PIX_PCA:    y0 = " << y0 << endl;
          cerr << "PIX_PCA:    yscale = " << yscale << endl;
          cerr << "PIX_PCA:    xx = " << xx << endl;
          cerr << "PIX_PCA:    yy = " << yy << endl;
          
        //}
        
        int nx;
        int ny;
        
        for ( pix = 0; pix < nvalid; ++pix ) {
          vals [ pix ] = 0.0;
        }
        
        for ( size_t co = 0; co < ncoeff_; ++co ) {
          nx = xexp_[ co ];
          ny = yexp_[ co ];
          
          //if ( specbin == 0 ) {
          //  cerr << "PIX_PCA:    nx = " << nx << " ny = " << ny << endl;
          //}
          
          for ( pix = 0; pix < nvalid; ++pix ) {
            vals [ pix ] += pow ( xx, nx ) * pow ( yy, ny ) * ( (psfimage_[ igroup ])[ co ] )[ pix ];
            //if ( specbin == 0 ) {
              //cerr << "PIX_PCA:     vals[" << pix << "] += " << pow ( xx, nx ) << " * " << pow ( yy, ny ) << " * " << ( (psfimage_[ igroup ])[ co ] )[ pix ] << " => " << vals[pix] << endl;
            //}
          }
        }
        
        //if ( ( spec == 50 ) && ( specbin == 50 ) ) {
          for ( pix = 0; pix < nvalid; ++pix ) {
            cerr << "PIX_PCA:     vals[" << pix << "] = " << vals[ pix ] << endl;
          }
        //}
    
        pix = 0;
    
        for ( size_t imgrow = startY; imgrow <= stopY; ++imgrow ) {
          rowoff = ( imgrow - firstY ) * ( lastX - firstX + 1 );
    
          for ( size_t imgcol = startX; imgcol <= stopX; ++imgcol ) {
      
            datarow = rowoff + ( imgcol - firstX );

            builder ( datarow, b ) += vals[pix];
        
            ++pix;
          }
        }
      
      } else if ( type_ == "GAUSS-HERMITE" ) {
        
        
        
      } else {
        MOAT_THROW( "unknown type in boss psf projection" );
      }
      
    }

    ++complete;
    
    int progfrac = (int) ( 10 * complete / nbins );
    /*
    char msg[256];
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
  
  data.reserve ( nonzeros );
  
  sparse_rowmat :: iterator1 itrow;
  sparse_rowmat :: iterator2 itcol;
  
  /*
  for ( itcol = builder.begin2(); itcol != builder.end2(); ++itcol ) {
    for ( itrow = itcol.begin(); itrow != itcol.end(); ++itrow ) {
      data( itrow.index1(), itcol.index2() ) = (*itrow);
    }
  }
  */
  
  //fprintf ( stderr, "  Memory Remap [          ]\r" );
  
  complete = 0;
  //char msg[256];
  int progfrac;
  lastfrac = 0;
  
  for ( itrow = builder.begin1(); itrow != builder.end1(); ++itrow ) {
    for ( itcol = itrow.begin(); itcol != itrow.end(); ++itcol ) {
      data( itrow.index1(), itcol.index2() ) = (*itcol);
      ++complete;
      //cout << "data[ " << itrow.index1() << ", " << itcol.index2() << " ] = " << (*itcol) << endl;
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


