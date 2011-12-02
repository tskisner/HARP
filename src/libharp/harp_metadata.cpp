// @COPYRIGHT@

#include <harp_internal.hpp>


using namespace std;
using namespace harp;

image_p harp::meta::parse_image ( boost::ptree & tree ) {
  
  image_p ret = 

}

  boost::ptree serialize_image ( image_p img );

  psf_p parse_psf ( boost::ptree & tree );

  boost::ptree serialize_psf ( psf_p pf );

  spec_p parse_spec ( boost::ptree & tree );

  boost::ptree serialize_spec ( spec_p sp );


  class psfgroup {

    public :
      psfgroup ( );
      psfgroup ( boost::ptree & tree );
      ~psfgroup ( ) { }

      void load ( boost::ptree & tree );
      boost::ptree dump ( );

      psf_p psf ( ) { return handle_; }
      std::list < image_p > images ( ) { return imgs_; }
      
    private :
      psf_p handle_;
      std::list < image_p > imgs_;

  };


