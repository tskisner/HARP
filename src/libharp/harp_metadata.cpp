// @COPYRIGHT@

#include <harp_internal.hpp>


using namespace std;
using namespace harp;


harp::meta::psfgroup


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


