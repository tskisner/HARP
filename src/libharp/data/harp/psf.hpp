/*
  High Performance Astrophysical Reconstruction and Processing (HARP)
  (c) 2014-2015, The Regents of the University of California, 
  through Lawrence Berkeley National Laboratory.  See top
  level LICENSE file for details.
*/

#ifndef HARP_PSF_HPP
#define HARP_PSF_HPP


namespace harp {

  class psf : public boost::enable_shared_from_this < psf > {

    friend class boost::serialization::access;
    
    public :

      psf ( ) {
        type_ = "";
      }

      psf ( std::string const & type, boost::property_tree::ptree const & props );

      virtual ~psf ( ) { }

      virtual size_t n_spec ( ) const {
        HARP_THROW( "fell through to virtual method" );
        return 0;
      }
      
      virtual size_t n_lambda ( ) const {
        HARP_THROW( "fell through to virtual method" );
        return 0;
      }

      virtual size_t img_rows ( ) const {
        HARP_THROW( "fell through to virtual method" );
        return 0;
      }

      virtual size_t img_cols ( ) const {
        HARP_THROW( "fell through to virtual method" );
        return 0;
      }

      virtual vector_double lambda ( ) const {
        HARP_THROW( "fell through to virtual method" );
        return vector_double();
      }

      virtual void extent ( size_t spec, size_t lambda, size_t & x_offset, size_t & y_offset, size_t & n_x, size_t & n_y ) const {
        HARP_THROW( "fell through to virtual method" );
        return;
      }

      virtual void response ( size_t spec, size_t lambda, size_t & x_offset, size_t & y_offset, matrix_double & patch ) const {
        HARP_THROW( "fell through to virtual method" );
        return;
      }

      virtual size_t response_nnz_estimate ( ) const {
        HARP_THROW( "fell through to virtual method" );
        return 0;
      }

      // The base class provides a default implementation of these 3 functions, so that derived classes
      // only need to implement the response() and extent() methods.

      virtual void extent_multi ( std::map < size_t, std::set < size_t > > const & speclambda, std::vector < size_t > & x_offset, std::vector < size_t > & y_offset, std::vector < size_t > & n_x, std::vector < size_t > & n_y ) const;
      
      virtual void project ( std::map < size_t, std::set < size_t > > const & speclambda, matrix_double & A ) const;

      virtual void project_transpose ( std::map < size_t, std::set < size_t > > const & speclambda, matrix_double_sparse & AT ) const;

      // These are convenience functions if you want the whole projection matrix

      void project ( matrix_double & A ) const;

      void project_transpose ( matrix_double_sparse & AT ) const;


      // convenience function
      static size_t total_bins ( std::map < size_t, std::set < size_t > > const & speclambda );
      
      std::string type ( ) const;
      
      template < class T >
      boost::shared_ptr < T > shared_ref ( ) {
        return boost::dynamic_pointer_cast < T, psf > ( shared_from_this ( ) );
      }
      
      template < class T >
      boost::weak_ptr < T > weak_ref ( ) {
        boost::shared_ptr < T > temp = boost::dynamic_pointer_cast < T, psf > ( shared_from_this() );
        return boost::weak_ptr < T > ( temp );
      }
      
    private :

      template < class Archive >
      void serialize ( Archive & ar, const unsigned int version ) {
        ar & BOOST_SERIALIZATION_NVP(type_);
        ar & BOOST_SERIALIZATION_NVP(props_);
        return;
      }
    
      std::string type_;
      boost::property_tree::ptree props_;
      
      
  };

  BOOST_SERIALIZATION_ASSUME_ABSTRACT(psf)

  BOOST_SERIALIZATION_SHARED_PTR(psf)
  
  typedef boost::shared_ptr < harp::psf > psf_p;
  typedef boost::weak_ptr < harp::psf > psf_wp;

}


#endif

