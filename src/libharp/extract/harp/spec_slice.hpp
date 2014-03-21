// @COPYRIGHT@

#ifndef HARP_SPEC_SLICE_HPP
#define HARP_SPEC_SLICE_HPP


namespace harp {


  // class for a single slice point

  class spec_slice_region {

    friend class boost::serialization::access;

    public :

      spec_slice_region ( ) :
        overlap_spec( 0 ),
        overlap_lambda( 0 ),
        first_spec( 0 ),
        first_lambda( 0 ),
        first_good_spec( 0 ),
        first_good_lambda( 0 ),
        n_spec( 0 ),
        n_lambda( 0 ),
        n_good_spec( 0 ),
        n_good_lambda( 0 ) { }

      ~spec_slice_region ( ) { }

      spec_slice_region ( spec_slice_region const & orig ) :
        overlap_spec( orig.overlap_spec ),
        overlap_lambda( orig.overlap_lambda ),
        first_spec( orig.first_spec ),
        first_lambda( orig.first_lambda ),
        first_good_spec( orig.first_good_spec ),
        first_good_lambda( orig.first_good_lambda ),
        n_spec( orig.n_spec ),
        n_lambda( orig.n_lambda ),
        n_good_spec( orig.n_good_spec ),
        n_good_lambda( orig.n_good_lambda ) { }

      spec_slice_region & operator= ( spec_slice_region const & rhs ) {
        if ( &rhs != this ) {
          overlap_spec = rhs.overlap_spec;
          overlap_lambda = rhs.overlap_lambda;
          first_spec = rhs.first_spec;
          first_lambda = rhs.first_lambda;
          first_good_spec = rhs.first_good_spec;
          first_good_lambda = rhs.first_good_lambda;
          n_spec = rhs.n_spec;
          n_lambda = rhs.n_lambda;
          n_good_spec = rhs.n_good_spec;
          n_good_lambda = rhs.n_good_lambda;
        }
        return *this;
      }

      size_t overlap_spec;
      size_t overlap_lambda;
      size_t first_spec;
      size_t first_lambda;
      size_t first_good_spec;
      size_t first_good_lambda;
      size_t n_spec;
      size_t n_lambda;
      size_t n_good_spec;
      size_t n_good_lambda;

    private :

      template < class Archive >
      void serialize ( Archive & ar, const unsigned int version ) {
        ar & BOOST_SERIALIZATION_NVP(overlap_spec);
        ar & BOOST_SERIALIZATION_NVP(overlap_lambda);
        ar & BOOST_SERIALIZATION_NVP(first_spec);
        ar & BOOST_SERIALIZATION_NVP(first_lambda);
        ar & BOOST_SERIALIZATION_NVP(first_good_spec);
        ar & BOOST_SERIALIZATION_NVP(first_good_lambda);
        ar & BOOST_SERIALIZATION_NVP(n_spec);
        ar & BOOST_SERIALIZATION_NVP(n_lambda);
        ar & BOOST_SERIALIZATION_NVP(n_good_spec);
        ar & BOOST_SERIALIZATION_NVP(n_good_lambda);
        return;
      }

  };

  BOOST_SERIALIZATION_SHARED_PTR(spec_slice_region)


  class spec_slice : public boost::enable_shared_from_this < spec_slice > {

    friend class boost::serialization::access;
    
    public :

      spec_slice ( ) :
        nworker_(0),
        overlap_spec_(0),
        overlap_lambda_(0),
        nspec_(0),
        nlambda_(0),
        chunk_nspec_(0),
        chunk_nlambda_(0) { }

      spec_slice ( size_t nworker, size_t nspec, size_t nlambda, size_t chunk_nspec, size_t chunk_nlambda, size_t overlap_spec, size_t overlap_lambda );
      
      ~spec_slice ( ) { }

      spec_slice_region const & full_region ( ) const { return full_region_; }

      std::vector < spec_slice_region > regions ( size_t const & worker ) const;
      
    private :

      size_t nworker_;
      size_t overlap_spec_;
      size_t overlap_lambda_;
      size_t nspec_;
      size_t nlambda_;
      size_t chunk_nspec_;
      size_t chunk_nlambda_;
      std::map < size_t, std::vector < spec_slice_region > > regions_;
      spec_slice_region full_region_;

      void calc ( size_t n, size_t chunk, size_t overlap, std::vector < size_t > & start, std::vector < size_t > & stop, std::vector < size_t > & good_start, std::vector < size_t > & good_stop );

      template < class Archive >
      void serialize ( Archive & ar, const unsigned int version ) {
        ar & BOOST_SERIALIZATION_NVP(nworker_);
        ar & BOOST_SERIALIZATION_NVP(overlap_spec_);
        ar & BOOST_SERIALIZATION_NVP(overlap_lambda_);
        ar & BOOST_SERIALIZATION_NVP(nspec_);
        ar & BOOST_SERIALIZATION_NVP(nlambda_);
        ar & BOOST_SERIALIZATION_NVP(chunk_nspec_);
        ar & BOOST_SERIALIZATION_NVP(chunk_nlambda_);
        ar & BOOST_SERIALIZATION_NVP(regions_);
        ar & BOOST_SERIALIZATION_NVP(full_region_);
        return;
      }  
      
  };

  BOOST_SERIALIZATION_SHARED_PTR(spec_slice)
  
  typedef boost::shared_ptr < harp::spec_slice > spec_slice_p;

}


#endif

