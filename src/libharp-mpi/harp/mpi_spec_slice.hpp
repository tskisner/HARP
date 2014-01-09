// @COPYRIGHT@

#ifndef HARP_MPI_SPEC_SLICE_HPP
#define HARP_MPI_SPEC_SLICE_HPP


namespace harp {

  /*


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

      std::vector < spec_slice_region > regions ( size_t const & worker ) const {
        if ( worker > nworker_ - 1 ) {
          HARP_THROW( "worker rank is out of range" );
        }
        std::map < size_t, std::vector < spec_slice_region > > :: const_iterator rit = regions_.find ( worker );
        return rit->second;
      }
      
    private :

      size_t nworker_;
      size_t overlap_spec_;
      size_t overlap_lambda_;
      size_t nspec_;
      size_t nlambda_;
      size_t chunk_nspec_;
      size_t chunk_nlambda_;
      std::map < size_t, std::vector < spec_slice_region > > regions_;

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
        return;
      }  
      
  };

  BOOST_SERIALIZATION_SHARED_PTR(spec_slice)
  
  typedef boost::shared_ptr < harp::spec_slice > spec_slice_p;

  */

}


#endif

