/*
  High Performance Astrophysical Reconstruction and Processing (HARP)
  (c) 2014-2015, The Regents of the University of California, 
  through Lawrence Berkeley National Laboratory.  See top
  level LICENSE file for details.
*/

#ifndef HARP_MPI_IMAGE_HPP
#define HARP_MPI_IMAGE_HPP


namespace harp {

  class mpi_image {

    public :

      mpi_image ( boost::mpi::communicator const & comm, std::string const & type, boost::property_tree::ptree const & props );
      
      ~mpi_image ( ) { }

      boost::mpi::communicator comm ( ) { return comm_; }
      
      size_t n_rows ( ) const;
      
      size_t n_cols ( ) const;
      
      void values ( vector_double & data ) const;

      void inv_variance ( vector_double & invvar ) const;
      
      void values ( matrix_double & data ) const;

      void inv_variance ( matrix_double & invvar ) const;

      void values ( elem_matrix_local & data ) const;

      void inv_variance ( elem_matrix_local & invvar ) const;

      void mask ( vector_mask & msk ) const;

      void mask ( matrix_mask & msk ) const;

      std::string type ( ) const;

      image_p local ();
      
    private :
    
      boost::mpi::communicator comm_;
      image_p local_;
      
  };
  
  typedef boost::shared_ptr < harp::mpi_image > mpi_image_p;
  typedef boost::weak_ptr < harp::mpi_image > mpi_image_wp;

}


#endif

