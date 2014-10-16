/*
  High Performance Astrophysical Reconstruction and Processing (HARP)
  (c) 2014-2015, The Regents of the University of California, 
  through Lawrence Berkeley National Laboratory.  See top
  level LICENSE file for details.
*/

#ifndef HARP_MPI_METADATA_HPP
#define HARP_MPI_METADATA_HPP

#include <list>


namespace harp {

  mpi_spec_p mpi_load_spec ( boost::mpi::communicator const & comm, boost::property_tree::ptree const & tree );

  mpi_image_p mpi_load_image ( boost::mpi::communicator const & comm, boost::property_tree::ptree const & tree );

  mpi_psf_p mpi_load_psf ( boost::mpi::communicator const & comm, boost::property_tree::ptree const & tree );

  mpi_targets_p mpi_load_targets ( boost::mpi::communicator const & comm, boost::property_tree::ptree const & tree );


  class mpi_group {

    public :
      
      mpi_group ( boost::mpi::communicator const & comm, boost::property_tree::ptree const & tree );

      void load ( boost::mpi::communicator const & comm, boost::property_tree::ptree const & tree );

      mpi_psf_p psf ( ) { return handle_; }
      std::list < mpi_image_p > images ( ) { return imgs_; }
      
    private :
      mpi_psf_p handle_;
      std::list < mpi_image_p > imgs_;

  };

  typedef boost::shared_ptr < mpi_group > mpi_group_p;

  std::list < mpi_group_p > mpi_load_groups ( boost::mpi::communicator const & comm, std::string const & path );

}

#endif

