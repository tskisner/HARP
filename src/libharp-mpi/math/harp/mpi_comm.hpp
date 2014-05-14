// @COPYRIGHT@

#ifndef HARP_MPI_COMM_HPP
#define HARP_MPI_COMM_HPP

#include <boost/iostreams/stream_buffer.hpp>
#include <boost/iostreams/stream.hpp>
#include <boost/iostreams/device/back_inserter.hpp>

#include <boost/archive/portable_oarchive.hpp>
#include <boost/archive/portable_iarchive.hpp>


namespace harp {

  typedef std::vector < char > mpi_comm_buffer_type;

  template < typename T >
  void mpi_comm_pack ( const T & data, mpi_comm_buffer_type & buf ) {
    buf.resize(0);
    buf.clear();

    boost::iostreams::stream < boost::iostreams::back_insert_device < mpi_comm_buffer_type > > output_stream ( buf );
    eos::portable_oarchive oa ( output_stream );

    oa << data;

    output_stream.flush();

    return;
  }

  template < typename T >
  void mpi_comm_unpack ( const mpi_comm_buffer_type & buf, T & data ) {

    boost::iostreams::basic_array_source < char > source ( &buf[0], buf.size() );
    boost::iostreams::stream < boost::iostreams::basic_array_source < char > > input_stream ( source );
    eos::portable_iarchive ia ( input_stream );
  
    ia >> data;

    return;
  }

  template < typename T >
  void mpi_comm_bcast ( boost::mpi::communicator const & comm, T & data, int root ) {

    mpi_comm_buffer_type buf;

    if ( comm.rank() == root ) {
      mpi_comm_pack ( data, buf );
    }

    boost::mpi::broadcast ( comm, buf, root );

    if ( comm.rank() != root ) {
      mpi_comm_unpack ( buf, data );
    }

    return;
  }

  template < typename T >
  void mpi_comm_send ( boost::mpi::communicator const & comm, T & data, int receiver, int tag ) {

    mpi_comm_buffer_type buf;

    mpi_comm_pack ( data, buf );

    comm.send ( receiver, tag, buf );

    return;
  }

  template < typename T >
  void mpi_comm_recv ( boost::mpi::communicator const & comm, T & data, int sender, int tag ) {

    mpi_comm_buffer_type buf;

    comm.recv ( sender, tag, buf );

    mpi_comm_unpack ( buf, data );

    return;
  }


}


#endif
