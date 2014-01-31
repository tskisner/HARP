// @COPYRIGHT@

#ifndef HARP_PLUGIN_HPP
#define HARP_PLUGIN_HPP

// This header file is included by external plugins, and includes all
// serial archive types that we support so that the external plugins 
// can instantiate serialization types with BOOST_CLASS_EXPORT.

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>


#endif
