#ifndef _UN_ORD_MAP_SELECTOR_HPP_
#define _UN_ORD_MAP_SELECTOR_HPP_


#include <boost/config.hpp> //< for BOOST_NO_CXX11_HDR_UNORDERED_MAP

#ifdef BOOST_NO_CXX11_HDR_UNORDERED_MAP
#include <boost/unordered/unordered_map.hpp>

using boost::unordered::unordered_map;

#else

#include <unordered_map>
using std::unordered_map;

#endif


#endif // _UN_ORD_MAP_SELECTOR_HPP_
