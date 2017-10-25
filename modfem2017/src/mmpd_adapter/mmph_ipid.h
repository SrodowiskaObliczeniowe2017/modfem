#ifndef MMPH_IPID_H_
#define MMPH_IPID_H_


#include "un_ord_map_selector.hpp"

#ifdef __cplusplus
extern "C"{
#endif

/* type IPID - InterProcessorID - main type for parallel modules */
struct mmpt_ipid{

  static const int size_in_ints=2;

  int owner,        /* owner ID (from 1 to nrproc) */
      id_at_owner;

  mmpt_ipid(const int Owner = -1, const int Id_at_owner =- 1)
      : owner(Owner),id_at_owner(Id_at_owner)
  {
  }

  operator const int*() const { return & owner; }

} ;

#ifdef __cplusplus
}
#endif

bool   operator==(const mmpt_ipid& one, const mmpt_ipid& two) ;
bool   operator!=(const mmpt_ipid& one, const mmpt_ipid& two) ;

#ifdef BOOST_NO_CXX11_HDR_UNORDERED_MAP

namespace boost {

template<>
struct hash<mmpt_ipid>
{
    typedef mmpt_ipid argument_type;
    typedef size_t result_type;

    result_type operator () (const argument_type& x) const {
        return boost::hash<int64_t>()( *reinterpret_cast<const int64_t*>(& x) );
    }
};

}

#else

/// Hashes for std::unordered_map
namespace std {

template<>
struct hash<mmpt_ipid>
{
    typedef mmpt_ipid argument_type;
    typedef size_t result_type;

    result_type operator () (const argument_type& x) const {
        return std::hash<int64_t>()( *reinterpret_cast<const int64_t*>(& x) );
    }
};

template<>
    struct __is_fast_hash<hash<mmpt_ipid>> : public std::false_type
    { };

}

#endif

#endif //MMPH_IPID_H_
