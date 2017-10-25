
#include <inttypes.h>
#include "mmph_ipid.h"

bool   operator==(const mmpt_ipid& one, const mmpt_ipid& two) {
    return ( *reinterpret_cast<const int64_t*>(& one) == *reinterpret_cast<const int64_t*>(& two));
}

bool   operator!=(const mmpt_ipid& one, const mmpt_ipid& two) {
    return ( *reinterpret_cast<const int64_t*>(& one) != *reinterpret_cast<const int64_t*>(& two));
}
