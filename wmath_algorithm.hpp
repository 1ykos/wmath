#include "wmath.hpp"

namespace wmath{
  /* interpolation search in half open interval [begin,end)
   * elements in range are assumed to be in order according to comp and 
   * the application of the ordering function yields uniformly distributed 
   * integers
   */
  template<class random_access_iterator,class comparator>
  auto interpolation_search(
      random_access_iterator begin,
      random_access_iterator end,
      const ordering& order, 
      const comparator& comp
      ) {
    while (true) {
      if (begin==end) return begin;

    }
  }
};
