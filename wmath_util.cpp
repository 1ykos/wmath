#ifndef WMATH_UTIL_H
#define WMATH_UTIL_H
#include "wmath_forward.hpp"
#include "wmath_bits.hpp"
namespace wmath{
  // random_order_iterator iterates range in random order
  template<class iterator,class generator>
  class random_order_iterator{
    private:
    const iterator begin;
    const iterator end;
    const size_t n;
    const size_t m;
    size_t i = 0;
    size_t a;
    size_t b;
    public:
    random_order_iterator(
        const iterator begin,
        const iterator end,
        generator& rng
        ): begin(begin),end(end),n(end-begin),m(roundup_2(n)) {
      uniform_int_distribution<size_t> distr(m);
      a = 2*distr(rng)+1;
      b = distr(rng);
    }
    random_order_iterator(
        const iterator begin,
        const iterator end,
        const size_t& n,
        const size_t& m,
        const size_t& i,
        const size_t& a,
        const size_t& b
        ): begin(begin),end(end),n(n),m(m),i(i),a(a),b(b) {
      uniform_int_distribution<size_t> distr(m);
      a = 2*distr(rng)+1;
      b = distr(rng);
    }
  }
  random_order_iterator& operator++(){ // prefix
    ++i;
    return *this;
  }
  random_order_iterator operator++(int){ // postfix
    random_order_iterator pre = *this;
    ++i;
    return pre;
  }
  random_order_iterator& operator--(){ // prefix
    --i;
    return *this;
  }
  random_order_iterator operator--(int){ // postfix
    random_order_iterator pre = *this;
    --i;
    return pre;
  }
}
#endif // WMATH_UTIL_H
