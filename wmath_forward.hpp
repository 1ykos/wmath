#ifndef WMATH_FORWARD_H
#define WMATH_FORWARD_H

#include <algorithm>
#include <array>
#include <cassert>
#include <cstdint>
#include <functional>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <limits>
#include <numeric>
#include <random>
#include <vector>
#include <cstring>
#if defined (__has_include) && (__has_include(<x86intrin.h>))
#include <x86intrin.h>
#endif
#include <bitset>

namespace wmath{
  using std::abs;
  using std::accumulate;
  using std::advance;
  using std::array;
  using std::bernoulli_distribution;
  using std::cerr;
  using std::conditional;
  using std::cout;
  using std::distance;
  using std::enable_if;
  using std::endl;
  using std::get;
  using std::index_sequence;
  using std::index_sequence_for;
  using std::initializer_list;
  using std::is_const;
  using std::is_floating_point;
  using std::is_integral;
  using std::is_same;
  using std::is_trivially_copyable;
  using std::is_unsigned;
  using std::isnan;
  using std::istream;
  using std::iter_swap;
  using std::iterator;
  using std::iterator_traits;
  using std::make_unsigned;
  using std::memcpy;
  using std::min_element;
  using std::numeric_limits;
  using std::pair;
  using std::random_access_iterator_tag;
  using std::remove_if;
  using std::result_of;
  using std::setw;
  using std::swap;
  using std::to_string;
  using std::tuple;
  using std::unique;
  using std::vector;
}
#endif // WMATH_FORWARD_H
