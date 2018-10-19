#ifndef WMATH_HASH_H
#define WMATH_HASH_H

#include <cstdint>
#include <limits>
#include <limits.h>
#include <numeric>
#include "wmath_forward.hpp"
#include "wmath_primes.hpp"
#include "wmath_bits.hpp"

namespace wmath{
  
  using std::enable_if;
  using std::false_type;
  using std::true_type;
  using std::is_trivially_copyable;
  using std::integral_constant;
  using std::is_fundamental;
  using std::numeric_limits;
  
  template<typename,typename=void>
  struct is_injective : false_type {};

  template<typename T>
  struct is_injective<T,typename enable_if<T::is_injective::value>::type>
  : true_type {};
  //: false_type {};

  template<typename,typename=void>
  struct has_std_hash : false_type {};

  template<typename T>
  struct has_std_hash<T,decltype(std::hash<T>()(std::declval<T>()),void())>
  : true_type {};

  template<typename T>
  typename
  enable_if<has_std_hash<T>::value&&(!is_fundamental<T>::value),size_t>::type
  constexpr hash(const T& v){
    return std::hash<T>()(v);
  }

  size_t constexpr hash(const size_t& v){
    return v;
  }
  
  size_t constexpr hash(const size_t& v0,const size_t& v1) {
    // This hash_combine is slightly better than what boost offers, but it could
    // be better too by using multiplication in a galois field like ghash
    const size_t i = 1+(numeric_limits<size_t>::digits*144+116)/233;
    const size_t m = numeric_limits<size_t>::digits-1;
    const size_t c = i&m;
    const size_t n = v0+v1;
    return ((n<<c)|(n>>((-c)&m)))^(v0-v1);
  }

  template<class K> typename enable_if<is_fundamental<K>::value,size_t>::type
  constexpr hash(const K& k){
    size_t h(k);
    const size_t n = sizeof(K)/sizeof(size_t);
    for (size_t i=sizeof(size_t);i<sizeof(K);i+=sizeof(size_t))
      h = hash(h,size_t(k>>(i*CHAR_BIT)));
    return h;
  }
 
  template <typename T,typename... Rest>
  size_t constexpr hash(const T& v,Rest... rest);
  
  template<typename T,size_t... I>
  size_t constexpr hash_tuple_impl(const T& t, index_sequence<I...>){
    return hash(std::get<I>(t)...);
  }

  template<typename... Ts>
  size_t constexpr hash(const tuple<Ts...>& t){
    return hash_tuple_impl(t,index_sequence_for<Ts...>{});
  }
  
  template <typename T, typename... Rest>
  size_t constexpr hash(const T& v, Rest... rest) {
    return hash(hash(v),hash(rest...));
  }

  template<class K,class enable = void>
  struct hash_functor{
    typedef typename false_type::type is_injective;
    size_t operator()(const K& k) const {
      return hash(k);
    }
  };
  
  template<class K>
  struct hash_functor<
    K,
    typename enable_if<is_fundamental<K>::value,void>::type
  >{
    // if size_t has at least as many digits as the hashed type the hash can
    // be injective and I will make it so
    typedef typename integral_constant<bool,sizeof(K)<=sizeof(size_t)>::type
      is_injective;
    size_t constexpr operator()(const K& k) const {
      return hash(k);
    }
  };

  /*
  template<> struct
  hash<typename enable_if<sizeof(bool)<=sizeof(size_t),bool>::type>
  {
    typedef typename true_type is_injective;
    uint8_t constexpr operator()(const bool& k) const {
      return k;
    }
  };
  
  template<> struct
  hash<typename enable_if<(sizeof(bool)>sizeof(size_t)),bool>::type>
  {
    uint8_t constexpr operator()(const bool& k) const {
      size_t h(k);
      const size_t n = sizeof(bool)/sizeof(size_t);
      for (size_t i=0;i!=n;++i){
        h = hash_combine(h,size_t(k>>(i*CHAR_BIT)));
      }
      return h;
    }
  };

  template<> struct
  hash<typename enable_if<(sizeof(char)<=sizeof(size_t)),char>::type>
  {
    typedef typename true_type is_injective;
    uint8_t constexpr operator()(const char& k) const {
      return k;
    }
  };
  
  template<> struct
  hash<typename enable_if<(sizeof(char)>sizeof(size_t)),char>::type>
  {
    uint8_t constexpr operator()(const char& k) const {
      return k;
    }
  };

  template<> struct
  hash<typename enable_if<(sizeof(int8_t)<=sizeof(size_t)),int8_t>::type>
  {
    typedef typename true_type is_injective;
    uint8_t constexpr operator()(const signed char& k) const {
      return k;
    }
  };
  
  template<> struct
  hash<typename enable_if<(sizeof(int8_t)>sizeof(size_t)),int8_t>::type>
  {
    uint8_t constexpr operator()(const signed char& k) const {
      return k;
    }
  };

  template<> struct
  hash<typename enable_if<(sizeof(uint8_t)<=sizeof(size_t)),uint8_t>::type>
  {
    typedef typename true_type is_injective;
    uint8_t constexpr operator()(const uint8_t& k) const {
      return k;
    }
  };
  
  template<> struct
  hash<typename enable_if<(sizeof(uint8_t)>sizeof(size_t)),uint8_t>::type>
  {
    uint8_t constexpr operator()(const uint8_t& k) const {
      return k;
    }
  };

  template<> struct
  hash<typename enable_if<(sizeof(char16_t)<=sizeof(size_t)),char16_t>::type>
  {
    typedef typename true_type is_injective;
    uint16_t constexpr operator()(const char16_t& k) const {
      return k;
    }
  };
  
  template<> struct
  hash<typename enable_if<(sizeof(char16_t)>sizeof(size_t)),char16_t>::type>
  {
    uint16_t constexpr operator()(const char16_t& k) const {
      return k;
    }
  };

  template<> struct
  hash<typename enable_if<(sizeof(char32_t)<=sizeof(size_t)),char32_t>::type>
  {
    typedef typename true_type is_injective;
    uint32_t constexpr operator()(const char32_t& k) const {
      return k;
    }
  };
  
  template<> struct
  hash<typename enable_if<(sizeof(char32_t)>sizeof(size_t)),char32_t>::type>
  {
    uint32_t constexpr operator()(const char32_t& k) const {
      return k;
    }
  };

  template<> struct
  hash<typename enable_if<(sizeof(wchar_t)<sizeof(size_t)),wchar_t>::type>
  {
    typedef typename true_type is_injective;
    constexpr typename std::make_unsigned<T>::type operator()(const wchar_t& k)
    const{
      return k;
    }
  };
  
  template<> struct
  hash<typename enable_if<(sizeof(wchar_t)>=sizeof(size_t)),wchar_t>::type>
  {
    constexpr typename std::make_unsigned<T>::type operator()(const wchar_t& k)
    const{
      return k;
    }
  };

  template<> struct
  hash<typename enable_if<(sizeof(int16_t)<sizeof(size_t)),int16_t>::type>
  {
    typedef typename true_type is_injective;
    uint16_t constexpr operator()(const int16_t& k) const {
      return k;
    }
  };
  
  template<> struct
  hash<typename enable_if<(sizeof(int16_t)>=sizeof(size_t)),int16_t>::type>
  {
    uint16_t constexpr operator()(const int16_t& k) const {
      return k;
    }
  };
  
  template<> struct
  hash<typename enable_if<(sizeof(uint16_t)<sizeof(size_t)),uint16_t>::type>
  {
    typedef typename true_type is_injective;
    uint16_t constexpr operator()(const uint16_t& k) const {
      return k;
    }
  };
  
  template<> struct
  hash<typename enable_if<(sizeof(uint16_t)>=sizeof(size_t)),uint16_t>::type>
  {
    uint16_t constexpr operator()(const uint16_t& k) const {
      return k;
    }
  };

  template<> struct
  hash<typename enable_if<(sizeof(int32_t)<sizeof(size_t)),int32_t>::type>
  {
    typedef typename true_type is_injective;
    uint32_t constexpr operator()(const int32_t& k) const {
      return k;
    }
  };
  
  template<> struct
  hash<typename enable_if<(sizeof(int32_t)>=sizeof(size_t)),int32_t>::type>
  {
    uint32_t constexpr operator()(const int32_t& k) const {
      return k;
    }
  };

  template<> struct
  hash<typename enable_if<(sizeof(uint32_t)<sizeof(size_t)),uint32_t>::type>
    typedef typename true_type is_injective;
    uint32_t constexpr operator()(const uint32_t& k) const {
      return k;
    }
  };
  
  template<> struct
  hash<typename enable_if<(sizeof(uint32_t)>=sizeof(size_t)),uint32_t>::type>
  {
    uint32_t constexpr operator()(const uint32_t& k) const {
      return k;
    }
  };

  template<> struct
  hash<typename enable_if<(sizeof(int64_t)<=sizeof(size_t)),int64_t>::type>
  {
    typedef typename true_type is_injective;
    uint64_t constexpr operator()(const int64_t& k) const {
      return k;
    }
  };
  
  template<> struct
  hash<typename enable_if<(sizeof(int64_t)>sizeof(size_t)),int64_t>::type>
  {
    uint64_t constexpr operator()(const int64_t& k) const {
      return k;
    }
  };

  template<> struct
  hash<typename enable_if<(sizeof(float)<=sizeof(size_t)),float>::type>
  {
    typedef typename true_type is_injective;
    uint32_t operator()(const float& k) const{
      return *reinterpret_cast<uint32_t*>(&k);
    }
  };
  
  template<> struct
  hash<typename enable_if<(sizeof(float)>sizeof(size_t)),float>::type>
  {
    uint32_t operator()(const float& k) const{
      return *reinterpret_cast<uint32_t*>(&k);
    }
  };

  template<> struct
  hash<typename enable_if<(sizeof(double)<=sizeof(size_t)),double>::type>
    typedef typename true_type is_injective;
    uint64_t operator()(const double& k) const{
      return *reinterpret_cast<uint32_t*>(&k);
    }
  };

  template<> struct
  hash<typename enable_if<(sizeof(double)>=sizeof(size_t)),double>::type>
    uint64_t operator()(const double& k) const{
      return *reinterpret_cast<uint32_t*>(&k);
    }
  };

  // template<> struct hash<long double>; TODO
  
  template<class T>
  hash<T*>
  {
    typedef typename true_type is_injective;
    size_t constexpr operator()(const T*& k){
      return k;
    }
  };
  */
}
#endif // WMATH_HASH_H
