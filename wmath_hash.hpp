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
  
  template <typename T>
  typename std::enable_if<std::is_unsigned<T>::value,T>::type
  constexpr distribute(const T& a); // mix the hash value good, clmul_circ
                                    // with odious integer is suitable

  uint8_t  const inline distribute(const uint8_t& a){
    return clmul_circ(a,uint8_t(0b01000101u));
    return (a+111)*97;
  } 
  uint16_t const inline distribute(const uint16_t& a){
    return clmul_mod(a,uint16_t(0b01000101'10100101u));
    return (a+36690)*43581;
  }
  uint32_t const inline distribute(const uint32_t& a){
    const uint32_t c0 = 3107070805ul;
    const uint32_t c1 = 3061963241ul;
    return rol((a+c1)*c0,16)*c1;
    return clmul_mod(uint32_t(a+3061963241ul),uint32_t(3107070805ul));
    return (a^(a>>16))*3061963241ul;
  }
  uint64_t const inline distribute(const uint64_t& a){
    const uint64_t c0 = 16123805160827025777ull;
    const uint64_t c1 =  3619632413061963241ull;
    return rol((a+c0)*c1,32)*c0;
    return clmul_mod(a,uint64_t(16123805160827025777ull))
          *16123805160827025777ull;
    return (a^(a>>32))*16123805160827025777ull;
    return (a+3619632413061963241ull)*16123805160827025777ull;
  }
  
  
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

  size_t constexpr hash(const size_t& v0,const size_t& v1) {
    // This hash_combine is slightly better than what boost offers, but it could
    // be better too by using multiplication in a galois field like ghash
    // or just any multiplication at all
    const size_t i = 1+(numeric_limits<size_t>::digits*144+116)/233;
    const size_t m = numeric_limits<size_t>::digits-1;
    const size_t c = i&m;
    const size_t n = v0+v1;
    return ((n<<c)|(n>>((-c)&m)))^(v0-v1);
  }

  template<class K>
  typename enable_if<(sizeof(K)>sizeof(size_t)
                   &&is_fundamental<K>::value),size_t>::type
  constexpr hash(const K& k){
    size_t h(k);
    const size_t n = sizeof(K)/sizeof(size_t);
    for (size_t i=sizeof(size_t);i<sizeof(K);i+=sizeof(size_t))
      h = hash(h,size_t(k>>(i*CHAR_BIT)));
    return h;
  }

  uint8_t constexpr hash(const uint8_t& v){
    return v;
  }

  uint8_t constexpr hash(const int8_t& v){
    return v;
  }

  uint16_t constexpr hash(const uint16_t& v){
    return v;
  }

  uint16_t constexpr hash(const int16_t& v){
    return v;
  }

  uint32_t constexpr hash(const uint32_t& v){
    return v;
  }

  uint32_t constexpr hash(const int32_t& v){
    return v;
  }

  uint64_t constexpr hash(const uint64_t& v){
    return v;
  }

  uint64_t constexpr hash(const int64_t& v){
    return v;
  }
 
  template <typename T,typename... Rest>
  size_t constexpr hash(const T& v,Rest... rest);

  uint16_t constexpr hash(const uint8_t& v0,const uint8_t& v1){
    return (uint16_t(v0)<<8)^uint16_t(v1);
  }

  uint32_t constexpr hash(const uint16_t& v0,const uint16_t& v1){
    return (uint32_t(v0)<<16)^(uint32_t(v1));
  }
  
  uint64_t constexpr hash(const uint32_t& v0,const uint32_t& v1){
    return (uint64_t(v0)<<32)^(uint64_t(v1));
  }
  
  uint64_t constexpr hash(const uint64_t& v0,const uint32_t& v1){
    return hash(uint64_t(v0),uint64_t(v1));
  }
  
  uint64_t constexpr hash(const uint32_t& v0,const uint64_t& v1){
    return hash(uint64_t(v0),uint64_t(v1));
  }

  template<typename T,size_t... I>
  size_t constexpr hash_tuple_impl(const T& t, index_sequence<I...>){
    return hash(std::get<I>(t)...);
  }

  template<typename... Ts>
  size_t constexpr hash(const tuple<Ts...>& t){
    return hash_tuple_impl(t,index_sequence_for<Ts...>{});
  }

  template<typename T,size_t n>
  size_t constexpr hash(const array<T,n> a){
    size_t h(0);
    for (size_t i=0;i!=n;++i) {
      if constexpr(sizeof(T)<=sizeof(size_t))
        if (i%(sizeof(size_t)/sizeof(T))==0) h = distribute(h); 
      h = rol(h,sizeof(T)*CHAR_BIT)^hash(a[i]);
      if constexpr(sizeof(T)>sizeof(size_t)) h = distribute(h);
    }
    return h;
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
    //typedef typename integral_constant<CHAR_BIT*sizeof(K)>::type bits; 
    auto constexpr operator()(const K& k) const {
      return hash(k);
    }
  };
  
  template<typename T,size_t n>
  struct hash_functor<
    array<T,n>,void>
  {
    // if size_t has at least as many digits as the hashed type the hash can
    // be injective and I will make it so
    typedef typename integral_constant<bool,n*sizeof(T)<=sizeof(size_t)>::type
      is_injective;
    //typedef typename integral_constant<CHAR_BIT*sizeof(K)>::type bits; 
    auto constexpr operator()(const array<T,n>& k) const {
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
