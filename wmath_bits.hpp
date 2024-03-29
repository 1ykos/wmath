#ifndef WMATH_BITS_H
#define WMATH_BITS_H
#include "wmath_forward.hpp"
#include "wmath_primes.hpp"

namespace wmath{
  
  template <typename T>
  constexpr size_t digits(const T& n=0){
    return std::numeric_limits<T>::digits;
  }
  
  template<typename T>
  typename std::enable_if<std::is_unsigned<T>::value,tuple<T,T>>::type
  constexpr long_mul0(const T a, const T b){
    const T N  = digits<T>()/2;
    const T t0 = (a>>N)*(b>>N);
    const T t1 = ((a<<N)>>N)*(b>>N);
    const T t2 = (a>>N)*((b<<N)>>N);
    const T t3 = ((a<<N)>>N)*((b<<N)>>N);
    const T t4 = t3+(t1<<N);
    const T r1 = t4+(t2<<N);
    const T r0 = (r1<t4)+(t4<t3)+(t1>>N)+(t2>>N)+t0;
    return {r0,r1};
  }
  
  //template<typename T>
  //constexpr long_mul(const T a, const T b);
  
  // calculate a * b = r0r1
  template<typename T>
  typename std::enable_if<std::is_unsigned<T>::value,tuple<T,T>>::type
  constexpr long_mul(const T a, const T b){
    const T N  = digits<T>()/2;
    const T t0 = (a>>N)*(b>>N);
    const T t1 = ((a<<N)>>N)*(b>>N);
    const T t2 = (a>>N)*((b<<N)>>N);
    const T t3 = ((a<<N)>>N)*((b<<N)>>N);
    const T t4 = t3+(t1<<N);
    const T r1 = t4+(t2<<N);
    const T r0 = (r1<t4)+(t4<t3)+(t1>>N)+(t2>>N)+t0;
    return {r0,r1};
  }
 
#ifdef __SIZEOF_INT128__
  //template<>
  tuple<uint64_t,uint64_t>
  constexpr long_mul(const uint64_t a, const uint64_t b){
    unsigned __int128 r = ((unsigned __int128)(a))*((unsigned __int128)(b));
    return {r>>64,r};
  }

  //template<>
  tuple<__uint128_t,__uint128_t>
  constexpr long_mul(const __uint128_t a, const __uint128_t b){
    const int         N  = 64;
    const __uint128_t t0 = (a>>N)*(b>>N);
    const __uint128_t t1 = ((a<<N)>>N)*(b>>N);
    const __uint128_t t2 = (a>>N)*((b<<N)>>N);
    const __uint128_t t3 = ((a<<N)>>N)*((b<<N)>>N);
    const __uint128_t t4 = t3+(t1<<N);
    const __uint128_t r1 = t4+(t2<<N);
    const __uint128_t r0 = (r1<t4)+(t4<t3)+(t1>>N)+(t2>>N)+t0;
    return {r0,r1};
  }
#endif

  template<typename T,size_t n,size_t m>
  typename std::enable_if<std::is_unsigned<T>::value,array<T,max(n,m)>>::type
  constexpr add(
      const array<T,n> a,
      const array<T,m> b
      ) {
    array<T,max(n,m)> r{};
    bool carry = false;
    for (size_t i=0;i!=min(n,m);++i) {
      r[i]=carry+a[i];
      carry=(r[i]<a[i]);
      r[i]+=b[i];
      carry|=r[i]<a[i];
    }
    if (n<m) {
      for (size_t i=n;i!=m;++i) {
        r[i]=carry+b[i];
        carry=(r[i]<b[i]);
      }
    }
    if (m<n) {
      for (size_t i=m;i!=n;++i) {
        r[i]=carry+a[i];
        carry=(r[i]<a[i]);
      }
    }
    return r;
  }
  
  // compute a-b
  template<typename T,size_t n,size_t m>
  typename std::enable_if<std::is_unsigned<T>::value,array<T,max(n,m)>>::type
  constexpr sub(
      const array<T,n> a,
      const array<T,m> b
      ) {
    array<T,max(n,m)> r{};
    bool carry = false;
    for (size_t i=0;i!=min(n,m);++i) {
      r[i]=a[i]-carry;
      carry=(r[i]>a[i]);
      r[i]-=b[i];
      carry|=r[i]>a[i];
    }
    if (n<m) {
      for (size_t i=n;i!=m;++i) {
        r[i]=T(0)-carry;
        r[i]-=b[i];
        carry|=(r[i]>b[i]);
      }
    }
    if (m<n) {
      for (size_t i=m;i!=n;++i) {
        r[i]=a[i];
        r[i]-=carry;
        carry=(r[i]>a[i]);
      }
    }
    return r;
  }
  
  template<class it0,class it1>
  void const add_with_carry(
      const it0 a_begin,
      const it0 a_end,
            it1 o_begin,
            it1 o_end
      ) {
    bool carry = false;
    auto a = a_begin;
    for (auto it=o_begin;it!=o_end;++it) {
      *it+=carry;
      carry=(*it)<carry;
      if (a!=a_end) {
        *it+=*a;
        carry|=((*it)<(*a));
        ++a;
      } else {
        if (!carry) break;
      }
    }
  }
  
  template<class it0,class it1>
  void const sub_with_carry(
      const it0 a_begin,
      const it0 a_end,
            it1 o_begin,
            it1 o_end
      ) {
    bool carry = false;
    auto a = a_begin;
    for (auto it=o_begin;it!=o_end;++it) {
      const auto t = *it;
      *it-=carry;
      carry=(*it)>t;
      if (a!=a_end) {
        const auto t = *it;
        *it-=*a;
        carry|=((*it)>t);
        ++a;
      } else {
        if (!carry) break;
      }
    }
  }

  template<typename T,size_t n,size_t m>
  typename std::enable_if<std::is_unsigned<T>::value,array<T,n+m>>::type
  constexpr long_mul(
      const array<T,n> a,
      const array<T,m> b
      ) {
    array<T,n+m> r{};
    for (size_t j=0;j!=m;++j) {
      for (size_t i=0;i!=n;++i) {
        const auto t = long_mul(a[i],b[j]);
        //cout << r[0] << " " << r[1] << " " << r[2] << " " << r[3] << endl;
        const array<T,2> a = {get<1>(t),get<0>(t)}; 
        add_with_carry(a.begin(),a.end(),r.begin()+i+j,r.end());
      }
    }
    return r;
  }
  
  // all the extra work means this is slower even though it saves one
  // multiplication... it's about 30% slower than long_mul
  template<typename T>
  tuple<T,T> constexpr karatsuba(const T a, const T b) {
    const T N  = digits<T>()/2;
    const T x1 = a>>N;
    const T x2 = a&((~T(0))>>N);
    const T y1 = b>>N;
    const T y2 = b&((~T(0))>>N);
    const T F  = x1*y1; // multiplication 1
    const T G  = x2*y2; // multiplication 2
    const T t0 = x1+x2;
    const T t1 = y1+y2;
    // this could have overflown if it weren't masked
    T H  = (t0&((~T(0))>>N))*(t1&((~T(0))>>N));  // multiplication 3
    // this is the overflow
    T o  = ((t0>>N)&&(t1>>N));
    H+=   ((t0>>N)?(t1<<N):0);
    o+= H<((t0>>N)?(t1<<N):0);
    H+=   ((t1>>N)?(t0<<N):0);
    o+= H<((t1>>N)?(t0<<N):0);
    T K  = H;
    o-= (K-F)>K;
    K =  K-F;
    o-= (K-G)>K;
    K =  K-G;
    T hi = F;
    T lo = G;
    lo += K<<N;
    hi += lo<(K<<N);
    hi += (K>>N)+(o<<N);
    return tuple<T,T>{hi,lo};
  }

  template<>
  tuple<uint8_t,uint8_t> constexpr long_mul(const uint8_t a,const uint8_t b){
    const int_fast16_t r = int_fast16_t(a)*int_fast16_t(b);
    return {uint8_t(r>>8),uint8_t(r)};
  }
 
  template<>
  tuple<uint16_t,uint16_t> constexpr long_mul(
      const uint16_t a,
      const uint16_t b){
    const uint_fast32_t r = uint_fast32_t(a)*uint_fast32_t(b);
    return {uint16_t(r>>16),uint16_t(r)};
  }
  
  template<>
  tuple<uint32_t,uint32_t> constexpr long_mul(
      const uint32_t a,
      const uint32_t b){
    const uint_fast64_t r = uint_fast64_t(a)*uint_fast64_t(b);
    return {uint32_t(r>>32),uint32_t(r)};
  }
  
  template<typename T>
  typename std::enable_if<std::is_unsigned<T>::value,tuple<T,T>>::type
  constexpr long_mul(
      const tuple<T,T> a,
      const tuple<T,T> b) {
      auto ab = long_mul(get<1>(a),get<1>(b));
      get<0>(ab) += get<0>(a)*get<1>(b)+get<1>(a)*get<0>(b);
      return ab;
  }
  
  template <typename T>
  constexpr size_t popcount(T n){
    size_t c=0;
    while(n) (n=n&(n-1),++c);
    return c;
  }

  constexpr size_t popcount(const uint32_t n){
    return __builtin_popcountl(n);
  }

  constexpr size_t popcount(const uint64_t n){
    return __builtin_popcountll(n);
  }
  
  template <typename T>
  constexpr size_t hamming_distance(const T a, const T b){
    return popcount(a^b);
  }

  // ( a < m ) && ( b < m)
  // ( a - b ) % m
  template <typename T>  
  typename std::enable_if<std::is_unsigned<T>::value,T>::type
  constexpr diff_mod(const T& a,const T& b,const T& m){
    if (a<b) return m-b+a;
    return a-b;
  }
  
  // ( a < m ) && ( b < m )
  template <typename T>  
  typename std::enable_if<std::is_unsigned<T>::value,T>::type
  constexpr add_mod(const T& a,const T& b,const T& m){
    T s = a + b;
    if (s<a) s-=m;
    return ((a%m)+(b%m))%m;
  }
  
  template <typename T,typename S>
  typename std::enable_if<std::is_unsigned<T>::value,T>::type
  constexpr shl(const T n, const S i){
    if (i>=digits<T>()) return 0;
    if (i>=0) return n<<i;
    else      return n>>i;
  }
  
  template <typename T,typename S>
  typename std::enable_if<std::is_unsigned<T>::value,T>::type
  constexpr shr(const T n, const S i){
    if (i>=digits<T>()) return 0;
    if (i>=0) return n>>i;
    else      return n<<i;
  }

  template <typename T,typename S>
  typename std::enable_if<std::is_unsigned<T>::value,T>::type
  constexpr ror(const T n, const S i){
    const T m = (std::numeric_limits<T>::digits-1);
    const T c = i&m;
    return (n>>c)|(n<<((-c)&m));
  }

  template <typename T,typename S>
  typename std::enable_if<std::is_unsigned<T>::value,T>::type
  constexpr rol(const T n, const S i){
    const T m = (std::numeric_limits<T>::digits-1);
    const T c = i&m;
    return (n<<c)|(n>>((-c)&m));
  } 

  template <typename T>
  typename std::enable_if<std::is_unsigned<T>::value,T>::type
  constexpr bitmask(const T& lower, const T& upper){ // exclusive upper limit
    return T(1)<<upper-T(1)<<lower;
  }

  template <typename T>
  typename std::enable_if<std::is_unsigned<T>::value,T>::type
  constexpr bitmask_upper(const T& upper){
    return (upper<std::numeric_limits<T>::digits)?
      T(1)<<upper-T(1):std::numeric_limits<T>::max();
  }

  template <typename T>
  typename std::enable_if<std::is_unsigned<T>::value,T>::type
  constexpr bitmask_lower(const T& lower){
    return T(0)-T(1)<<lower;
  }

  template <typename T>
  typename std::enable_if<std::is_unsigned<T>::value,T>::type
  constexpr alternating_bitmask(const size_t step){
    T mask(0);
    for (size_t i=0;i<digits<T>();i+=2*step){
      mask|=(~T(0)>>(digits<T>()-step))<<i;
    }
    return mask;
  }

  template<class T, class U = typename std::make_unsigned<T>::type>
  constexpr T bswap(T n){
    for (size_t i=digits<unsigned char>();i<digits<T>();i*=2){
      n = ((n&(~(alternating_bitmask<T>(i))))>>i)|
          ((n&( (alternating_bitmask<T>(i))))<<i);
    }
    return n;
  }

#if __GNUC__ > 3 || __clang__
  constexpr uint16_t bswap(uint16_t n){
    return __builtin_bswap16(n);
  }

  constexpr uint32_t bswap(uint32_t n){
    return __builtin_bswap32(n);
  }

  constexpr uint64_t bswap(uint64_t n){
    return __builtin_bswap64(n);
  }
  
  constexpr int16_t bswap(int16_t n){
    return __builtin_bswap16(n);
  }

  constexpr int32_t bswap(int32_t n){
    return __builtin_bswap32(n);
  }

  constexpr int64_t bswap(int64_t n){
    return __builtin_bswap64(n);
  }
#endif

  template<class T,size_t step=1,class U = typename std::make_unsigned<T>::type>
  constexpr U reflect(T n);
  
  template<>
  uint8_t inline reflect(const uint8_t n){
    return (n * 0x0202020202ull & 0x010884422010ull) % 1023;
    // this would be faster on a machine without 64 bit multiplication:
    switch (n){
      case 0b00000000: return 0b00000000;
      case 0b00000001: return 0b10000000;
      case 0b00000010: return 0b01000000;
      case 0b00000011: return 0b11000000;
      case 0b00000100: return 0b00100000;
      case 0b00000101: return 0b10100000;
      case 0b00000110: return 0b01100000;
      case 0b00000111: return 0b11100000;
      case 0b00001000: return 0b00010000;
      case 0b00001001: return 0b10010000;
      case 0b00001010: return 0b01010000;
      case 0b00001011: return 0b11010000;
      case 0b00001100: return 0b00110000;
      case 0b00001101: return 0b10110000;
      case 0b00001110: return 0b01110000;
      case 0b00001111: return 0b11110000;

      case 0b00010000: return 0b00001000;
      case 0b00010001: return 0b10001000;
      case 0b00010010: return 0b01001000;
      case 0b00010011: return 0b11001000;
      case 0b00010100: return 0b00101000;
      case 0b00010101: return 0b10101000;
      case 0b00010110: return 0b01101000;
      case 0b00010111: return 0b11101000;
      case 0b00011000: return 0b00011000;
      case 0b00011001: return 0b10011000;
      case 0b00011010: return 0b01011000;
      case 0b00011011: return 0b11011000;
      case 0b00011100: return 0b00111000;
      case 0b00011101: return 0b10111000;
      case 0b00011110: return 0b01111000;
      case 0b00011111: return 0b11111000;

      case 0b00100000: return 0b00000100;
      case 0b00100001: return 0b10000100;
      case 0b00100010: return 0b01000100;
      case 0b00100011: return 0b11000100;
      case 0b00100100: return 0b00100100;
      case 0b00100101: return 0b10100100;
      case 0b00100110: return 0b01100100;
      case 0b00100111: return 0b11100100;
      case 0b00101000: return 0b00010100;
      case 0b00101001: return 0b10010100;
      case 0b00101010: return 0b01010100;
      case 0b00101011: return 0b11010100;
      case 0b00101100: return 0b00110100;
      case 0b00101101: return 0b10110100;
      case 0b00101110: return 0b01110100;
      case 0b00101111: return 0b11110100;

      case 0b00110000: return 0b00001100;
      case 0b00110001: return 0b10001100;
      case 0b00110010: return 0b01001100;
      case 0b00110011: return 0b11001100;
      case 0b00110100: return 0b00101100;
      case 0b00110101: return 0b10101100;
      case 0b00110110: return 0b01101100;
      case 0b00110111: return 0b11101100;
      case 0b00111000: return 0b00011100;
      case 0b00111001: return 0b10011100;
      case 0b00111010: return 0b01011100;
      case 0b00111011: return 0b11011100;
      case 0b00111100: return 0b00111100;
      case 0b00111101: return 0b10111100;
      case 0b00111110: return 0b01111100;
      case 0b00111111: return 0b11111100;

      case 0b01000000: return 0b00000010;
      case 0b01000001: return 0b10000010;
      case 0b01000010: return 0b01000010;
      case 0b01000011: return 0b11000010;
      case 0b01000100: return 0b00100010;
      case 0b01000101: return 0b10100010;
      case 0b01000110: return 0b01100010;
      case 0b01000111: return 0b11100010;
      case 0b01001000: return 0b00010010;
      case 0b01001001: return 0b10010010;
      case 0b01001010: return 0b01010010;
      case 0b01001011: return 0b11010010;
      case 0b01001100: return 0b00110010;
      case 0b01001101: return 0b10110010;
      case 0b01001110: return 0b01110010;
      case 0b01001111: return 0b11110010;

      case 0b01010000: return 0b00001010;
      case 0b01010001: return 0b10001010;
      case 0b01010010: return 0b01001010;
      case 0b01010011: return 0b11001010;
      case 0b01010100: return 0b00101010;
      case 0b01010101: return 0b10101010;
      case 0b01010110: return 0b01101010;
      case 0b01010111: return 0b11101010;
      case 0b01011000: return 0b00011010;
      case 0b01011001: return 0b10011010;
      case 0b01011010: return 0b01011010;
      case 0b01011011: return 0b11011010;
      case 0b01011100: return 0b00111010;
      case 0b01011101: return 0b10111010;
      case 0b01011110: return 0b01111010;
      case 0b01011111: return 0b11111010;

      case 0b01100000: return 0b00000110;
      case 0b01100001: return 0b10000110;
      case 0b01100010: return 0b01000110;
      case 0b01100011: return 0b11000110;
      case 0b01100100: return 0b00100110;
      case 0b01100101: return 0b10100110;
      case 0b01100110: return 0b01100110;
      case 0b01100111: return 0b11100110;
      case 0b01101000: return 0b00010110;
      case 0b01101001: return 0b10010110;
      case 0b01101010: return 0b01010110;
      case 0b01101011: return 0b11010110;
      case 0b01101100: return 0b00110110;
      case 0b01101101: return 0b10110110;
      case 0b01101110: return 0b01110110;
      case 0b01101111: return 0b11110110;

      case 0b01110000: return 0b00001110;
      case 0b01110001: return 0b10001110;
      case 0b01110010: return 0b01001110;
      case 0b01110011: return 0b11001110;
      case 0b01110100: return 0b00101110;
      case 0b01110101: return 0b10101110;
      case 0b01110110: return 0b01101110;
      case 0b01110111: return 0b11101110;
      case 0b01111000: return 0b00011110;
      case 0b01111001: return 0b10011110;
      case 0b01111010: return 0b01011110;
      case 0b01111011: return 0b11011110;
      case 0b01111100: return 0b00111110;
      case 0b01111101: return 0b10111110;
      case 0b01111110: return 0b01111110;
      case 0b01111111: return 0b11111110;

      case 0b10000000: return 0b00000001;
      case 0b10000001: return 0b10000001;
      case 0b10000010: return 0b01000001;
      case 0b10000011: return 0b11000001;
      case 0b10000100: return 0b00100001;
      case 0b10000101: return 0b10100001;
      case 0b10000110: return 0b01100001;
      case 0b10000111: return 0b11100001;
      case 0b10001000: return 0b00010001;
      case 0b10001001: return 0b10010001;
      case 0b10001010: return 0b01010001;
      case 0b10001011: return 0b11010001;
      case 0b10001100: return 0b00110001;
      case 0b10001101: return 0b10110001;
      case 0b10001110: return 0b01110001;
      case 0b10001111: return 0b11110001;

      case 0b10010000: return 0b00001001;
      case 0b10010001: return 0b10001001;
      case 0b10010010: return 0b01001001;
      case 0b10010011: return 0b11001001;
      case 0b10010100: return 0b00101001;
      case 0b10010101: return 0b10101001;
      case 0b10010110: return 0b01101001;
      case 0b10010111: return 0b11101001;
      case 0b10011000: return 0b00011001;
      case 0b10011001: return 0b10011001;
      case 0b10011010: return 0b01011001;
      case 0b10011011: return 0b11011001;
      case 0b10011100: return 0b00111001;
      case 0b10011101: return 0b10111001;
      case 0b10011110: return 0b01111001;
      case 0b10011111: return 0b11111001;

      case 0b10100000: return 0b00000101;
      case 0b10100001: return 0b10000101;
      case 0b10100010: return 0b01000101;
      case 0b10100011: return 0b11000101;
      case 0b10100100: return 0b00100101;
      case 0b10100101: return 0b10100101;
      case 0b10100110: return 0b01100101;
      case 0b10100111: return 0b11100101;
      case 0b10101000: return 0b00010101;
      case 0b10101001: return 0b10010101;
      case 0b10101010: return 0b01010101;
      case 0b10101011: return 0b11010101;
      case 0b10101100: return 0b00110101;
      case 0b10101101: return 0b10110101;
      case 0b10101110: return 0b01110101;
      case 0b10101111: return 0b11110101;

      case 0b10110000: return 0b00001101;
      case 0b10110001: return 0b10001101;
      case 0b10110010: return 0b01001101;
      case 0b10110011: return 0b11001101;
      case 0b10110100: return 0b00101101;
      case 0b10110101: return 0b10101101;
      case 0b10110110: return 0b01101101;
      case 0b10110111: return 0b11101101;
      case 0b10111000: return 0b00011101;
      case 0b10111001: return 0b10011101;
      case 0b10111010: return 0b01011101;
      case 0b10111011: return 0b11011101;
      case 0b10111100: return 0b00111101;
      case 0b10111101: return 0b10111101;
      case 0b10111110: return 0b01111101;
      case 0b10111111: return 0b11111101;

      case 0b11000000: return 0b00000011;
      case 0b11000001: return 0b10000011;
      case 0b11000010: return 0b01000011;
      case 0b11000011: return 0b11000011;
      case 0b11000100: return 0b00100011;
      case 0b11000101: return 0b10100011;
      case 0b11000110: return 0b01100011;
      case 0b11000111: return 0b11100011;
      case 0b11001000: return 0b00010011;
      case 0b11001001: return 0b10010011;
      case 0b11001010: return 0b01010011;
      case 0b11001011: return 0b11010011;
      case 0b11001100: return 0b00110011;
      case 0b11001101: return 0b10110011;
      case 0b11001110: return 0b01110011;
      case 0b11001111: return 0b11110011;

      case 0b11010000: return 0b00001011;
      case 0b11010001: return 0b10001011;
      case 0b11010010: return 0b01001011;
      case 0b11010011: return 0b11001011;
      case 0b11010100: return 0b00101011;
      case 0b11010101: return 0b10101011;
      case 0b11010110: return 0b01101011;
      case 0b11010111: return 0b11101011;
      case 0b11011000: return 0b00011011;
      case 0b11011001: return 0b10011011;
      case 0b11011010: return 0b01011011;
      case 0b11011011: return 0b11011011;
      case 0b11011100: return 0b00111011;
      case 0b11011101: return 0b10111011;
      case 0b11011110: return 0b01111011;
      case 0b11011111: return 0b11111011;

      case 0b11100000: return 0b00000111;
      case 0b11100001: return 0b10000111;
      case 0b11100010: return 0b01000111;
      case 0b11100011: return 0b11000111;
      case 0b11100100: return 0b00100111;
      case 0b11100101: return 0b10100111;
      case 0b11100110: return 0b01100111;
      case 0b11100111: return 0b11100111;
      case 0b11101000: return 0b00010111;
      case 0b11101001: return 0b10010111;
      case 0b11101010: return 0b01010111;
      case 0b11101011: return 0b11010111;
      case 0b11101100: return 0b00110111;
      case 0b11101101: return 0b10110111;
      case 0b11101110: return 0b01110111;
      case 0b11101111: return 0b11110111;

      case 0b11110000: return 0b00001111;
      case 0b11110001: return 0b10001111;
      case 0b11110010: return 0b01001111;
      case 0b11110011: return 0b11001111;
      case 0b11110100: return 0b00101111;
      case 0b11110101: return 0b10101111;
      case 0b11110110: return 0b01101111;
      case 0b11110111: return 0b11101111;
      case 0b11111000: return 0b00011111;
      case 0b11111001: return 0b10011111;
      case 0b11111010: return 0b01011111;
      case 0b11111011: return 0b11011111;
      case 0b11111100: return 0b00111111;
      case 0b11111101: return 0b10111111;
      case 0b11111110: return 0b01111111;
      case 0b11111111: return 0b11111111;
      
      default: return 0;
    }
  }

  template<class T,size_t step,class U>
  constexpr U reflect(T n) {
    // I thought this was faster but it is not...
    //U m(0);
    //for (size_t i=0;i!=digits<U>();i+=digits<unsigned char>())
    //  m|=U(reflect(uint8_t(n>>i)))<<i;
    //return bswap(m);*/
    for (size_t i=step;i<digits<unsigned char>();i*=2){ 
      n = ((n&(~(alternating_bitmask<T>(i))))>>i)|
        ((n&( (alternating_bitmask<T>(i))))<<i);
    }
    return bswap(n);
  }

  
  /* there are m! possible bijections for integers modulo m
   * for a 64 bit integer this is a lot...
   * 2**(2**(69.9671))
   * = 2**(64*1.8e19)
   *
   * addition with another integer
   *
   * addition with carry
   *
   * xor with another integer
   *
   * negation
   *
   * multiplication with an odd integer
   *
   * rotations
   *
   * multiplication modulo prime and leaving numbers greater than modulus
   * untouched
   *
   * clmul: carryless multiplication and
   * reduction with an irreducible polynomial on GF(2^k)
   *
   */

  template<typename T>
  typename std::enable_if<std::is_unsigned<T>::value,T>::type
  constexpr clmul(const T& a,const T& b){
    T result(0);
    for (size_t i=0;i!=digits(a);++i){
      if (a&(T(1)<<i)) result^=b<<i;
    }
    return result;
  }

  /* circmul with an odious integer is a bijection on the integers modulo m
   * odious := odd number of 1s in binary
   * popcount(b)%2==1
   * 2**m-1
   * can be represented as sum of rotations of a (?)
   */
  template<typename T>
  typename std::enable_if<std::is_unsigned<T>::value,T>::type
  constexpr circmul(const T& a, const T& b){
    /*const size_t s = digits(a)/2;
    const T m  = (~T(0))>>s;
    const T t0 = (a>>s)*(b>>s);
    const T t1 =  (a&m)*(b>>s);
    const T t2 = (a>>s)*(b&m);
    const T t3 =  (a&m)*(b&m);
    T        r = circadd(t0,t3);
             r = circadd(r ,ror(t1,s));
             r = circadd(r ,ror(t2,s));
    return r;*/
    auto r = long_mul(a,b);
    return get<0>(r)+get<1>(r);
    //return circadd(get<0>(r),get<1>(r));
  }
  
  uint8_t constexpr clmul_mod(uint8_t a, uint8_t b){
    uint8_t p=0;
    for (uint8_t i=0;i!=8;++i){
      p^=-(b&1)&a;
      uint8_t m = -((a>>7)&1);
      a = (a<<1)^(27&m);
      b>>=1;
    }
    return p;
  }
  
  // 01234567
  // 11000101
  // 0157----

  template<typename T>
  typename std::enable_if<std::is_unsigned<T>::value,T>::type
  const inline pext(const T& n, const T& mask){
    T result(0);
    for (size_t i=0,j=0;i!=digits(n);++i){
      const T b = T(1)<<i;
      if (mask&b) result|=b&(n>>(i-j));
      else ++j;
    }
    return result;
  }
  
  // 01234567
  // 11000101
  // 01---5-7

  template<typename T>
  typename std::enable_if<std::is_unsigned<T>::value,T>::type
  const inline pdep(const T& n, const T& mask){
    T result(0);
    for (size_t i=0;i!=digits(n);++i){
      const T b = T(1)<<i;
      if (mask&b) result|=b&n;
    }
    return result;
  }


#if defined (__has_include) && (__has_include(<x86intrin.h>))
  
#if defined (__BMI2__)
  uint16_t const inline pdep(const uint16_t& n,const uint16_t& mask){
    return _pdep_u32(n,mask);
  }

  uint32_t const inline pdep(const uint32_t& n,const uint32_t& mask){
    return _pdep_u32(n,mask);
  }
  
  uint64_t const inline pdep(const uint64_t& n,const uint64_t& mask){
    return _pdep_u64(n,mask);
  }
  
  uint16_t const inline pext(const uint16_t& n,const uint16_t& mask){
    return _pext_u32(n,mask);
  }

  uint32_t const inline pext(const uint32_t& n,const uint32_t& mask){
    return _pext_u32(n,mask);
  }
  
  uint64_t const inline pext(const uint64_t& n,const uint64_t& mask){
    return _pext_u64(n,mask);
  }
#endif // __BMI2__

  uint16_t const inline clmul(const uint16_t& i,const uint16_t& j){
    __m128i I{}; I[0]^=i;
    __m128i J{}; J[0]^=j;
    __m128i X = _mm_clmulepi64_si128(I,J,0);
    return X[0];
  }
  
  uint32_t const inline clmul(const uint32_t& i,const uint32_t& j){
    __m128i I{}; I[0]^=i;
    __m128i J{}; J[0]^=j;
    __m128i X = _mm_clmulepi64_si128(I,J,0);
    return X[0];
  }
  
  uint64_t const inline clmul(const uint64_t& i,const uint64_t& j){
    __m128i I{};I[0]^=i;
    __m128i J{};J[0]^=j;
    __m128i X = _mm_clmulepi64_si128(I,J,0);
    return X[0];
  }

  uint16_t const inline clmul_mod(const uint16_t& i,const uint16_t& j){
    __m128i I{}; I[0]^=i;
    __m128i J{}; J[0]^=j;
    __m128i M{}; M[0]^=29189u;
    __m128i X = _mm_clmulepi64_si128(I,J,0);
    __m128i X0{};X0[0]^=X[0]&(~uint16_t(0));
    __m128i A = _mm_clmulepi64_si128(X0,M,0);
    __m128i A0{};A0[0]^=A[0]&(~uint16_t(0));
    __m128i B = _mm_clmulepi64_si128(A0,M,0);
    return A[0]^(A[0]>>16)^(B[0]>>16)^X[0]^(X[0]>>16);
  }
  
  uint32_t const inline clmul_mod(const uint32_t& i,const uint32_t& j){
    __m128i I{}; I[0]^=i;
    __m128i J{}; J[0]^=j;
    __m128i M{}; M[0]^=1073877003u;
    __m128i X = _mm_clmulepi64_si128(I,J,0);
    __m128i X0{};X0[0]^=X[0]&(~uint32_t(0));
    __m128i A = _mm_clmulepi64_si128(X0,M,0);
    __m128i A0{};A0[0]^=A[0]&(~uint32_t(0));
    __m128i B = _mm_clmulepi64_si128(A0,M,0);
    return A[0]^(A[0]>>32)^(B[0]>>32)^X[0]^(X[0]>>32);
  }
  
  uint64_t const inline clmul_mod(const uint64_t& i,const uint64_t& j){
    __m128i I{};I[0]^=i;
    __m128i J{};J[0]^=j;
    __m128i M{};M[0]^=0xb000000000000000ull;
    __m128i X = _mm_clmulepi64_si128(I,J,0);
    __m128i A = _mm_clmulepi64_si128(X,M,0);
    __m128i B = _mm_clmulepi64_si128(A,M,0);
    return A[0]^A[1]^B[1]^X[0]^X[1];
  }
  
  template<typename T>
  typename std::enable_if<std::is_unsigned<T>::value,T>::type
  const inline clmul_mod_pow(
            T  i,
      const T& j,
       size_t  k){
    uint64_t r = 0;
    while(k){
      if (k&size_t(1)) r^=i;
      i=clmul_mod(i,i);
      k>>=1;
    }
    return r;
  }

  template<typename T>
  typename std::enable_if<std::is_unsigned<T>::value,T>::type
  const inline clmul_mod_inverse(
      const T& A
      ) {
    return clmul_mod_pow(A,T(0)-T(2));
  }
  
  uint8_t const clmul_circ(const uint8_t& a,const uint8_t& b){
    int64_t _a = a;
    int64_t _b = b;
    __m128i ma{_a,0ull};
    __m128i mb{_b,0ull};
    auto t = _mm_clmulepi64_si128(ma,mb,0);
    return t[0]^(t[0]>>8);
  }
  
  uint16_t const clmul_circ(const uint16_t& a,const uint16_t& b){
    int64_t _a(0);
    int64_t _b(0);
    _a^=a;
    _b^=b;
    __m128i ma{_a,0ull};
    __m128i mb{_b,0ull};
    auto t = _mm_clmulepi64_si128(ma,mb,0);
    return t[0]^(t[0]>>16);
  }
  
  uint32_t const clmul_circ(const uint32_t& a,const uint32_t& b){
    int64_t _a(0);
    int64_t _b(0);
    _a^=a;
    _b^=b;
    __m128i ma{_a,0ull};
    __m128i mb{_b,0ull};
    auto t = _mm_clmulepi64_si128(ma,mb,0);
    return t[0]^(t[0]>>32);
  }
  
  uint64_t const clmul_circ(const uint64_t& a,const uint64_t& b){
    int64_t _a(0);
    int64_t _b(0);
    _a^=a;
    _b^=b;
    __m128i ma{_a,0ull};
    __m128i mb{_b,0ull};
    auto t = _mm_clmulepi64_si128(ma,mb,0);
    return t[0]^t[1];
  }

#else

  uint16_t constexpr clmul_mod(uint16_t a, uint16_t b){
    uint16_t p=0;
    for (uint16_t i=0;i!=16;++i){
      p^=-(b&1)&a;
      uint16_t m = -((a>>15)&1);
      a = (a<<1)^(29189u&m);
      b>>=1;
    }
    return p;
  }
  
  uint32_t constexpr clmul_mod(uint32_t a, uint32_t b){
    uint32_t p=0;
    for (uint32_t i=0;i!=32;++i){
      p^=-(b&1)&a;
      uint32_t m = -((a>>31)&1);
      a = (a<<1)^(1073877003u&m);
      b>>=1;
    }
    return p;
  }

  uint64_t constexpr clmul_mod(uint64_t a, uint64_t b){
    uint64_t p=0;
    for (uint64_t i=0;i!=64;++i){
      p^=-(b&1)&a;
      uint64_t m = -((a>>63)&1);
      a = (a<<1)^(0b11011ul&m);
      b>>=1;
    }
    return p;
  }
  
  template <typename T>
  typename std::enable_if<std::is_unsigned<T>::value,T>::type
  constexpr clmul_circ(const T& i,const T& j){
    T r=0;
    for (size_t n=0;n!=digits<T>();++n) if(j&(T(1)<<n)) r^=wmath::rol(i,n);
    return r;
  }
#endif

  template <typename T>
  typename std::enable_if<std::is_unsigned<T>::value,T>::type
  constexpr clz(const T x,const T lower=0,const T upper=digits<T>()){
    return 
      (upper-lower==T(1))
        ?digits<T>()-upper
        :(x&(T(0)-T(1)<<((upper+lower)/2))?
           clz(x,(upper+lower)/2,upper):
           clz(x,lower,(upper+lower)/2));
  }
 
  template <typename T>
  typename std::enable_if<std::is_unsigned<T>::value,T>::type
  constexpr ctz(const T x,const T lower=0,const T upper=digits<T>()){
    return
      (upper-lower==T(1))
        ?lower
        :(x&(T(0)-T(1)<<((upper+lower)/2))?
          ctz(x,(upper+lower)/2,upper):
          ctz(x,lower,(upper+lower)/2));
  }


  template <typename T>
  typename std::enable_if<std::is_unsigned<T>::value,T>::type
  constexpr log2(const T x,const T lower=0,const T upper=digits<T>()){
    return (upper-lower==T(1))?lower:(x&(T(0)-T(1)<<((upper+lower)/2))?
           log2(x,(upper+lower)/2,upper):
           log2(x,lower,(upper+lower)/2));
  }

#if __GNUC__ > 3 || __clang__
  uint32_t constexpr clz(const uint32_t x){
    return x==0?32:__builtin_clz(x);
  }

  uint32_t constexpr ctz(const uint32_t x){
    return x==0?0:__builtin_ctz(x);
  }
 
  uint64_t constexpr clz(const uint64_t x){
    return x==0?64:__builtin_clzll(x);
  }

  uint64_t constexpr ctz(const uint64_t x){
    return x==0?0:__builtin_ctzll(x);
  }

  uint32_t constexpr log2(const uint32_t x){
    return 31-__builtin_clz(x);
  }
  
  uint64_t constexpr log2(const uint64_t x){
    return 63-__builtin_clzll(x);
  }
#ifdef __SIZEOF_INT128__
  __uint128_t constexpr log2(const __uint128_t x){
    if (x>>64) return log2(uint64_t(x>>64));
    return log2(uint64_t(x&(~uint64_t(0))));
  }
#endif
#endif

  /*
  // 2**digits(n)/n
  template <typename T>
  typename std::enable_if<std::is_unsigned<T>::value,T>::type
  constexpr inverse(const T& n){
    if (n==0) return n;
    return (~T(0))/n+(popcount(n)==1);
  }
  */
  
  // modular multiplicative inverse for power two
  // modular inverses only exist for uneven d
  template<typename T>
  T constexpr modular_inverse_power2(const T& d,const T& m = T(0)) {
    T x(1); // minimally more elegant than x = d
    //   1 + 2 + 4 + 8 + 16 + 32 + 64 ... // new correct digits starting with 1
    //   1   3   7  15   31   63  127 ... // total correct digits
    //j= 0   1   2   3    4    5    6 ... -> log2(n) is the correct bound
    for (size_t j=0;j!=wmath::log2(digits(m)-clz(m-1));++j) x*=T(2u)-x*d;
    return x&(m-1);
  }
  
  // lookup tabe is just not worth it
  template<typename T>
  T constexpr fixpoint_integer_inverse(const T& d) {
    T x = T(1)<<clz(d-1);
    for (size_t j=0;j!=wmath::log2(digits<T>());++j) 
      x+=get<0>(long_mul(x,T(1)-(x*d)));
    return x;
  }
 
  template<typename T,size_t n>
  array<T,n> constexpr fixpoint_integer_inverse(const array<T,n>& d) {
    array<T,n> d1 = sub(d,array<T,1>{1}); // d1 = d-1
    array<T,n> x{};
    size_t p = 0;
    for (p=n-1;p!=0;--p) if (d1[p]) break;
    //cout << p << " " << n-p-1 << endl;
    x[n-p-1] = T(1)<<clz(d1[p]);
    for (size_t j=0;j!=wmath::log2(n*digits<T>());++j) {
      //cout << "x= ";
      //for (auto it=x.begin();it!=x.end();++it) {
      //  cout << *it << " ";
      //}
      //cout << endl;
      //cout << "x*d= ";
      const auto xd = long_mul(x,d);
      //for (auto it=xd.begin();it!=xd.end();++it) {
      //  cout << *it << " ";
      //}
      //cout << endl;
      array<T,n> lo;copy(xd.begin(),xd.begin()+n,lo.begin());
      lo = sub(array<T,1>{1},lo);
      //for (auto it=lo.begin();it!=lo.end();++it) *it=~*it;
      const auto dx = long_mul(x,lo);
      //cout << "dx= ";
      //for (auto it=dx.begin();it!=dx.end();++it) {
      //  cout << *it << " ";
      //}
      //cout << endl;
      add_with_carry(dx.begin()+n,dx.end(),x.begin(),x.end());
    }
    return x;
  }
  
  /*
  template<typename T>
  tuple<T,T> constexpr fixpoint_integer_inverse(const tuple<T,T>& d) {
    tuple<T,T> x;
    if (get<0>(d1)) {
      get<1>(x)=T(1)<<clz(get<0>(d1));
    } else {
      get<0>(x)=T(1)<<clz(
    }
    for (size_t j=0;j!=wmath::log2(digits<T>());++j) 
      x+=get<0>(long_mul(x+1,T(0)-(x*d)));
    return x;
  }
  */
  
  template <typename T>
  typename std::enable_if<std::is_unsigned<T>::value,T>::type
  constexpr addmod(const T& a,const T& b,const T& p){
    return ((a+b)<a?(a+b-p)%p:(a+b)%p);
  }

  // ( a * b ) % m with the help of i , this is Barrett reduction
  // 1 multiplication n ,  n
  // 1 multiplication n , 2n
  // -> 2n log(2n) + 3n log(3n)
  // = n ( 5 log(n) + 2 log(2) + 3 log(3) )
  // ~ 5 n log(n)
  // Montgommery uses two multiplication n , n
  // -> 4 n log(4 n)
  // = n ( 4 log(n) + 4 log(4) )
  // ~ 4 n log(n)
  // -> Barrett is 25% less efficient
  // Montgomery would be 20% more efficient
  // below is one (2n,2n) multiplication, where the lower half is not needed
  // and one (2n,n) multiplication, where the top most part is guaranteed to be 0
  // if m is close to 2**n i[1] is 0
  // and only the middle part of the result is needed -> 2 n log(4 n)
  // -> Barrett equivalent to Montgomery?
  // yes. and it can be even faster...
  // Wiliam Hasenplaugh, Gunnar Gaubatz and Vinodh Gopal
  // Mongomery is dead, long live Barrett
  template <typename T>
  typename std::enable_if<std::is_unsigned<T>::value,T>::type
  constexpr mulmod(const T& a,const T& b,const T& m,const array<T,2>& i){
    auto ab = long_mul(a,b);
    // cout << i[0] << "+2**64*" << i[1] << endl;
    // quotient, maybe one too little
    array<T,4> q = long_mul(array<T,2>{get<1>(ab),get<0>(ab)},i);
    // auto q = get<1>(long_mul(i[0],get<0>(ab))); // if clz(m)==0
    //cout << q[0] << " " << q[1] << " " << q[2] << " " << q[3] << endl;
    // trial to see if too little
    array<T,3> t  = long_mul(array<T,2>{q[2],q[3]},array<T,1>{m});
    // auto t = long_mul(q,m);
    //cout << t[0] << " " << t[1] << " " << t[2] << endl;
    array<T,2> r {get<1>(ab),get<0>(ab)};            // remainder
    //cout << r[0] << " " << r[1] << endl;
    sub_with_carry(t.begin(),t.end(),r.begin(),r.end());
    //cout << r[0] << " " << r[1] << endl; 
    const T c = ((r[1])||(r[0]>=m))?m:T(0);
    r[0]-=c;      // no if there, so it is less efficient, but side channel
    r[1]-=r[0]>c; // attacks. If this were to be used in a cryptographic context.
    return r[0];
  }

  template <typename T>
  typename std::enable_if<std::is_unsigned<T>::value,T>::type
  constexpr mulmod(const T& a,const T& b,const T& m){
    return mulmod(a,b,m,fixpoint_integer_inverse(array<T,2>{m}));
  }

  /* Modular multiplication made into a bijection on the fixed size integers */ 
  template <typename T>
  typename std::enable_if<std::is_unsigned<T>::value,T>::type
  constexpr bijective_mulmod(const T& a, const T& b){
    const T p = largest_prime(T(0));
    return a>p?a:mulmod(a,b,p);
  }

}
#endif // WMATH_BITS_H
