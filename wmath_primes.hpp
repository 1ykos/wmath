#ifndef WMATH_PRIMES_H
#define WMATH_PRIMES_H
#include "wmath_forward.hpp"

namespace wmath{
 
  /* depends on a good potentiation with modulus reduction
  template <typename T>
  typename std::enable_if<std::is_integral<T>::value,bool>::type
  miller_rabin(
      const T& n,
      const T& a
      )
  {
    if (n%2==0) return 0;
    auto d = n-1;
    // if there was a count trailing zeros intrinsic I'd use that instead
    while (d&1u) d>>=1;
    const auto r = pow_mod(a,d,n);
    return (r==n-1)||(r==1);
  }
  */

  /* deterministic miller-rabin, probably not worth it
  Input: n > 1, an odd integer to be tested for primality
  Output: “composite” if n is composite, “prime” otherwise

  write n as 2r·d + 1 with d odd (by factoring out powers of 2 from n − 1)
  WitnessLoop: for all a in the range [2, min(n−2, ⌊2(ln n)2⌋)]:
    x ← ad mod n
    if x = 1 or x = n − 1 then
        continue WitnessLoop
    repeat r − 1 times:
        x ← x2 mod n
        if x = n − 1 then
            continue WitnessLoop
    return “composite”
  return “prime”
  */
  /*
  is_prime(const T& n) {
    // miller rabin with  2  3     sufficient below  1373653
    // miller rabin with 31 73     sufficient below  9080191
    // miller rabin with  2  3  5  sufficient below 25326001
    // miller rabin with  2 13 23 1662803 below 1122004669633
    // miller rabin with  2  3  5  7 11 below 2152302898747
    // 2 3 5 7 11 13 17 19 23 29 31 37 below 318'665'857'834'031'151'167'461
  }
  */

  template <typename T>
  typename std::enable_if<std::is_integral<T>::value,bool>::type
  constexpr is_prime(const T& n) {
    if (n<2) return false;
    if (n%2 ==0)  return n==2;
    if (n%3 ==0)  return n==3;
    if (n%5 ==0)  return n==5;
    if (n<49) return true;
    const T m = n%30;
    switch(m) {
      case 1:
      case 7:
      case 11:
      case 13:
      case 17:
      case 19:
      case 23:
      case 29:
        break;
      default:
        return false;
    }
    for (T i=7;i*i<=n;i+=2){
      if (n%i==0) return false;
    }
    return true;
  }
  
  template <typename T>
  typename std::enable_if<std::is_integral<T>::value,vector<T>>::type
  const prime_factorization(const T& n){
    if (n<2) return {};
    if (n%2 ==0){
      T m = n/2;
      while (m%2==0) m/=2;
      auto v = prime_factorization(m);
      v.push_back(2);
      return v;
    }
    if (n%3 ==0){
      T m = n/3;
      while (m%3==0) m/=3;
      auto v = prime_factorization(m);
      v.push_back(3);
      return v;
    }
    if (n%5 ==0){
      T m = n/5;
      while (m%5==0) m/=5;
      auto v = prime_factorization(m);
      v.push_back(5);
      return v;
    }
    for (size_t k=7,i=1;k*k<=n;){
      if (n%k==0){
        T m = n/k;
        while (m%k==0) m/=k;
        auto v = prime_factorization(m);
        v.push_back(k);
        return v;
      }
      switch(i){
        case 0: k+=6; i=1; break; // 1
        case 1: k+=4; i=2; break; // 7
        case 2: k+=2; i=3; break; //11
        case 3: k+=4; i=4; break; //13
        case 4: k+=2; i=5; break; //17
        case 5: k+=4; i=6; break; //19
        case 6: k+=6; i=7; break; //23
        case 7: k+=2; i=0; break; //29
      }
    }
    return {n};
  }
  
  template <typename T>
  typename std::enable_if<std::is_unsigned<T>::value,T>::type
  constexpr largest_prime(const T& i=0);

  template <>
  uint8_t constexpr largest_prime(const uint8_t& i){
    switch (i){
      case 0:
        return 251;
      case 1:
        return 241;
      case 2:
        return 239;
      case 3:
        return 233;
      default:
        return 1;
    }
  }

  template <>
  uint16_t constexpr largest_prime(const uint16_t& i){
    switch (i){
      case 0:
        return uint16_t(0)-uint16_t(15);
      case 1:
        return uint16_t(0)-uint16_t(17);
      case 2:
        return uint16_t(0)-uint16_t(39);
      case 3:
        return uint16_t(0)-uint16_t(57);
      default:
        return 1;
    }
  }

  template <>
  uint32_t constexpr largest_prime(const uint32_t& i){
    switch (i){
      case 0:
        return uint32_t(0)-uint32_t(5);
      case 1:
        return uint32_t(0)-uint32_t(17);
      case 2:
        return uint32_t(0)-uint32_t(65);
      case 3:
        return uint32_t(0)-uint32_t(99);
      default:
        return 1;
    }
  }

  template <>
  uint64_t constexpr largest_prime(const uint64_t& i){
    switch (i){
      case 0:
        return uint64_t(0)-uint64_t(59);
      case 1:
        return uint64_t(0)-uint64_t(83);
      case 2:
        return uint64_t(0)-uint64_t(95);
      case 3:
        return uint64_t(0)-uint64_t(179);
      default:
        return 1;
    }
  }

  template <typename T>
  typename std::enable_if<std::is_unsigned<T>::value,T>::type
  constexpr random_prime(const T& i=0);

  template <>
  uint8_t constexpr random_prime(const uint8_t& i){
    uint8_t primes[] = {127u,
                        131u,
                        251u,
                        223u,
                        191u,
                        193u,
                        47u,
                        97u};
    return primes[i];
  }
  
  template <>
  uint16_t constexpr random_prime(const uint16_t& i){
    uint16_t primes[] = {42307u,
                         52313u,
                         51307u,
                         11317u,
                         60317u,
                         60337u,
                         60037u,
                         30137u};
    return primes[i];
  }

  template <>
  uint32_t constexpr random_prime(const uint32_t& i){
    uint32_t primes[] = {4184867191u,
                         4184864197u,
                         4184411197u,
                         3184410197u,
                         2184200197u,
                          728033399u,
                         1061068399u,
                         3183208117u};
    return primes[i];
  }
  
  template <>
  uint64_t constexpr random_prime(const uint64_t& i){
    uint64_t primes[] = {15112557877901478707ul,
                         18446744073709503907ul,
                          5819238023569667969ul,
                         17457704070697003907ul,
                         14023704271282629773ul,
                         15457704070697023907ul,
                         12023704271182029287ul,
                         10023704271182029357ul,
                          8023704271998834967ul};
    return primes[i];
  }
}
#endif // WMATH_PRIMES_H
