#ifndef WMATH_ENCODE_H
#define WMATH_ENCODE_H
namespace wmath{
  // levensthein coding in the range [0,2^51]
  uint64_t constexpr universal_encode_uint64(uint64_t n){
    uint64_t log2star = 0;
    uint64_t r=0;
    while (n){
      r>>=log2(n);
      r|=n<<(64-log2(n)); //prepend everything but the msb from n to r
      n=log2(n);
      ++log2star;
    }
    r>>= 1+log2star;
    r |= (~uint64_t(0))<<(64-log2star);
    return r;
  }

  uint64_t constexpr universal_encode_uint64(uint64_t n,uint64_t& l){
    if (n==0){
      l=1;
      return 0;
    }
    if (n==1){
      l=2;
      return 1ull<<63;
    }
    uint64_t log2star = 0;
    uint64_t r=0;
    while (n){
      r>>=log2(n);
      l+=log2(n);
      r|=n<<(64-log2(n)); //prepend everything but the msb from n to r
      n=log2(n);
      ++log2star;
    }
    r>>=1+log2star;
    l+=1+log2star;
    r |= (~uint64_t(0))<<(64-log2star);
    return r;
  }
  
  uint64_t constexpr universal_decode_uint64(uint64_t i){
    uint64_t log2star = 0;
    while (i&(uint64_t(1)<<(63-log2star))) ++log2star; // could be done with count leading zeroes
    i<<=1+log2star;
    uint64_t n = log2star>0;
    while (log2star>1){
      const uint64_t _n=n;
      n  = (i>>(64-n))+(uint64_t(1)<<n);
      i<<= _n;
      --log2star;
    }
    return n;
  }

  uint64_t constexpr universal_decode_uint64(uint64_t i,uint64_t& l){
    if (i==0){
      l=1;
      return 0;
    }
    uint64_t log2star = 0;
    while (i&(uint64_t(1)<<(63-log2star))) ++log2star; // could be done with count leading zeroes
    i<<=1+log2star;
    l+=1+log2star;
    uint64_t n = log2star>0;
    while (log2star>1){
      const uint64_t _n=n;
      n  = (i>>(64-n))+(uint64_t(1)<<n);
      i<<= _n;
      l+=_n;
      --log2star;
    }
    return n;
  }

  uint64_t inline length_of_universal_code(uint64_t n){
    uint64_t r=0;
    while (n){
      r+=log2(n);
      n=log2(n);
      ++r;
    }
    ++r;
    return r;
  }

  /*
  //f(a,x)=(x<0)?f(a,x**a)-log(x)/log(a):(x<1)?1-(0.25-(x-0.5)**2)/(exp(1)-1):log(x)/log(a)+f(a,log(x)/log(a));
  double constexpr slog(const double& x){
    return (x<1)?0.5*(x*x-x)+1:log(x)+slog(log(x));
  }

  
  template<typename T>
  typename std::enable_if<std::is_unsigned<T>::value,double>::type
  constexpr universal_number_distribution(const T& n){
    return exp(-slog(n+1)-exp(-1));
  }
  */

  
  //Zigzag decoding and encoding for uint32_t and uint64_t
  //encoding:
  //0 -> 0; -1 -> 1; 1 -> 2; -2 ->  3; 2 -> 4 ...
  //decoding:
  //0 -> 0;  1 ->-1; 2 -> 1;  3 -> -2; 4 -> 2 ...
  template<typename T,class U = typename std::make_unsigned<T>::type>
  U constexpr zigzag_encode(const T &n){
    return (U(n)<<1)^(n<0?~U(0):U(0));
  }

  template<typename T,class S = typename std::make_signed<T>::type>
  S constexpr zigzag_decode(const T &n){
      return (n&1)?-S(n>>1)-1:(n>>1);
  }
  
  template<typename T>
  typename std::enable_if<std::is_unsigned<T>::value,T>::type
  constexpr log2c1p(const T& d){
    //log2c1p(x)=(22.0*x-(22.0-16)*x**2)/16.0
    // (22*d-6*d*d)/16
    const T t0 = get<0>(long_mul(d,T(22)));
    const T t1 = get<1>(long_mul(d,T(22)));
    const T t2 = pow(d>>(digits<T>()/2),2);
    const T t3 = get<0>(long_mul(t2,T(6)));
    const T t4 = get<1>(long_mul(t2,T(6)));
    const T t5 = t1-t4;
    const T t6 = t0-t3-(t5>t1);
    const T t7 = ror(t6,4)+(t5>>4);
    //const T t7 = (t6<<(digits<T>()-4))+(t5>>4);
    return t7;
  }
  
  template<typename T>
  typename std::enable_if<std::is_unsigned<T>::value,tuple<T,T>>::type
  constexpr log2c(const T& n,const T& d = 0){
    //log2c(x)=(x>=1)?flog2(x)+log2c1p((x-2.0**flog2(x))/(2.0**flog2(x))):-log2c(1.0/x);
    if (n){ // x>=1
      const T t = log2(n);
      return {t,log2c1p((d>>t)+ror(n,t))};
    }
    return {~T(0),~T(0)};
  }
  
  template<typename T>
  typename std::enable_if<std::is_unsigned<T>::value,T>::type
  constexpr universal_cumulative_distribution(const T& n,const T& d = 0){
    //ESLOG2(x)=(x>1)?ESLOG2(log(x)/log(2))/2+PESLOG2(1):9*x/16.0-x**2/16.0;
    //ESLOG2(x)=(x>1)?ESLOG2(log(x)/log(2))/2+PESLOG2(1):x/2;
    //ESLOG2(x)=(x>1)?ESLOG2(log2c(x))/2.0+0.5:(9.0*x-x**2)/16.0;
    if (n) return (T(1)<<(digits(n)-1))
          +(universal_cumulative_distribution(
                std::get<0>(log2c(n,d)),
                std::get<1>(log2c(n,d)))>>1);
    return d>>1;
  }

  const double universal_distribution(const double& x){
    // eslog2(x)=(x<=1)?9.0/16-0.125*x:eslog2(log(x)/log(2))/(x*2*log(2))
    if (x>1) return universal_distribution(std::log2(x))/(x*2*log(2));
    return 9.0/16-0.125*x;
  }
  
  // starting with 0.lo , 1.hi where lo is T(0) and hi is T(0)
  // ranges are meant [lo,hi) i.e. including lo, excluding hi
  // TODO test
  // TODO write decoder... oO
  template<typename T>
  typename std::enable_if<std::is_unsigned<T>::value,tuple<T,T,T,T>>::type
  constexpr arithmetic_encode(const T& lo,const T& hi,
                              const T& c0,const T& c1
      ){
    // lo <- lo + (hi-lo)*c0;
    // hi <- hi + (hi-lo)*c1;
    const T diff = hi-lo; // 1.diff
    auto [lo0,lo1] = long_mul(diff,c0);
    auto [hi0,hi1] = long_mul(diff,c1);
    T t0,t1;
    t0 = lo0+c0;
    t1 = t0+lo;
    lo0= t1;
    const T lb0 = (t0<c0)||(t1<t0);
    t0 = hi0+c1;
    t1 = t0+lo;
    hi0 = t1;
    const T hb0 = (t0<c1)||(t1<t0);
    if (lb0!=hb0) return {0,0,lo0,hi0};
    const auto n = clz(~(lo0^hi0)); // n+1 bits are shifted out
    const auto s = digits<T>()-n;   // s-1 bits stay 
    return {n,(lb0<<(s+1))|shr(lo0,s),
      (lo0<<(n+1))|(lo1>>(s-1)),
      (hi0<<(n+1))|(hi1>>(s-1))
    };
  }

  /*
  // returns the expanded fraction obtained by multiplying the two fractions
  // given by the digits in lo0 hi0 and lo1 hi1
  template<typename T>
  typename std::enable_if<std::is_unsigned<T>::value,tuple<T,T,T,T>>::type
  constexpr arithmetic_encode_expand(const T& lo0,const T& hi0,
                                     const T& lo1,const T& hi1){
    const T N = digits<T>();
    const T n = N/2;
    const T diff = hi0-lo0;
    const T rlo0 = lo0+get<0>(long_mul(diff,lo1));
    const T rlo1 =     get<1>(long_mul(diff,lo1));
    const T t0   = hi1+get<1>(long_mul(diff,hi1));
    const T rhi1 = t0+diff;
    const T rhi0 = lo0+get<0>(long_mul(diff,hi1))+(t0<hi1)+(rhi1<t0);
    return {rlo0,rlo1,rhi0,rhi1};
  }
  
  // returns the number of equal bits to shift out, the bits to shift out in
  // two values of type T and the remaining fraction
  template<typename T>
  typename std::enable_if<std::is_unsigned<T>::value,tuple<T,T,T,T,T>>::type
  constexpr arithmetic_encode_compress(const T& lo0,const T& lo1,
                                       const T& hi0,const T& hi1){
    const T N = digits<T>();
    const T shift = clz(~(lo0^hi0));
    if (shift==N){
      const T shift = clz(~(lo1^hi1));
      return {shift+N,lo0,lo1>>(N-shift),hi0,h1>>(N-shift)};
    }
    return {shift,T(0),
      lo0>>(N-shift),(lo0<<shift)|(lo1>>(N-shift)),
      hi0>>(N-shift),(hi0<<shift)|(hi1>>(N-shift))};
  }

  template<typename T,typename B>
  typename std::enable_if<std::is_unsigned<T>::value,void>::type
  const inline arithmetic_encode(B& bitstream,T& lo0,T& hi0,
      const T& lo1,const T& hi1){
    const T N = digits<T>();
    const auto t0 = arithmetic_encode_expand
      (lo0,hi0,lo1,hi1);
    const auto t1 = arithmetic_encode_compress
      (get<0>(t0),get<1>(t0),get<2>(t0),get<3>(t0));
    if (get<0>(t1)>N){
      bitstream.write(&get<1>(t1),get<0>(t1)-N);
      bitstream.write(&get<2>(t1),N);
      return;
    }
    bitstream.write(&get<2>(t1),get<0>(t1));
  }

  // state is the remaining probability that could not be shifted to output
  // [low/max,high/max] is the range in the cumulative distribution where the
  // probability of the symbol to encode lies
  // returns {low,high}
  template<typename T>
  typename std::enable_if<std::is_unsigned<T>::value,tuple<T,T>>::type
  constexpr arithmetic_encode(const T& lo0,const T& hi0,
                              const T& lo1,const T& hi1){
    const T diff = hi0-lo0;
    const T 
    (hi0-lo0)*lo0
    const T shift = clz(~lo);
    ///////////////////////////////////////////////////////////////
    // find the number of equal bits in range, then shift them out
    const T shift = clz(~(lo^hi));
    // output lo>>(digits(lo)-shift)
    // update: lo<<shift; (hi<<shift)+(~T(0)>>(digits(hi)-shift)
    return {low<<shift};
  }
  */
}
#endif // WMATH_ENCODE_H
