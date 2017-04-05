#ifndef WMATH_H
#define WMATH_H

#include <algorithm>
#include <array>
#include <cstdint>
#include <functional>
#include <iostream>
#include <iterator>
#include <limits>
#include <numeric>
#include <vector>
#include <random>
#include <iomanip>
//#include <dlib/matrix.h>
//#include <Eigen/Dense>
#include <bitset>

#define VERBOSE false

/*
namespace std{
  template<typename T, typename = void>
  struct is_iterator{
    static constexpr bool value = false;
  };
  template<typename T>
  struct is_iterator<T,
    typename std::enable_if<
      !std::is_same<
        typename std::iterator_traits<T>::value_type,
        void
      >::value
    >::type
  >{
    static constexpr bool value = true;
  };
}
*/

namespace wmath{

  using std::abs;
  using std::accumulate;
  using std::advance;
  using std::array;
  using std::cerr;
  using std::cout;
  using std::distance;
  using std::enable_if;
  using std::endl;
  using std::is_floating_point;
  using std::is_integral;
  using std::is_unsigned;
  using std::iterator;
  using std::iterator_traits;
  using std::min_element;
  using std::remove_if;
  using std::swap;
  using std::unique;
  using std::vector;
  using std::setw;
  using std::iter_swap;
  using std::random_access_iterator_tag;
  using std::bernoulli_distribution;
  using std::numeric_limits;
  typedef array<array<double,3>,3> Matrix3d;

  template <typename T>
  constexpr size_t digits(){
    return size_t(std::numeric_limits<T>::digits);
  }

  template <typename T>
  constexpr size_t popcount(const T n){
    size_t c=0;
    while(n) (n&=(n-1),++c);
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


  /*struct popcount{
    const uint64_t m1  = 0x5555555555555555; //binary: 0101...
    const uint64_t m2  = 0x3333333333333333; //binary: 00110011..
    const uint64_t m4  = 0x0f0f0f0f0f0f0f0f; //binary:  4 zeros,  4 ones ...
    const uint64_t m8  = 0x00ff00ff00ff00ff; //binary:  8 zeros,  8 ones ...
    const uint64_t m16 = 0x0000ffff0000ffff; //binary: 16 zeros, 16 ones ...
    const uint64_t m32 = 0x00000000ffffffff; //binary: 32 zeros, 32 ones
    const uint64_t hff = 0xffffffffffffffff; //binary: all ones
    const uint64_t h01 = 0x0101010101010101; //the sum of 256 to the power of
                                             //0,1,2,3...
    // the last line is obvously a bug TODO
    
    //This uses fewer arithmetic operations than any other known  
    //implementation on machines with slow multiplication.
    //It uses 17 arithmetic operations.
    uint64_t operator()(uint64_t x) {
      x -= (x >> 1) & m1;             //put count of each 2 bits into those 2 bits
      x = (x & m2) + ((x >> 2) & m2); //put count of each 4 bits into those 4 bits 
      x = (x + (x >> 4)) & m4;        //put count of each 8 bits into those 8 bits 
      x += x >>  8;  //put count of each 16 bits into their lowest 8 bits
      x += x >> 16;  //put count of each 32 bits into their lowest 8 bits
      x += x >> 32;  //put count of each 64 bits into their lowest 8 bits
      return x & 0x7f;
    }*/
    /*
    //This uses fewer arithmetic operations than any other known  
    //implementation on machines with fast multiplication.
    //It uses 12 arithmetic operations, one of which is a multiply.
    int popcount_3(uint64_t x) {
      x -= (x >> 1) & m1;             //put count of each 2 bits into those 2 bits
      x = (x & m2) + ((x >> 2) & m2); //put count of each 4 bits into those 4 bits 
      x = (x + (x >> 4)) & m4;        //put count of each 8 bits into those 8 bits 
      return (x * h01)>>56;  //returns left 8 bits of x + (x<<8) + (x<<16) + (x<<24) + ... 
    }*/
  /*};*/
  /*
  struct bitstream{
    uintmax_t digits;
    uintmax_t* data;
    integer& push_back(const integer& i){
      return this;
    }
    integer& push_front(const integer& i){
      return this;
    }
    integer& push_back(const bool b){
      return this;
    }
    integer& push_back(const bool b,size_t n){
      return this;
    }
    integer& push_front(const bool b){
      return this;
    }
    integer& push_front(const bool b,size_t n){
      return this;
    }
    integer& push_back(const uint8_t c){
      if ((log2+8)/numeric_limits<uintmax_t>::digits>log2/numeric_limits<uintmax_t>){
        data = realloc(data,
          (log2+8)/numeric_limits<uintmax_t>::digits*numeric_limits<uintmax_t>digits/8);
        if (!data){
          exit(EXIT_FAILURE);
        }
      }
      return this;
    }
    integer& push_front(const uint8_t c){
      return this;
    }
    integer& push_front(const size_t v){
      return this;
    }
    integer& push_front(const size_t v,const size_t n){
      return this;
    }
    integer& push_back(const size_t v){
      return this;
    }
    bitstream& push_back(const uintmax_t v,const uintmax_t n){
      return this;
    }
    uintmax_t uintmax_subset(const size_t i,const size_t j){

    }
    bool operator[](const size_t i){
      return (data[i/numeric_limits<uintmax_t>::digits]
             >>i%numeric_limits<uintmax_t>::digits)&uintmax_t(1);
    }
    bitstream(){
      log2=0;
      data = new uintmax_t[];
    }
    ~bitstream(){
      log2=0;
      delete[] data;
      data=0;
    }
  }

  // i.push_back(0,log2(log2(n)+1));
  // i.push_back(log2(n)+1,log2(log2(n)+1));
  // i.push_back(n,log2(n));
  template<typename T>
  typename std::enable_if<std::is_unsigned<T>::value,bitstream>::type
  const elias_delta_encode(const T n){
    const uintmax_t len = log2(n);
    const uintmax_t lol = log2(len+1);
    bitstream b;
    b.push_back(0,lol);
    b.push_back(len+1,lol);
    b.push_back(n,len);
    return b;
  }  
  // len = i<<clz(i)&(1<<(clz(i)+1)-1);
  // n   = i<<(clz(i)*2)
  const uintmax_t elias_delta_decode_as_uintmax_t(const bitstream& b){
    const uintmax_t n = 1;
    const uintmax_t len = 1;
    const uintmax_t lol = 0;
    size_t i;
    for(i=0;!b[i];++i){
      ++lol;
    }
    for(;i<lol;
    return 0;
  }
*/
  //typedef for a fraction
  template<typename T>
  using fraction = array<std::enable_if<std::is_integral<T>::value,T>,2>;
 
  /* only if cpu branches slowly, maybe there should be a constexpr version
  //fast way to calculate the log2 rounded down of a uint32_t
  const uint32_t log2(const uint32_t& n){
    uint32_t shift;
    uint32_t r;
    uint32_t v=n;
    r =     (v > 0xFFFF) << 4; v >>= r;
    shift = (v > 0xFF  ) << 3; v >>= shift; r |= shift;
    shift = (v > 0xF   ) << 2; v >>= shift; r |= shift;
    shift = (v > 0x3   ) << 1; v >>= shift; r |= shift;
                                            r |= (v >> 1);
    return r;
  }
  
  //fast way to calculate the log2 rounded down of a uint64_t
  const uint64_t log2(const uint64_t& n){
    uint64_t shift;
    uint64_t r;
    uint64_t v=n;
    r =     (v > 0xFFFFFFFF) << 5; v >>= r;
    shift = (v > 0xFFFF    ) << 4; v >>= shift; r |= shift;
    shift = (v > 0xFF      ) << 3; v >>= shift; r |= shift;
    shift = (v > 0xF       ) << 2; v >>= shift; r |= shift;
    shift = (v > 0x3       ) << 1; v >>= shift; r |= shift;
                                                r |= (v >> 1);
    return r;
  }
  */

  template <typename T>
  typename std::enable_if<std::is_unsigned<T>::value,T>::type
  constexpr ror(const T n, const T i){
    const T m = (std::numeric_limits<T>::digits-1);
    const T c = i&m;
    return (n>>c)|(n<<((-c)&m));
  }

  template <typename T>
  typename std::enable_if<std::is_unsigned<T>::value,T>::type
  constexpr rol(const T n, const T i){
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
/* repeated sqaring and testing if > 2 could potentially give a way to get
 * the binary logarithm of a number
 * uint64_t x;
 * uint64_t r;
 * for (size_t i=0; i<64; ++i){
 *   x*=x;
 *   if (x<1<<32) continue;
 *   r^=1<<(63-i);
 *   x>>=32;
 * }
 * can it be inverted?
 *
 * 0000 -> 
 * 0001 -> 0000
 * 0010 -> 0100
 * 0011 -> 0110 
 * 0100 -> 1000
 * 0101 -> 1001
 * 0110 -> 1010
 * 0111 -> 1011
 * 1000 -> 1100
 * 1001 -> 1100
 * 1010 -> 1101
 * 1011 -> 1101
 * 1100 -> 1110
 * 1101 -> 1110
 * 1110 -> 1111
 * 1111 -> 1111
 *  
 *
  template <typename T>
  array<typename std::enable_if<std::is_unsigned<T>::value,T>::type,2>
  constexpr log2_(const T& x,const T& lower=0,const T& upper=digits<T>()){
    const T floor_log2 = floor_log2(x);
    return {{floor_log2,
  }
*/

  template <typename T>
  typename std::enable_if<std::is_unsigned<T>::value,T>::type
  constexpr bswap(const T n){
    return 0;// This is wrong, TODO
  }
#if defined __builtin_bswap16
  constexpr uint16_t bswap(const uint16_t n){
    return __builtin_bswap16(n);// This is wrong
  }
#endif
#if defined __builtin_bswap32
  constexpr uint32_t bswap(const uint32_t n){
    return __builtin_bswap32(n);// This is wrong
  }
#else
  constexpr uint32_t bswap(const uint32_t n) // TODO generalize
  {
    uint32_t
    x = (((x & 0xaaaaaaaa) >> 1) | ((x & 0x55555555) << 1));
    x = (((x & 0xcccccccc) >> 2) | ((x & 0x33333333) << 2));
    x = (((x & 0xf0f0f0f0) >> 4) | ((x & 0x0f0f0f0f) << 4));
    x = (((x & 0xff00ff00) >> 8) | ((x & 0x00ff00ff) << 8));
    return((x >> 16) | (x << 16));
  }
#endif
#if defined __builtin_bswap64
  constexpr uint64_t bswap(const uint64_t n){
    return __builtin_bswap64(n);// This is wrong
  }
#endif
  template <typename T>
  typename std::enable_if<std::is_unsigned<T>::value,T>::type
  constexpr clz(const T x,const T lower=0,const T upper=digits<T>()){
    return (upper-lower==T(1))?digits<T>()-upper:(x&(T(0)-T(1)<<((upper+lower)/2))?
           clz(x,(upper+lower)/2,upper):
           clz(x,lower,(upper+lower)/2));
  }
#if defined __builtin_clzl
  uint32_t constexpr clz(const uint32_t x){
    return __builtin_clzl(x);
  }
#endif 
#if defined __builtin_clzll
  uint64_t constexpr clz(const uint64_t x){
    return __builtin_clzl(x);
  }
#endif 
  template <typename T>
  typename std::enable_if<std::is_unsigned<T>::value,T>::type
  constexpr log2(const T x,const T lower=0,const T upper=digits<T>()){
    return (upper-lower==T(1))?lower:(x&(T(0)-T(1)<<((upper+lower)/2))?
           log2(x,(upper+lower)/2,upper):
           log2(x,lower,(upper+lower)/2));
  }
#if defined __builtin_clzl
  uint32_t constexpr log2(const uint32_t x){
    return x==0?0:31-__builtin_clzl(x);
  }
#endif
#if defined __builtin_clzll
  uint64_t constexpr log2(const uint64_t x){
    return x==0?0:63-__builtin_clzll(x);
  }
#endif
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

  const inline uint64_t length_of_universal_code(uint64_t n){
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
  class bitstream{//please be little endian...
    private:
      deque<intmax_t> data;
      size_t digits=0;
      size_t foffs=0;
      size_t boffs=0;
    public:
      void push_back(const bool b){ //msb
      }
      void push_back(const intmax_t b,const size_t n){
      }
      void push_front(const bool b){ //lsb
      }
      void push_front(const intmax_t b,const size_t n){
      }
      void align(){
      }
  }
  
  class bigint{
    private:
      int digits;
      intmax_t * data;
    public:
      bigint(){
        data=new intmax_t[1];
        digits=0;
      }
      bigint& operator+=(const bigint o){

      }
      bigint& operator+=(const intmax_t o){

      }
      bigint& operator-=(const bigint o){

      }
      bigint& operator-=(const intmax_t o){

      }
      bigint& operator*=(const bigint o){

      }
      bigint& operator*=(const intmax_t o){

      }
      bigint& operator<<=(const size_t n){

      }
      bigint& operator>>=(const size_t n){

      }
      ~bigint(){
        digits=0;
        delete[] data;
        data=0;
      }
  }
*/
  /* todo, mantissa has one 1 too much...
  const inline uint16_t minifloat_encode(double v){
    const int l = v!=0?log2(abs(v)):0;
    uint16_t r;
    r = (v<0)<<15; // sign bit
    r|= l<<10&((~uint16_t(0))>>1);// exponent;
    double m = abs(v)*pow(2,l);
    r|= ((m>255.0)?(~uint16_t(0)):uint16_t(m+0.5))&(~uint16_t(0)>>6);
    return r;
  }

  const inline uint16_t minifloat_decode(uint16_t v){
    double r = 
  }
  */

  const unsigned int b[] = {0x2, 0xC, 0xF0, 0xFF00, 0xFFFF0000};
  const unsigned int S[] = {1, 2, 4, 8, 16};
  const uint32_t log2_v2(uint32_t v){
    unsigned int r = 0; // result of log2(v) will go here
    for (uint32_t i = 4; i!=std::numeric_limits<uint32_t>::max(); --i){
      // unroll for speed...
      if (v & b[i]){
        v >>= S[i];
        r |= S[i];
      } 
    }
    return r;
  }
  

  //Zigzag decoding and encoding for uint32_t and uint64_t
  //encoding:
  //0 -> 0; -1 -> 1; 1 -> 2; -2 ->  3; 2 -> 4 ...
  //decoding:
  //0 -> 0;  1 ->-1; 2 -> 1;  3 -> -2; 4 -> 2 ...
  template<typename T>
  typename std::make_unsigned<typename std::enable_if<std::is_integral<T>::value,T>::type>::type
  constexpr zigzag_encode(const T n){
    return (n<<1)^(n>>(std::numeric_limits<T>::digits-1));
  }

  template<typename T>
  typename std::make_signed<typename std::enable_if<std::is_unsigned<T>::value,T>::type>::type
  constexpr zigzag_decode(const T n){
    return (n>>1)^(-(n&1));
  }

  // TODO: Somehow is broken for negative integers
  // Function to determine lowest common denominator
  // This is a generalization of Euclids algorithm to a set of numbers
  template <typename T, size_t N>
  typename std::enable_if<std::is_integral<T>::value,T>::type
  const inline gcd(array<T,N> a){
    //for (auto it=a.begin(); it!=a.end();++it){
    //  *it=abs(*it);
    //}
    auto end=a.end();
    while(true){
      T num=0;
      auto newend=a.begin();
      for (auto it=a.begin();it!=end;++it){
        if (*it==0) continue;
        if (it!=newend) *newend=*it;
        ++newend;
      }
      end=newend;
      auto min=min_element(
          a.begin(),
          end,
          [](const T&x,const T&y){return abs(x)<abs(y);}
          ); 
      for (auto it=a.begin();it!=end;++it){
        if (it==min) continue;
        if (*it==0) continue;
        *it%=*min;
        ++num;
      }
      if (num==0) return *min;
    }
    return 0;
  }

  // is faster when a<b ...
  template <typename T>
  typename std::enable_if<std::is_unsigned<T>::value,T>::type
  constexpr gcd(const T& a,const T&b){
    return a?gcd(b%a,a):b;
  }

  template <class It>
  typename iterator_traits<It>::value_type
  const gcd(It begin, It end){
    size_t num;
    while (true){
      num=0;
      end=remove_if(begin,end,
          [](const typename iterator_traits<It>::value_type& x){return x==0;});
      auto min=min_element(begin,end);
      for (auto it=begin;it!=end;++it){
        if (it==min) continue;
        *it%=*min;
        ++num;
      }
      if (num==0) return *min;
    }
  }

  template <typename T>
  typename std::enable_if<std::is_integral<T>::value,T>::type
  constexpr lcm(const T& a, const T& b){
    return a/gcd(a,b)*b;
  }

  // TODO least common multiple, may not be needed for NTT
  template <class It>
  typename iterator_traits<It>::value_type
  constexpr inline lcm(It begin,
                       It end){
    return accumulate(begin,end,
                      typename iterator_traits<It>::value_type(1),
                      lcm<typename iterator_traits<It>::value_type>);
  }

  
  template<typename T>
  typename std::enable_if<std::is_integral<T>::value,T>::type
  constexpr count_trailing_zeroes(const T& n){
    if (!n) return std::numeric_limits<T>::digits();
    n=!n;
    T i(0);
    while (n&1){
      n>>=1;
      ++i;
    }
    return i;
  }
  
  template<typename T>
  typename std::enable_if<std::is_integral<T>::value,T>::type
  const factor_out_power_of_2(T& n){
    if (n==0) return 0;
    n=~n;
    T i(0);
    while (n&1){
      n>>=1;
      ++i;
    }
    n=~n;
    return i;
  }
  
  template<typename T>
  typename std::enable_if<std::is_signed<T>::value,T>::type
  const inline jacobi_symbol(T a, T n){
    if (n%2==0) return 0; // not defined for even n
    if (n<1) return 0;    // not defined for negative n
    while(a<0) a+=n;      // for negative a
    T sign = 1;
    while (true){
      a=a%n;
      T f = factor_out_power_of_2(a);
      if (f&1){ // if uneven factors
        switch (n%8){
          case 3:
          case 5:
            sign=-sign;
            break;
        }
      }
      if (a==1) return sign;
      if (a==0) return 0;
      if (gcd(a,n)!=1) return 0;
      swap(a,n);
      if ((a%4==3)&&(n%4==3)) sign=-sign;
    }
  }

  template<typename T>
  typename std::enable_if<std::is_integral<T>::value>::type
  const inline normalize(fraction<T>& x){
    const T g = gcd<T,2>(x);
    x[0]/=g;
    x[0]/=g;
  }

  template<typename T>
  typename std::enable_if<std::is_integral<T>::value,T>::type 
  const inline div_to_nearest(const T& n,const T& d) {
    return (n+(((n>0)^(d<0))?(d/2):(-d/2)))/d;
  }
  
  template<typename T>
  typename std::enable_if<std::is_integral<T>::value,T>::type
  const inline mod_to_nearest(const T& n,const T& d) {
    return (n+(((n>0)^(d<0))?(d/2):(-d/2)))%d+(((n>0)^(d<0))?(-d/2):(d/2));
  }

  template<typename T>
  typename std::enable_if<std::is_integral<T>::value>::type
  const inline div_mod_to_nearest(const T& n,const T& d, T& div, T& mod){
    const T t0 = ((n>0)^(d<0))?(d/2):(-d/2);
    const T t1 = n+t0;
    div = t1/d;
    mod = t1%d-t0;
  }

/* numerically stable sum of T (being double or float)
 * with ints this does not make sense and will not compile
 * start with
 * T c=0.0, sum=0.0;
 * and then call
 * kahan_sum(x,c,sum) for each x */
  template<typename T>
  typename std::enable_if<std::is_floating_point<T>::value>::type
  const inline kahan_sum(
      const T &x,
      T &c,
      T &sum
      ){
    const T y=x-c;
    const T t=sum+y;
    c = (t-sum)-y;
    sum=t;
  }

  template<typename T>
  typename std::enable_if<is_integral<T>::value>::type
  const inline average(
      const T& x,
      T& cum_avg,
      T& cum_rem,
      T& n
      )
  {
    ++n;
    const T addendum= x-cum_avg+cum_rem;
    T div;
    div_mod_to_nearest<T>(addendum,n,div,cum_rem);
    cum_avg+= div;
  }

  template<typename T>
  typename enable_if<is_integral<T>::value>::type
  const inline average(
      const T& x,
      const T& w,
      T& cum_avg,
      T& cum_rem,
      T& sumw
      )
  {
    sumw+=w;
    const T addendum= x*w-cum_avg*w+cum_rem;
    T div;
    div_mod_to_nearest<T>(addendum,sumw,div,cum_rem);
    cum_avg+=div;
  }

  template<typename T>                   //avg,sumw,rema
  typename enable_if<is_integral<T>::value,array<T,3>>::type
  const inline average(
      const array<T,3>& a, //essentially these are fractions
      const array<T,3>& b //with the integer part factored out
      )
  {
    // if (a[2]<b[2]) swap(a,b); could be done, is not necessarey, aye
    array<T,3> r;
    r[1] = a[1]+b[1]; // sumw
    const T diff = a[0]-b[0];
    //avg  = a[0]+((a[0]-b[0])*b[2]+a[3]+b[3])/r[2]; // *
    //avg  = b[0]+((b[0]-a[0])*a[2]+a[3]+b[3])/r[2];
    //avg  = (a[0]*a[2]+a[3]+b[0]*b[2]+b[3])/r[2];
    //avg  = (a[0]*a[2]+b[0]*b[2])/r[2]+(a[3]+b[3])/r[2];
    div_mod_to_nearest((a[0]-b[0])*b[1]+a[2]+b[2],r[1],r[0],r[2]);
    r[0]+=a[0];
    return r;
  }
  
/* power of x to e modulo m via exponentiation by squaring
 */ 
  template<typename T>
  typename std::enable_if<std::is_integral<T>::value,T>::type
  constexpr power(T x,T e){
    T result(1);
    while (e){
      if (e&1ull) result=(x*result);
      x=(x*x);
      e>>=1;
    }
    return result;
  }

/* power of x to e modulo m via exponentiation by squaring
 */ 
  template<typename T>
  typename std::enable_if<std::is_integral<T>::value,T>::type
  power_mod(T x,T e,const T& m){
    T result(1);
    while (e){
      if (e&1ull) result=(x*result)%m;
      x=(x*x)%m;
      e>>=1;
    }
    return result;
  }
  
  template<typename T>
  typename std::enable_if<std::is_integral<T>::value,T>::type
  power_mod_inplace(T& x,T& e,const T& m){
    T result(1);
    while (e){
      if (e&1ull) result=(x*result)%m;
      x=(x*x)%m;
      e>>=1;
    }
    return result;
  }
  
/* floor of square root of n
 */
  template<typename T>
  typename std::enable_if<std::is_integral<T>::value,T>::type
  const floor_square_root(const T& n){
    if (n==0) return 0;
    T l = log2(n)/T(2);
    T t0 = ((n>>l)+(T(1)<<l))/T(2);
    T t1,t2; 
    while (true){
      t2 = t1;
      t1 = t0;
      t0 = (t0+n/t0)/T(2);// Newton Heron
      //cout << t0 << endl;
      if (t0==t1) return t0;
      if (t0==t2) return t0<t1?t0:t1;
    }
    return 0;
  }

  /*template<int dim>
  const inline void mean_variance(
      const Eigen::Matrix<double,dim,  1>& x,
      const Eigen::Matrix<double,dim,dim>& w,
            Eigen::Matrix<double,dim,dim>& sumw,
            Eigen::Matrix<double,dim,  1>& mean,
            Eigen::Matrix<double,dim,dim>& M2
      ){
    const Eigen::Matrix<double,dim,dim> temp  = w+sumw;
    const Eigen::Matrix<double,dim,  1> delta = x-mean;
    const Eigen::Matrix<double,dim,  1> R     = (w*temp.inverse())*delta;
    mean += R;
    M2   += (sumw*delta)*R.transpose();
    sumw = temp;
  }*/
/*
  template<size_t d>
  const inline void mean_variance(
      const dlib::matrix<double,dim,  1>& x,
      const dlib::Matrix<double,dim,dim>& w,
            dlib::Matrix<double,dim,dim>& sumw,
            dlib::Matrix<double,dim,  1>& mean,
            dlib::Matrix<double,dim,dim>& M2
      ){
    const Eigen::Matrix<double,dim,dim> temp  = w+sumw;
    const Eigen::Matrix<double,dim,  1> delta = x-mean;
    const Eigen::Matrix<double,dim,  1> R     = (w*temp.inverse())*delta;
    mean += R;
    M2   += (sumw*delta)*R.transpose();
    sumw = temp;
  }

  template<typename T1,typename T2, typename T3, typename T4, typename T5>
  const inline void mean_variance(
      const dlib::matrix_exp<T1>& x,
      const dlib::matrix_exp<T2>& w,
            dlib::matrix_exp<T3>& sumw,
            dlib::matrix_exp<T4>& mean,
            dlib::matrix_exp<T5>& M2
      ){
    M2+=((x-mean)*dlib::trans(x-mean))*(w*sumw*inv(w+sumw));
    //M2  +=(sumw*(x-mean))*((w*dlib::inv(w+sumw))*dlib::inv((x-mean)));
    mean+=(w*inv(w+sumw))*(x-mean);
    sumw+=w;
  }
*/
/* numerically stable and incremental mean and variance
 * start with T sumw=0.0, mean=0.0, M2=0.0;
 * or with T sumw=w_1, mean=x_1, M2=?;
 * and then call
 * mean_variance(x,w,sumw,mean,M2)
 * for each pair of x and w */
  template<typename T>
  typename enable_if<is_floating_point<T>::value>::type
  const inline mean_variance(
      const T &x,
      const T &w,
      T &sumw,
      T &mean,
      T &M2
      ){
    // if weight too small then what? It will be your fault if it breaks...
    // if (w<numeric_limits<T>::epsilon()) return;
    const T temp = w+sumw;
    const T delta = x-mean;
    const T R = delta*w/temp;
    mean += R;
    M2   += sumw*delta*R; // sumw*(x-mean)*(x-mean)*w/(sumw+w)
    sumw  = temp;
  }



/* non corrected sample variance, correction factor is n*(n-1) 
 * to be used in conjunction with mean_variance */
  template<typename T>
  inline const T variance(
      const T& M2,
      const T& sumw
      ){
    return M2/sumw;
  }


/* non corrected sample variance, correction factor is n*(n-1)
 * to be used in conjunction with mean_variance
 * n being the number of x for which mean_variance was called.*/
  template<typename T>
  inline const T variance(
      const T &M2,
      const T &sumw,
      const size_t &n
      ){
    return M2*n/(sumw*(n-1)); //TODO: or do I want M2/(sumw*(n-1)) ?? 
  }

  template<typename T>
  struct mean_variance_calculator{
    T M2=0;
    T sumw=0;
    size_t n=0;
    T mean=0;
    void push(const T& x, const T& w){
      ++n;
      mean_variance(x,w,sumw,mean,M2);
    } 
    T variance() const {
      return wmath::variance(M2,sumw,n);
    }
  };

  template <class InputIterator>
  typename iterator_traits<InputIterator>::value_type
  inline const mean_variance(InputIterator beg,InputIterator end){
    typedef typename iterator_traits<InputIterator>::value_type::value_type V;
    array<V,2> mv;
    if (distance(beg,end)==1){
      mv = *beg;
      mv[1]=1/mv[1];
      return mv;
    }
    V sumw=0;
    for (InputIterator it=beg;it!=end;++it){
      mean_variance((*it)[0],(*it)[1],sumw,mv[0],mv[1]);
    }
    mv[1]/=sumw*(distance(beg,end)-1);
    return mv;
  }

  template <class InputIterator>
  typename iterator_traits<InputIterator>::value_type
  inline const ransack_mean_variance(InputIterator begin,InputIterator end){
    using std::bernoulli_distribution;
    typedef typename iterator_traits<InputIterator>::value_type V;
    std::mt19937 mt=std::mt19937{std::random_device{}()};
    bernoulli_distribution distr(0.5);
    const size_t k=1<<12;
    V sumw,mean,M2,n;
    array<array<V,2>,k> samples;
    for (size_t i=0;i<k;){
      sumw=mean=M2=n=0.0;
      for (auto it=begin;it!=end;++it){
        if (distr(mt)) continue;
        mean_variance(
            (*it)[0],
            (*it)[1],
            sumw,
            mean,
            M2
            );
        ++n;
      }
      if (n>0.5){
        samples[i][0]=mean;
        samples[i][1]=n;
        ++i;
      }
    }
    sumw=mean=M2=0.0;
    for (size_t i=0;i<k;++i){
      mean_variance(
          samples[i][0],
          samples[i][1],
          sumw,
          mean,
          M2
          );
    }
    return {{mean,sqrt(0.5*M2/sumw)}};
  }

  template<class RandomAccessIterator>
  inline const void gnome_sort(RandomAccessIterator beg, RandomAccessIterator end){
    if (distance(beg,end)<2) return;
    for (RandomAccessIterator it=beg+1; it!=end; ++it)
      for (RandomAccessIterator pos = it;(pos!=beg)&&*(pos-1)>*(pos);--pos) iter_swap((pos-1),pos);
  }

/*
  template <class RandomAccessIterator>
  inline const void inplace_merge(RandomAccessIterator beg, RandomAccessIterator mid, RandomAccessIterator end){
    typedef typename std::iterator_traits<RandomAccessIterator>::iterator_category category;
    return __inplace_merge(beg, mid, end, category());
  }

  // broken, needs heap for elements [mid,up)
  template <class RandomAccessIterator>
  inline const void __inplace_merge(RandomAccessIterator beg, RandomAccessIterator mid, RandomAccessIterator end, random_access_iterator_tag){
    for (auto it=beg;it!=end; ++it){
      cerr << *it << " ";
    }
    cerr << endl;
    RandomAccessIterator it=beg;
    RandomAccessIterator up=mid;
    while (*it<*up){
      ++it;
      if (it==end) return;
    }
    iter_swap(it,up);
    ++it;
    // ## ## ## it ## ## ## ## ## -- -- qp -- qe -- -- up ## ## ## ##
    RandomAccessIterator qp=up,qe=up;
    ++up;
    while (it!=end){
      for (auto it=beg;it!=end; ++it){
        cerr << setw(2) << distance(beg,it) << " ";
      }
      cerr << endl;
      for (auto it=beg;it!=end; ++it){
        cerr << setw(2) << *it << " ";
      }
      cerr << "   it=" << distance(beg,it) << "  qp=" << distance(beg,qp) << "  up=" << distance(beg,up);
  }
      if (*up<*qp){
        iter_swap(up,it);
        if(qp==it) qp=up;
        ++up;
      }else{
        iter_swap(qp,it);
        ++qp;
        if(qp==qe){
          qp=max(it,mid);
          qe=up;
          cerr << "   wrap";
        }
      }
      ++it;
      if (qp<it){
        ++qp;
        cerr << "  ++qp ";
      }
      cerr << endl << endl;
    }
  }

  template <class RandomAccessIterator>
  inline const void inplace_mergesort(RandomAccessIterator beg, RandomAccessIterator end){
    typedef typename std::iterator_traits<RandomAccessIterator>::iterator_category category;
    if (distance(beg,end)<9) return gnome_sort(beg,end);
    RandomAccessIterator mid=beg;
    advance(mid,distance(beg,end)/2);
    inplace_mergesort(beg,mid);
    inplace_mergesort(mid,end);
    __inplace_merge<RandomAccessIterator>(beg,mid,end,category());
    beg+=3;
  }  

  template<class RandomAccessIterator>
  inline const void gnome_sort_strided(RandomAccessIterator beg, RandomAccessIterator end,size_t stride){
   // cerr << "strided gnomesort begin " << stride << endl;
    if (distance(beg,end)<2) return;
    for (RandomAccessIterator it=beg+stride; it<end; it+=stride)
      for (RandomAccessIterator pos = it;(pos!=beg)&&*(pos-1)>*(pos);pos-=stride) iter_swap((pos-1),pos);
   // cerr << "strided gnomesort end" << endl;
  }

  template<class RandomAccessIterator>
  inline const void __merge_sort_strided_helper(RandomAccessIterator fir, RandomAccessIterator las,size_t stride){
   // cerr << "helper" << endl;
   // cerr << distance(fir,las) << endl;
    las=fir+distance(fir,las)/(2*stride)*(2*stride);
   // cerr << distance(fir,las) << endl;
    fir+=stride;
    while (fir<las){
      if (*las<*fir) iter_swap(fir,las);
      fir+=stride;
      las-=stride;
    }
  //  cerr << "helper end" << endl;
  }
  // still todo
  template <class RandomAccessIterator>
  inline const void __inplace_merge_strided(RandomAccessIterator beg, RandomAccessIterator end, size_t stride, random_access_iterator_tag){
    RandomAccessIterator heap_beg=beg; 
    RandomAccessIterator heap_end=end;
    RandomAccessIterator it0=beg;
    RandomAccessIterator it1=beg+stride;
    while (it0<end){
      //compare it1 with min element of heap
      //if less then iter_swap it0 it1 and it1+=stride
      //if larger then pop_heap and continue  
      it0+=stride;
    }
  }

  template <class RandomAccessIterator>
  const void inplace_merge(RandomAccessIterator beg,RandomAccessIterator mid, RandomAccessIterator end){
    RandomAccessIterator it0=beg;
    RandomAccessIterator it1=mid;
    RandomAccessIterator it2=mid;
    while (it1!=end){
      if (*it1<it0){
        iter_swap(it1,it0);
        break;
      }
      ++it0;
    }
    while (it0!=end){
      if (*it2<*it1){
        iter_swap(it0,it2);
        ++it2;
        // range(it1,it2) is still heap as new last element is max element
      }else{
        iter_swap(it0,it1);
        heap_push_max_front(it1,it2,it1); // O(log(n)) moves
      }
      ++it0;
    }
  }

  template <class RandomAccessIterator>
  const void inplace_mergesort(RandomAccessIterator beg, RandomAccessIterator end, random_access_iterator_tag){
    if (distance(beg,end)<16) return gnome_sort(beg,end);
    const RandomAccessIterator mid=beg;
    advance(mid,(distance(beg,end)+1)/2);
    inplace_mergesort(beg,mid);
    inplace_mergesort(beg,mid);
    inplace_merge(beg,end);
  }


  template <class RandomAccessIterator>
  inline const void __inplace_mergesort_strided(RandomAccessIterator beg, RandomAccessIterator end,const size_t stride,random_access_iterator_tag){
    typedef typename std::iterator_traits<RandomAccessIterator>::iterator_category category;
  //  cerr << "anfang stride = " << stride<< endl;
    if (stride>distance(beg,end)) return;
    if (stride>distance(beg,end)/2){
  //    cerr << "test " << endl;
      RandomAccessIterator it0 = beg;
      RandomAccessIterator it1 = beg+stride;
      while (it1!=end){
        if (*it1<*it0) iter_swap(it1,it0);
        ++it1,++it0;
      }
    //  cerr << "this return ... " << endl;
      return;
    }
    __inplace_mergesort_strided(beg,end,stride<<1,category());
    //if (stride==1) return gnome_sort(beg,end);
    for (size_t i=0; i<stride; ++i){
      //__merge_sort_strided_helper(beg+i,end,stride);
      //gnome_sort_strided(beg+i,end,stride);
      __inplace_merge_strided(beg+i,end,stride,category());
    }
  //  cerr << "ende" << endl;
  }

  template <class RandomAccessIterator>
  inline const void inplace_mergesort_strided(RandomAccessIterator beg, RandomAccessIterator end){
    typedef typename std::iterator_traits<RandomAccessIterator>::iterator_category category;
    __inplace_mergesort_strided<RandomAccessIterator>(beg,end,1,category());
  }

  // still todo
  template <class RandomAccessIterator>
  inline const void partition(RandomAccessIterator beg, RandomAccessIterator end, RandomAccessIterator pivot){ 
    --end;
sfup0:
    if(*pivot<*beg) goto sfdw1;
    ++beg;
    if(beg==end) goto end;
sfdw0:
    if(*end<*pivot) goto sfup1;
    --end;
    if(beg==end) goto end;
    goto sfup0;
sfup1:
    if(*pivot<*beg) goto swap;
    if(distance(beg,end)==1){
      mid = end+1;//todo
      goto push;
    }
    ++beg;
    goto sfup1;
sfdw1:
    if(*end<*pivot) goto swap;
    if(distance(beg,end)==1){
      mid = beg;//todo
      goto push;
    }
    --end;
    goto sfdw1;
swap:
    iter_swap(beg,end);
    ++beg;
    if(beg!=end){
      --end;
      if (beg!=end){
        goto sfup0;
      }
    }
end:
    if (pivot<beg){
      while (*pivot<*beg){
       --beg;
      }
      iter_swap(beg,pivot);
    }else{
      while(*beg<*pivot){
        ++beg;
      }
      iter_swap(beg,pivot);
    }
    return pivot;
  */

  // Normalize weights so that their average is 1.0
  template <class InputIterator> 
  inline const void normalize_av_weight_to_1(InputIterator begin,InputIterator end){
    double mean=0,sumw=0,M2=0;
    for(InputIterator it=begin; it!=end; ++it){
      mean_variance(
          (*it)[1],
          1.0,
          sumw,
          mean,
          M2
          );
    }
    for (InputIterator it=begin; it!=end; ++it){
      (*it)[1]/=mean;
    }
  }

  // Assume exponential distribution on weights to correct them, set c to 1/n, n being sample size
  template <class InputIterator> 
  inline const void correct_weights(InputIterator begin,InputIterator end,const double c){ 
    double mean=0,sumw=0,M2=0;
    for(InputIterator it=begin; it!=end; ++it){
      mean_variance(
          (*it)[1],
          1.0,
          sumw,
          mean,
          M2
          );
    }
    for (InputIterator it=begin; it!=end; ++it){
      const double x = (*it)[1]/mean;
      (*it)[1]=mean*(c+x*exp(-x))/(c+x*exp(-x));
    }
  }

  template<class RandomAccessIterator>
  inline const void sort2(RandomAccessIterator fir, RandomAccessIterator las){
    if (*fir>*las) swap(*fir,*las);
  }

  template<class RandomAccessIterator>
  inline const void sort3(RandomAccessIterator fir, RandomAccessIterator las){ 
    const RandomAccessIterator sec=fir+1;
    if (*fir>*sec) swap(*fir,*sec);
    if (*sec<*las) return;
    swap(*sec,*las);
    if (*fir>*sec) swap(*fir,*sec);
  }


  template<class RandomAccessIterator>
  inline const void sort4(RandomAccessIterator fir, RandomAccessIterator las){
    const RandomAccessIterator sec=fir+1;
    const RandomAccessIterator thr=las-1; 
    if (*fir>*sec) swap(*fir,*sec);
    if (*thr>*las) swap(*thr,*las);
    if (*fir>*thr){
      swap(*fir,*thr);
      if (*thr>*las){
        swap(*sec,*thr);
      }else{
        swap(*sec,*thr);
      }
      return;
    }else{
      if (*sec>*thr){
        swap(*sec,*thr);
        if (*thr>*las) swap(*thr,*las);
      }
    }
  }



  template<class RandomAccessIterator>
  inline const void bubblesort(RandomAccessIterator beg, RandomAccessIterator end){
    for (RandomAccessIterator it0=end-1; it0!=beg; --it0){
      for (RandomAccessIterator it1=beg; it1!=it0; ++it1){
        if (*it1<*(it1+1)) swap(*it1,*(it1+1));
      }
    }
  }

  //todo heapsort, introsort
  /*
  template <class RandomAccessIterator>
  inline const void heapify(RandomAccessIterator beg, RandomAccessIterator end){
    RandomAccessIterator i = end;
    --i;
    do {
      //assert(isheap(H, left(i)) && isheap(H, right(i));
      RandomAccessIterator min = i; 
      if (left(i) < H.size && H.key(left(i)) < H.key(min))
      min = left(i)
      if (right(i) < H.size && H.key(right(i)) < H.key(min))
      min = right(i)
      if (min == i) break
      H.swap(i, min)
      i = min
    } while(true)
  }
  */

  //todo selber machen
  template <typename RandomAccessIterator>
  inline const void heapsort(RandomAccessIterator beg,RandomAccessIterator end){
    std::make_heap(beg,end);
    std::sort_heap(beg,end);
  }

  // first element is pivot, second is less than or equal, last is greater than or equal pivot
  template <typename RandomAccessIterator>
  inline const RandomAccessIterator __unguarded_partition(RandomAccessIterator beg,RandomAccessIterator end, RandomAccessIterator pivot){
    //++beg;
    //advance(end,-2);
    --end;
    while(true){
      while (*beg<*pivot) ++beg;
      while (*pivot<*end) --end;
      if (!(beg<end)) return beg;
      iter_swap(beg,end);
      --end;
      ++beg;
    }
  }

  template <typename RandomAccessIterator>
  inline const void move_median_to_first(RandomAccessIterator mov,RandomAccessIterator fir,RandomAccessIterator sec,RandomAccessIterator thr){
    if (*sec<*fir) swap(*fir,*sec);
    if (*thr<*sec){
      iter_swap(sec,thr);
      if (*fir>*sec) swap(*fir,*sec);
    }
    iter_swap(mov,sec);
  }


  template<class RandomAccessIterator>
  typename iterator_traits<RandomAccessIterator>::value_type::value_type // median
  inline const weighted_slow_median(RandomAccessIterator beg, RandomAccessIterator end){
    if (distance(beg,end)<3){
      if (distance(beg,end)==2){
        const RandomAccessIterator sec = beg+1;
        return ((*beg)[0]*(*beg)[1]+(*sec)[0]*(*sec)[1])/((*beg)[1]+(*beg)[1]);
      }
      if (beg==end) return 0.0/0.0;
      return (*beg)[0];
    }
    typedef typename iterator_traits<RandomAccessIterator>::value_type::value_type V;
    std::sort(beg,end,
        [](const typename iterator_traits<RandomAccessIterator>::value_type& a,
           const typename iterator_traits<RandomAccessIterator>::value_type& b) -> bool{
          return a[0]<b[0];});
    RandomAccessIterator upper=end-1;
    RandomAccessIterator lower=beg;
    V weight = 0;
    for (auto it=beg;it!=end;++it){
      weight+=(*it)[1];
    }
    V w0=(*lower)[1],w1=(*upper)[1];
    while (w0+(*(lower+1))[1]<weight/2){
      w0+=(*lower)[1];
      ++lower;
      if(lower+1<end) continue;
      break;
    }
    w0=weight/2-w0;
    w1=weight/2-w1;
    while (w1+(*(upper-1))[1]<weight/2){
      w1+=(*upper)[1];
      --upper;
      if(upper>beg) continue;
      break;
    }
    return (w1*(*lower)[0]+w0*(*upper)[1])/(w0+w1); 
  }

  template<class RandomAccessIterator>
  typename iterator_traits<RandomAccessIterator>::value_type::value_type // median 
  inline const weighted_quick_median(RandomAccessIterator beg, RandomAccessIterator end){
    //cerr << "weighted_quick_median, " << distance(beg,end) << endl;
    typedef typename iterator_traits<RandomAccessIterator>::value_type::value_type V;
    V w0=0;
    V w1=0;
    for (auto it=beg;it!=end;++it){
      w0+=(*it)[1];
    }
    //cerr << "ok" << endl;
    while(distance(beg,end)>6){
      //cerr << beg-beg << " " << 1 << " " << (end-beg)/2 << " " << end-beg-1 << endl;
      move_median_to_first(beg,beg+1,beg+(end-beg)/2,end-1);
      RandomAccessIterator it0=beg+2,it1=end-2;
      //cerr << "hä" << endl;
      while (it0<it1){
        while((*it0)[0]<(*beg)[0]) ++it0;
        while((*beg)[0]<(*it1)[0]) --it1;
        //cerr << "testing" << endl;
        //cerr << distance(it0,it1);
        iter_swap(it0,it1);
        ++it0;--it1;
      }
      //cerr << "lived up until here" << endl;
      iter_swap(beg,it1);
      for (auto it=beg;it!=it1;++it){
        w0+=(*it)[1];
      }
      //cerr << "hmm" << endl;
      //cerr << distance(beg,end) << " " << w0 << " " << w1 << endl;
      if (w0/2<w1){
        end=it1;
        w0=w1;
        continue;
      }
      if (w1<w0/2){
        beg=it1;
        w0-=w1;
        continue;
      }
      return (*it1)[0];
    }
    //cerr << "quickfind loop passed " << endl;
    if (distance(beg,end)<2) return (*beg)[0];
    for (auto it=beg+1; it!=end; ++it)
      for (auto pos = it;(pos!=beg)&&(*(pos-1))[0]>(*pos)[0];--pos)
        iter_swap((pos-1),pos);
    w1=0;
    for (auto it=beg;it!=end;++it){
      w1+=(*it)[1];
      if(w1>w0/2) return (*it)[0]; 
    }
    return (*beg)[0];
  }

  template <class RandomAccessIterator>
  typename iterator_traits<RandomAccessIterator>::value_type //&{{median,MAD}}
  inline const weighted_quick_mad_median(RandomAccessIterator beg, RandomAccessIterator end){
    typedef typename iterator_traits<RandomAccessIterator>::value_type::value_type V;
    V median=weighted_quick_median(beg,end);
    vector<array<V,2>> ad;
    ad.reserve(distance(beg,end));
    for (auto it=beg;it!=end;++it){
      ad.push_back({{median>(*it)[0]?median-(*it)[0]:(*it)[0]-median,(*it)[1]}});
    }
    V mad = 1.4826*weighted_quick_median(ad.begin(),ad.end());
    return {{median,mad}};
  }

  template <class RandomAccessIterator>
  typename iterator_traits<RandomAccessIterator>::value_type //&{{median,MAD}}
  inline const weighted_slow_mad_median(RandomAccessIterator beg, RandomAccessIterator end){
    typedef typename iterator_traits<RandomAccessIterator>::value_type::value_type V;
    V median=weighted_quick_median(beg,end);
    vector<array<V,2>> ad;
    ad.reserve(distance(beg,end));
    for (auto it=beg;it!=end;++it){
      ad.push_back({{median>(*it)[0]?median-(*it)[0]:(*it)[0]-median,(*it)[1]}});
    }
    V sumw = 0;
    for (auto it=beg;it!=end;++it){
      sumw+=(*it)[1]*(*it)[1];
    }
    V mad = 1.4826*weighted_slow_median(ad.begin(),ad.end());
    V mv  = 1/sumw;
    return {{median,sqrt((mad*mad+mv*mv)/distance(beg,end))}};
  }

  template <class RandomAccessIterator>
  typename iterator_traits<RandomAccessIterator>::value_type
  inline const ransack_weighted_median(RandomAccessIterator beg,RandomAccessIterator end){ 
    //cerr << distance(beg,end) << endl;
    typedef typename iterator_traits<RandomAccessIterator>::value_type::value_type V;
    if (distance(beg,end)==0) return {{1.0/0.0,numeric_limits<V>::max()}};
    if (distance(beg,end)==1) return {{(*beg)[0],sqrt(1/(*beg)[1])}};
    std::mt19937 mt=std::mt19937{std::random_device{}()};
    bernoulli_distribution distr(0.5);
    const size_t k=1<<12;
    array<array<V,2>,k> samples;
    vector<array<V,2>> subsamples(distance(beg,end));
    //cerr << "ransack_weighted_median is doing fine" << endl;
    V minsd=0;
    for(size_t j,i=0;i<k;){
      //cerr << "begin first loop" << endl;
      j=0;
      for (auto it=beg;it!=end;++it){
        if (distr(mt)){
          subsamples[j]=*it;
          ++j;
        }
      }
      //cerr << "samples found" << endl;
      if (j){
        samples[i] = weighted_slow_mad_median(&subsamples[0],&subsamples[j]);
        //samples[i][0]=weighted_slow_median(&subsamples[0],&subsamples[j]);
        minsd+=1/(samples[i][1]*samples[i][1]);
        samples[i][1]=j;
        ++i;
      }
    }
    //cerr << "first loop done" << k << endl;
    V sumw=0,mean=0,M2=0;
    for (size_t i=0;i<k;++i){
      mean_variance(
          samples[i][0],
          samples[i][1],
          sumw,
          mean,
          M2
          );
    }
    //cerr << "next: return" << endl;
    return array<V,2>{{mean,sqrt(0.5*M2/sumw+0.5/minsd)}};
  }

  template <typename RandomAccessIterator>
  inline const RandomAccessIterator __unguarded_partition_pivot(RandomAccessIterator beg,RandomAccessIterator end){
    RandomAccessIterator mid = beg;
    advance(mid,distance(beg,end)/2);
    move_median_to_first(beg,beg+1,mid,end-1);
    return __unguarded_partition(beg+1,end,beg);
    //return __unguarded_partition(beg,end,beg);
  }

  template <typename RandomAccessIterator>
  const void __introsort_loop(RandomAccessIterator beg,RandomAccessIterator end,size_t depth_limit){
    while(distance(beg,end)>16){
      if (depth_limit==0) return heapsort(beg,end);
      RandomAccessIterator cut = __unguarded_partition_pivot(beg,end);
      iter_swap(cut,beg);
      --depth_limit;
      const size_t d0 = distance(beg,cut);
      const size_t d1 = distance(cut+1,end);
      if (d1<d0){
        __introsort_loop(cut+1,end,depth_limit);
        end=cut;
      }else{
        __introsort_loop(beg,cut,depth_limit);
        beg=cut+1;
      }
    }
  }

  template <typename RandomAccessIterator>
  inline const void sort(RandomAccessIterator beg, RandomAccessIterator end){
    if (beg!=end){
      __introsort_loop(beg,end,2*log2(size_t(distance(beg,end))));
      //gnome_sort(beg,end);
    }
  }

  // quicksort seems ok
  template <class RandomAccessIterator>
  inline const void quicksort(RandomAccessIterator beg, RandomAccessIterator end){
    const size_t c = 17;
    //typedef typename iterator_traits<RandomAccessIterator>::value_type V;
    typename iterator_traits<RandomAccessIterator>::value_type pivot;
    RandomAccessIterator* stack;
    RandomAccessIterator mid;
    RandomAccessIterator itu;
    RandomAccessIterator itd;
    RandomAccessIterator BEG=beg,END=end;
    size_t d0,d1;
    size_t p = ~size_t(0);
    size_t _p=0;
    size_t level=2*log2(size_t(distance(beg,end)));
    //if (distance(beg,end)<c) goto end;
    stack = new RandomAccessIterator[2*(log2(size_t(distance(beg,end)))-5)];
    goto partition;
push:
//    if (VERBOSE) cerr << "push" << endl;
//    if (VERBOSE){ cerr << ">>";
//      for (auto it=beg;it!=end;++it){
//        if (it==mid) cerr << "|";
//        else cerr << " ";
//        cerr << setw(2) << *it;
//      }
//      cerr << endl;
//    }
    d0 = distance(beg,mid);
    d1 = distance(mid,end);
    //cerr << d0 << " " << d1 << endl;
    //if (d0==0||d1==0) return;
    if (d0<d1){
      if (d0<c){
//        if (VERBOSE) stupid_sort(beg,mid);
        if (d1<c){
//          if (VERBOSE) stupid_sort(mid,end);
          goto pop;
        }else{
          beg=mid;
          goto partition;
        }
      }else{
//        if (VERBOSE) cerr << "pushing " << distance(mid,end) << endl;
        ++p;
        stack[2*p]  =mid;
        stack[2*p+1]=end;
        end=mid;
      }
    }else{
      if (d1<c){
//        if (VERBOSE) stupid_sort(mid,end);
        if (d0<c){
//          if (VERBOSE) stupid_sort(beg,mid);
          goto pop;
        }else{
          end=mid;
          goto partition;
        }
      }else{
//        if (VERBOSE) cerr << "pushing " << distance(beg,mid) << endl;
        ++p;
        stack[2*p]  =beg;
        stack[2*p+1]=mid;
        beg=mid;
      }
    } 
    goto partition;
pop:
//    if (VERBOSE) cerr << "pop " << endl;
    if (p==~size_t(0)) goto end;
    beg=stack[2*p];
    end=stack[2*p+1];
    --p;
//    if (VERBOSE) cerr << "popped " << distance(beg,end) << endl; 
partition:
//    if(VERBOSE){cerr << "partition" << endl;
//    cerr << distance(beg,end) << endl;}
    mid=beg;
    advance(mid,(distance(beg,end))/2);
    itu = beg;
    itd = end-1;
    if (*beg>*mid) iter_swap(beg,mid);
    if (*mid>*itd){
      iter_swap(mid,itd);
      if (*beg>*mid) iter_swap(beg,mid);
    }
    pivot=*mid;
    //cerr << "pivot = " << pivot << endl;
    --itd;++itu;
//    if (VERBOSE){ cerr << "<< ";
//    for (auto it=beg;it!=end;++it){
//      cerr << *it << " ";
//    }
//    cerr << endl;
//    cerr << "itu = " << distance(beg,itu) << ", " << *itu << " ;  itd =" << distance(beg,itd) << ", " << *itd << endl;
//    }
    //umbrellaswap(beg,end);
sfup0:
//    if (VERBOSE) cerr << "sfup0" << endl;
    if(*itu>pivot) goto sfdw1;
    ++itu;
//    if (VERBOSE) cerr << "itu = " << distance(beg,itu) << endl;
    if(itu==itd){
      mid = itu;
      if (*mid<pivot) ++mid;
      //goto findmid;
      goto push;
    }
sfdw0:
//    if (VERBOSE) cerr << "sfdw0"  << endl;
    if(*itd<pivot) goto sfup1;
    --itd;
//    if (VERBOSE) cerr << "itd = " << distance(beg,itd) << endl;
    if(itu==itd){
      mid = itu;
      //goto findmid;
      if (*mid<pivot) ++mid;
      goto push;
    }
    goto sfup0;
sfup1:
//    if (VERBOSE) cerr << "sfup1" << endl;
    if(*itu>pivot) goto swap;
    if(distance(itu,itd)==1){
      mid = itd+1;
      //goto findmid;
      goto push;
    }
    ++itu;
//    if (VERBOSE) cerr << "itu = " << distance(beg,itu) << endl;
    goto sfup1;
sfdw1:
//    if (VERBOSE) cerr << "sfdw1" << endl;
    if(*itd<pivot) goto swap;
    if(distance(itu,itd)==1){
      mid = itu;
      //goto findmid;
      goto push;
    }
    --itd;
//    if (VERBOSE) cerr << "itd = " << distance(beg,itd) << endl;
    goto sfdw1;
swap:
//    if (VERBOSE) cerr << "swap(" << distance(beg,itu) << "," << distance(beg,itd) << ")(" << *itu << "," << *itd << ")" << endl; 
    iter_swap(itu,itd);
//    if (VERBOSE){ cerr << ":: ";
//    for (auto it=beg;it!=end;++it){
//      cerr << *it << " ";
//    }
//    cerr << endl;
//    }
    ++itu;
    if(itu!=itd){
      --itd;
      if (itu!=itd){
        goto sfup0;
      }
    }
//    if (VERBOSE) cerr << "findmid" << endl;
    mid=itu;
//    if (VERBOSE) cerr << "mid = " << distance(beg,mid) << ", " << *mid << endl;
    if (*mid<pivot) ++mid;
    //while (*mid<pivot&&mid!=end) ++mid; 
    //if (*mid>pivot&&mid!=beg) --mid;
    //d0 = distance(beg,mid);
    //d1 = distance(mid,end);
    //if (*mid==pivot){
    //  if (d0<d1) ++mid;
    //  else --mid;
    //}
    goto push;
end:
//    cerr << "end" << endl;
//    if (VERBOSE){ cerr << "end" << endl;
//    cerr << p << endl;
//    }
    delete[] stack;
    gnome_sort(BEG,END);
    return;
  }

  /*
TODO: numerically stable implementation
calculate all moments with compensated summation P.S.: comensated summation is stupid, the only truth is integer math

  template<typename T, size_t N>
  inl<F7>ine const void mean_moments(
      const T& x,
      const T& w,
      T& sumw,
      array<T,N> M)
  {

    const T temp = w+sumw;
    const T delta = x-mean;
    const T R = delta*w/temp;
    mean += R;
    M2 += sumw*delta*R;
    sumw = temp;
  }

  template<typename T,0>
  inline const void mean_moments(
      const T& x,
      const T& w,
      T& sumw,
      array<T,1> M)
  {
  }

  template<typename T,1>
  inline const void mean_moments(
      const T& x,
      const T& w,
      T& sumw,
      array<T,1> M)
  {
    mean+=(x-mean)*w/
  }

*/

  // Pearsons product moment correlation coefficient
  template<typename T>
  class cc_pearson{
    private:
      T M2_x   =0;
      T M2_y   =0;
      T M2_xy  =0;
    public:
      size_t n =0;
      T mean_x =0;
      T sumw_x =0;
      T mean_y =0;
      T sumw_y =0;
      T mean_xy=0;
      T sumw_xy=0;
      inline void push(const T& x, const T &y){
        push(x,static_cast<T>(1),y,static_cast<T>(1));
      }
      // TODO: test it weights really work, especially for covariance
      inline void push(const T& x, const T& wx, const T& y, const T& wy){
        mean_variance(x,wx,sumw_x,mean_x,M2_x);
        M2_xy+=wy*wy*(x-mean_x)*(y-mean_y);
        sumw_xy+=wx*wy;
        mean_variance(y,wy,sumw_y,mean_y,M2_y);
        ++n;
      }
      T operator()() const {
        return covar()/sqrt(var_x()*var_y());
      }
      inline T covar() const {
        return variance(M2_xy,sumw_xy);
      }
      inline T var_x() const {
        return variance(M2_x ,sumw_x );
      }
      inline T var_y() const {
        return variance(M2_y ,sumw_y );
      }
      inline void clear(){
        M2_x   =0;
        M2_y   =0;
        M2_xy  =0;
        n      =0;
        mean_x =0;
        sumw_x =0;
        mean_y =0;
        sumw_y =0;
        mean_xy=0;
        sumw_xy=0;
      }
  };

  template<typename T>
  class cc_fisher{
    public:
      size_t n=0;
      T var_x=0;
      T var_y=0;
      T covar=0;
      cc_fisher& push(const T& x, const T&y){
        ++n;
        // TODO: numerically stable version a la mean_variance
        // -> calculate mean variance instead of Sum of squares
        var_x+=x*x;
        var_y+=y*y;
        covar+=x*y;
      }
      T operator()(){
        return covar/sqrt(var_x*var_y);
      }
      inline void clear(){
        var_x=0;
        var_y=0;
        covar=0;
      }
  };

  template <typename T>
  const void transpose(T m[], size_t h, size_t w)
  {
    for (size_t start = 0; start <= w * h - 1; ++start)
    {
      size_t next = start;
      size_t i = 0;
      do
      {
        ++i;
        next = (next % h) * w + next / h;
      } while (next > start);

      if (next >= start && i != 1)
      {
        const T tmp = m[start];
        next = start;
        do
        {
          i = (next % h) * w + next / h;
          m[next] = (i == start) ? tmp : m[i];
          next = i;
        } while (next > start);
      }
    }
  }

  template<class target>
  const inline double golden_section_search(double  a, double  va,
                                          double  b, double  vb,
                                          target function){
    //cerr << "entering golden section search" << endl;
    //cerr << a << " " << b << endl;
    //cerr << va << " " << vb << endl;
    const double gr = (sqrt(5.0) + 1.0) * 0.5;
    const double epsilon = 1e-6;
    while (true){
      const double c = b-(b-a)/gr;
      const double d = a+(b-a)/gr;
      //cerr << a << " " << c << " " << d << " " << b << endl;
      if (d-c<epsilon) break;
      const double vc=function(c);
      const double vd=function(d);
      if (vc<vd){
        b=d;
        vb=vd;
      }else{
        a=c;
        va=vc;
      }
      //cerr << a << " " << b << endl;
      //cerr << va << " " << vb << endl;
    }
    return a+(b-a)*0.5;
  }

  template<class target>
  const inline double linesearch(double c1,double v1,target function){
    double c0=0.0;
    double c2=2.0*c1;
    double v0,v2;
    v0=function(c0);
    v1=function(c1);
    v2=function(c2);
    for(double i=c1;;(i*=2.0,c0-=i,c2+=i)){
      v0=function(c0);
      v2=function(c2);
      if ((v0<v1)||(v2<v1)){
        if (v0<v1){
          v1=v0;
          c1=c0;
        }
        if (v2<v1){
          v1=v2;
          c1=c2;
        }
        continue;
      }else{
        break;
      }
    }
    //cerr << "first linesearch loop passed" << endl;
    return golden_section_search(c0,v0,c2,v2,function);
  }

  /*
  template<size_t N>
  class boolarray{
    uint_fast8_t[] data;
    boolarray(){
      data = new uint_fast8_t[N]();
    }
    boolarray(bool v=0){
      if (v){
        data = new uint_fast8_t[N];
        std::fill_n(data,N,~uint_fast8_t(0)); 
      }else{
        data =
          new uint_fast8_t[(N+digits<uint_fast8_t>-1)/digits<uint_fast8_t>]();
      }
    }
    const uint_fast8_t get(size_t x) const {
      return (uint_fast8_t(1))
        &(data[n/digits<uint_fast8_t>]>>(n%digits<uint_fast8_t>));
    }
    const uint_fast8_t set(size_t x,bool v){

    }
    uint_fast8_t operator[](size_t x){
      return (uint_fast8_t(1))
        &(data[n/digits<uint_fast8_t>]>>(n%digits<uint_fast8_t>));
    }
    ~boolarray(){
      delete[] data;
    }
  }

  template <typename T>
  T[] inline convolute2d(T[] a, const boolarray& m,
                                const size_t nx,const size_t ny,
                   const T[] kernel,
                                const size_t kx,const size_t ky){
    const size_t sx = kx/2;
    const size_t sy = ky/2;
    T[] b = new T[nx*ny]();
    for (size_t y=0;y<ny;++y){
      for (size_t x=0;x<nx;++x){
        for (size_t i=0;i<ky;++ky){
          for (size_t j=0;j<kx;++kx){
            const size_t k = y*nx+x;
            const size_t l = i*kx+j;
            const size_t o = (y+i-sy)*nx+(x+j-sx);
            b[k]+=kernel[l]*m[o]*a[o];
            n+=kernel[l]*m[o];
          }
        }
        b/=n;
      }
    }
  }
  
  template <typename T>
  T[] inline convolute_horizontal(T[] a, const boolarray& m,
                                         const size_t nx,const size_t ny,
                            const T[] kernel,
                                         const size_t kx){
    const size_t sx = kx/2;
    T[] b = new T[nx*ny]();
    for (size_t y=0;y<ny;++y){
      for (size_t x=0;x<nx;++x){
        for (size_t j=0;j<kx;++kx){
          const size_t k = y*nx+x;
          const size_t o = y*nx+x+j-sx;
          b[k]+=kernel[j]*m[o]*a[o];
          n+=kernel[j]*m[o];
        }
        b/=n;
      }
    }
  }

  template <typename T>
  T[] inline convolute_vertical(T[] a, const boolarray& m,
                                       const size_t nx,const size_t ny,
                          const T[] kernel,
                                       const size_t ky){
    const size_t sy = ky/2;
    T[] b = new T[nx*ny]();
    for (size_t y=0;y<ny;++y){
      for (size_t x=0;x<nx;++x){
        for (size_t i=0;i<ky;++ky){
            const size_t k = y*nx+x;
            const size_t o = (y+i-sy)*nx+x;
            b[k]+=kernel[i]*m[o]*a[o];
            n+=kernel[i]*m[o];
        }
        b/=n;
      }
    }
  }
  */
  template<typename T>
  typename std::enable_if<std::is_integral<T>::value,T>::type
  constexpr factorial(const T n){
    return n<2?1:n*factorial(n-1);
  }

  uint64_t power2_fixpoint(const array<uint64_t,2>){
    
  }

  array<uint64_t,2> log2_fixpoint_slow(const uint64_t n){
    uint64_t s = 0x6A09E667F3BCC909; // (sqrt(2)-1)*2**64
    uint64_t c = 0xB504F333F9DE6485; // 2**64/sqrt(2)
    uint64_t l = log2(n);
    uint64_t r = n<<(64-l);
    uint64_t m = 0u;
    uint64_t o = (1ull<<32)-1ull;
    for (size_t i=0;(i<64)&&r;++i){
      //cout << r << " " << s << " " << c << endl;
      if(r>=s){
        m|=uint64_t(1)<<(63-i);
        r=c+(c>>32)*(r>>32)+(((o&c)*(r>>32))>>32)+(((o&r)*(c>>32))>>32);
      }
      r=(r<<1)+(r>>32)*(r>>32)+(((o&r)*(r>>32))>>31);
    }
    return {{ l, m }};
  }

  array<uint64_t,2> log2_fixpoint_taylor(const uint64_t n){
    const uint64_t c = 0x71547652b82fe177;
    //const uint64_t c =   0x71547655ff2fe177;
    const uint64_t l = log2(n);
    const uint64_t o = 1ull<<31;
    const uint64_t r1 = n<<(64-l);
    const uint64_t r2 = ((r1+o)>>32)*((r1+o)>>32);
    const uint64_t r3 = ((r1+o)>>32)*((r2+o)>>32);
    const uint64_t r4 = ((r2+o)>>32)*((r2+o)>>32);
    const uint64_t r5 = ((r2+o)>>32)*((r3+o)>>32);
    const uint64_t r6 = ((r3+o)>>32)*((r3+o)>>32);
    const uint64_t r7 = ((r3+o)>>32)*((r4+o)>>32);
    const uint64_t r8 = ((r4+o)>>32)*((r4+o)>>32);
    /*const uint64_t r8 = (r4>>32)*(r4>>32);
    const uint64_t r9 = (r4>>32)*(r5>>32);
    const uint64_t ra = (r5>>32)*(r5>>32);
   const uint64_t rb = (r5>>32)*(r6>>32);
    const uint64_t rc = (r6>>32)*(r6>>32);
    const uint64_t rd = (r6>>32)*(r7>>32);
    const uint64_t re = (r7>>32)*(r7>>32);
    const uint64_t rf = (r7>>32)*(r8>>32);
    //uint64_t m = r;*/
    uint64_t m = r1;//+(c>>32)*(r>>32);
    //cout << r1 << " " << m << endl;
    m-=r2/2;
    //cout << r2 << " " << m << endl;
    m+=r3/3;
    //cout << r3 << " " << m << endl;
    m-=r4/4;
    //cout << r4 << " " << m << endl;
    m+=r5/5;
    //cout << r5 << " " << m << endl;
    m-=r6/6;
    //cout << r6 << " " << m << endl;
    m+=r7/7;
    //cout << r7 << " " << m << endl;
    /*
    m-=(r8+4)/8;
    m+=(r9+4)/9;
    cout << r8 << " " << m << endl;
    m-=(ra+5)/10;
    cout << r9 << " " << m << endl;
    m+=(rb+5)/11;
    cout << ra << " " << m << endl;
    m-=(rc+6)/12;
    cout << rb << " " << m << endl;
    m+=(rd+6)/13;
    cout << rc << " " << m << endl;
    m-=(re+7)/14;
    cout << rd << " " << m << endl;
    m+=(rf+7)/15;*/
    m+=(m>>32)*(c>>32);
    return {{ l, m }};
  }

/*
  array<uint64_t,2> log2_fixpoint(const uint64_t n){
    array<uint64_t,256> s = {{
        0x0000000000000000,
        0xFF4ECB59511EC8A5,
        0xFE9E115C7B8F884C,
        0xFDEDD1B496A89F35,
        0xFD3E0C0CF486C175,
        0xFC8EC01121E447BB,
        0xFBDFED6CE5F09C49,
        0xFB3193CC4227C3F4,
        0xFA83B2DB722A033A,
        0xF9D64A46EB939F35,
        0xF92959BB5DD4BA74,
        0xF87CE0E5B2094D9C,
        0xF7D0DF730AD13BB9,
        0xF7255510C4288239,
        0xF67A416C733F846E,
        0xF5CFA433E6537290,
        0xF5257D152486CC2C,
        0xF47BCBBE6DB9FDDF,
        0xF3D28FDE3A641A5B,
        0xF329C9233B6BAE9C,
        0xF281773C59FFB13A,
        0xF1D999D8B7708CC1,
        0xF13230A7AD094509,
        0xF08B3B58CBE8B76A,
        0xEFE4B99BDCDAF5CB,
        0xEF3EAB20E032BC6B,
        0xEE990F980DA3025B,
        0xEDF3E6B1D418A491,
        0xED4F301ED9942B84,
        0xECAAEB8FFB03AB41,
        0xEC0718B64C1CBDDC,
        0xEB63B74317369840,
        0xEAC0C6E7DD24392F,
        0xEA1E4756550EB27B,
        0xE97C38406C4F8C57,
        0xE8DA9958464B42AB,
        0xE8396A503C4BDC68,
        0xE798AADADD5B9CBF,
        0xE6F85AAAEE1FCE22,
        0xE658797368B3A717,
        0xE5B906E77C8348A8,
        0xE51A02BA8E26D681,
        0xE47B6CA0373DA88D,
        0xE3DD444C46499619,
        0xE33F8972BE8A5A51,
        0xE2A23BC7D7D91226,
        0xE2055AFFFE83D369,
        0xE168E6CFD3295D23,
        0xE0CCDEEC2A94E111,
        0xE031430A0D99E627,
        0xDF9612DEB8F04420,
        0xDEFB4E1F9D1037F2,
        0xDE60F4825E0E9124,
        0xDDC705BCD378F7F0,
        0xDD2D818508324C20,
        0xDC9467913A4F1C92,
        0xDBFBB797DAF23755,
        0xDB63714F8E295255,
        0xDACB946F2AC9CC72,
        0xDA3420ADBA4D8704,
        0xD99D15C278AFD7B6,
        0xD9067364D44A929C,
        0xD870394C6DB32C84,
        0xD7DA67311797F56A,
        0xD744FCCAD69D6AF4,
        0xD6AFF9D1E13BA2FE,
        0xD61B5DFE9F9BCE07,
        0xD5872909AB75D18A,
        0xD4F35AABCFEDFA1F,
        0xD45FF29E0972C561,
        0xD3CCF099859AC379,
        0xD33A5457A3029054,
        0xD2A81D91F12AE45A,
        0xD2164C023056BCAB,
        0xD184DF6251699AC6,
        0xD0F3D76C75C5DB8D,
        0xD06333DAEF2B2595,
        0xCFD2F4683F94EEB5,
        0xCF4318CF191918C1,
        0xCEB3A0CA5DC6A55D,
        0xCE248C151F8480E4,
        0xCD95DA6A9FF06445,
        0xCD078B86503DCDD2,
        0xCC799F23D11510E5,
        0xCBEC14FEF2727C5D,
        0xCB5EECD3B38597C9,
        0xCAD2265E4290774E,
        0xCA45C15AFCC72624,
        0xC9B9BD866E2F27A3,
        0xC92E1A9D517F0ECC,
        0xC8A2D85C8FFE2C45,
        0xC817F681416452B2,
        0xC78D74C8ABB9B15D,
        0xC70352F04336C51E,
        0xC67990B5AA245F79,
        0xC5F02DD6B0BBC3D9,
        0xC5672A115506DADD,
        0xC4DE8523C2C07BAA,
        0xC4563ECC5334CB33,
        0xC3CE56C98D21B15D,
        0xC346CCDA24976407,
        0xC2BFA0BCFAD907C9,
        0xC238D2311E3D6673,
        0xC1B260F5CA0FBB33,
        0xC12C4CCA66709456,
        0xC0A6956E8836CA8D,
        0xC0213AA1F0D08DB0,
        0xBF9C3C248E2486F8,
        0xBF1799B67A731083,
        0xBE935317FC378238,
        0xBE0F6809860993E2,
        0xBD8BD84BB67ED483,
        0xBD08A39F580C36BF,
        0xBC85C9C560E7B269,
        0xBC034A7EF2E9FB0D,
        0xBB81258D5B704B6F,
        0xBAFF5AB2133E45FB,
        0xBA7DE9AEBE5FEA09,
        0xB9FCD2452C0B9DEB,
        0xB97C143756844DBF,
        0xB8FBAF4762FB9EE9,
        0xB87BA337A1743834,
        0xB7FBEFCA8CA41E7C,
        0xB77C94C2C9D725E9,
        0xB6FD91E328D17791,
        0xB67EE6EEA3B22B8F,
        0xB60093A85ED5F76C,
        0xB58297D3A8B9F0D2,
        0xB504F333F9DE6484,
        0xB487A58CF4A9C180,
        0xB40AAEA2654B9841,
        0xB38E0E38419FAE18,
        0xB311C412A9112489,
        0xB295CFF5E47DB4A4,
        0xB21A31A66618FE3B,
        0xB19EE8E8C94FEB09,
        0xB123F581D2AC2590,
        0xB0A957366FB7A3C9,
        0xB02F0DCBB6E04584,
        0xAFB51906E75B8661,
        0xAF3B78AD690A4375,
        0xAEC22C84CC5C9465,
        0xAE493452CA35B80E,
        0xADD08FDD43D01491,
        0xAD583EEA42A14AC6,
        0xACE0413FF83E5D04,
        0xAC6896A4BE3FE929,
        0xABF13EDF162675E9,
        0xAB7A39B5A93ED337,
        0xAB0386EF48868DE1,
        0xAA8D2652EC907629,
        0xAA1717A7B5693979,
        0xA9A15AB4EA7C0EF8,
        0xA92BEF41FA77771B,
        0xA8B6D5167B320E09,
        0xA8420BFA298F70D1,
        0xA7CD93B4E965356A,
        0xA7596C0EC55FF55B,
        0xA6E594CFEEE86B1E,
        0xA6720DC0BE08A20C,
        0xA5FED6A9B15138EA,
        0xA58BEF536DBEB6EE,
        0xA5195786BE9EF339,
        0xA4A70F0C95768EC5,
        0xA43515AE09E6809E,
        0xA3C36B345991B47C,
        0xA3520F68E802BB93,
        0xA2E102153E918F9E,
        0xA27043030C496819,
        0xA1FFD1FC25CEA188,
        0xA18FAECA8544B6E4,
        0xA11FD9384A344CF7,
        0xA0B0510FB9714FC2,
        0xA041161B3D0121BE,
        0x9FD228256400DD06,
        0x9F6386F8E28BA651,
        0x9EF5326091A111AE,
        0x9E872A276F0B98FF,
        0x9E196E189D472420,
        0x9DABFDFF6367A2AA,
        0x9D3ED9A72CFFB751,
        0x9CD200DB8A0774CB,
        0x9C6573682EC32C2D,
        0x9BF93118F3AA4CC1,
        0x9B8D39B9D54E5539,
        0x9B218D16F441D63D,
        0x9AB62AFC94FF864A,
        0x9A4B13371FD166CA,
        0x99E0459320B7FA65,
        0x9975C1DD47518C77,
        0x990B87E266C189AA,
        0x98A1976F7597E996,
        0x9837F0518DB8A96F,
        0x97CE9255EC4357AB,
        0x97657D49F17AB08E,
        0x96FCB0FB20AC4BA3,
        0x96942D3720185A00,
        0x962BF1CBB8D97560,
        0x95C3FE86D6CC7FEF,
        0x955C5336887894D5,
        0x94F4EFA8FEF70961,
        0x948DD3AC8DDB7ED3,
        0x9426FF0FAB1C04B6,
        0x93C071A0EEF94BC1,
        0x935A2B2F13E6E92C,
        0x92F42B88F673AA7C,
        0x928E727D9531F9AC,
        0x9228FFDC10A051AD,
        0x91C3D373AB11C336,
        0x915EED13C89689D3,
        0x90FA4C8BEEE4B12B,
        0x9095F1ABC540CA6B,
        0x9031DC431466B1DC,
        0x8FCE0C21C6726481,
        0x8F6A8117E6C8E5C4,
        0x8F073AF5A2013520,
        0x8EA4398B45CD53C0,
        0x8E417CA940E35A01,
        0x8DDF042022E69CD6,
        0x8D7CCFC09C50E2F8,
        0x8D1ADF5B7E5BA9E6,
        0x8CB932C1BAE97A95,
        0x8C57C9C4646F4DDE,
        0x8BF6A434ADDE0085,
        0x8B95C1E3EA8BD6E7,
        0x8B3522A38E1E1032,
        0x8AD4C6452C728924,
        0x8A74AC9A79896E47,
        0x8A14D575496EFD9A,
        0x89B540A7902557A4,
        0x8955EE03618E5FDD,
        0x88F6DD5AF155AC6B,
        0x88980E8092DA8527,
        0x88398146B919F1D4,
        0x87DB357FF698D792,
        0x877D2AFEFD4E256C,
        0x871F61969E8D1010,
        0x86C1D919CAEF5C88,
        0x8664915B923FBA04,
        0x86078A2F23642A9F,
        0x85AAC367CC487B15,
        0x854E3CD8F9C8C95D,
        0x84F1F656379C1A29,
        0x8495EFB3303EFD30,
        0x843A28C3ACDE4046,
        0x83DEA15B9541B132,
        0x8383594EEFB6EE37,
        0x83285071E0FC4547,
        0x82CD8698AC2BA1D7,
        0x8272FB97B2A5894C,
        0x8218AF4373FC25EC,
        0x81BEA1708DDE6056,
        0x8164D1F3BC030773,
        0x810B40A1D81406D4,
        0x80B1ED4FD999AB6C,
        0x8058D7D2D5E5F6B1
    }};
    array<uint64_t,256> c={{
        0x0000000000000000,
        0x01709C46D7AAC774,
        0x02DFCA16DDE10A2F,
        0x044D8C45EA5EC312,
        0x05B9E5A170B48A62,
        0x0724D8EEA143E199,
        0x088E68EA899A0976,
        0x09F6984A342D1310,
        0x0B5D69BAC77EC398,
        0x0CC2DFE1A4A8CA30,
        0x0E26FD5C8555AF7A,
        0x0F89C4C19929CFD0,
        0x10EB389FA29F9AB3,
        0x124B5B7E135A3C89,
        0x13AA2FDD27F1C2D8,
        0x1507B836033BB6D4,
        0x1663F6FAC913167C,
        0x17BEEE96B8A2813C,
        0x1918A16E46335AAE,
        0x1A7111DF348493EB,
        0x1BC84240ADABBA63,
        0x1D1E34E35B82DA4D,
        0x1E72EC117FA5B21C,
        0x1FC66A0F0B00A490,
        0x2118B119B4F3C72C,
        0x2269C369120C5BFE,
        0x23B9A32EAA56F6BD,
        0x250852960F4C7F1A,
        0x2655D3C4F15C343E,
        0x27A228DB3514C7E0,
        0x28ED53F307EE9A62,
        0x2A375720F4B91491,
        0x2B803473F7AD0F3F,
        0x2CC7EDF592262CF9,
        0x2E0E85A9DE04FE53,
        0x2F53FD8FA0BBBDE2,
        0x309857A05E0765FB,
        0x31DB95D06A56D76B,
        0x331DBA0EFCE1BE05,
        0x345EC6464170D594,
        0x359EBC5B69D927DF,
        0x36DD9E2EBF2BD2D7,
        0x381B6D9BB29BDC81,
        0x39582C78EE1B912D,
        0x3A93DC9864B2DF91,
        0x3BCE7FC762901DC3,
        0x3D0817CE9CD4998F,
        0x3E40A672411E4E8E,
        0x3F782D7204D01447,
        0x40AEAE8934198EEC,
        0x41E42B6EC0C025BC,
        0x4318A5D550AA3A93,
        0x444C1F6B4C2DD72C,
        0x457E99DAEC23FD65,
        0x46B016CA47C1C14A,
        0x47E097DB62384CAF,
        0x49101EAC381CE609,
        0x4A3EACD6CC9A0D8C,
        0x4B6C43F1366ABDBC,
        0x4C98E58DACA0D66B,
        0x4DC4933A9337B366,
        0x4EEF4E828773EA56,
        0x501918EC6C1125D6,
        0x5141F3FB753F0E53,
        0x5269E12F346E2BF9,
        0x5390E203A3EDA7EE,
        0x54B6F7F1325ACDF7,
        0x55DC246CCDE32AC4,
        0x570068E7EF5A1E7E,
        0x5823C6D0A522B65B,
        0x59463F919DEE9B94,
        0x5A67D4923352E1A3,
        0x5B8887367433795E,
        0x5CA858DF2F060A50,
        0x5DC74AE9FBECEF91,
        0x5EE55EB146AB115D,
        0x6002958C587150CA,
        0x611EF0CF6186371F,
        0x623A71CB82C89692,
        0x635519CED70DC6AC,
        0x646EEA247C5C22D2,
        0x6587E4149D026E33,
        0x66A008E4788CBCD2,
        0x67B759D66C977E23,
        0x68CDD829FD814275,
        0x69E3851BDEFBD14A,
        0x6AF861E5FC7D2386,
        0x6C0C6FBF8190D199,
        0x6D1FAFDCE20A8290,
        0x6E32236FE219E65C,
        0x6F43CBA79E40C2AD,
        0x7054A9B0932B9700,
        0x7164BEB4A56D59F9,
        0x72740BDB291ECF4D,
        0x73829248E961F325,
        0x749053202FC9F548,
        0x759D4F80CBA83BF8,
        0x76A98888193EE404,
        0x77B4FF5108D9313A,
        0x78BFB4F425CA6037,
        0x79C9AA879D534831,
        0x7AD2E11F456F394F,
        0x7BDB59CCA38881F4,
        0x7CE3159EF3150343,
        0x7DEA15A32C1B3B38,
        0x7EF05AE409A0288E,
        0x7FF5E66A0FFE6AE7,
        0x80FAB93B9326FF91,
        0x81FED45CBCCBF99C,
        0x830238CF927591F4,
        0x8404E793FB81EA92,
        0x8506E1A7C70FDE00,
        0x86082806B1D532C4,
        0x8708BBAA6BE0889D,
        0x88089D8A9E4753D8,
        0x8907CE9CF0C0396A,
        0x8A064FD50F2A1CF0,
        0x8B042224AF00303C,
        0x8C01467B94BB5275,
        0x8CFDBDC799210B8B,
        0x8DF988F4AE806F1D,
        0x8EF4A8ECE5DD30C3,
        0x8FEF1E9874093212,
        0x90E8EADDB6ACD18A,
        0x91E20EA1393E4040,
        0x92DA8AC5B9E822B4,
        0x93D2602C2E5FC02B,
        0x94C98FB3C8AB028F,
        0x95C01A39FBD6879F,
        0x96B6009A809C0325,
        0x97AB43AF59F930A1,
        0x989FE450D9B791CF,
        0x9993E355A4E53643,
        0x9A874192B83EC743,
        0x9B79FFDB6C8B1202,
        0x9C6C1F017AE84940,
        0x9D5D9FD5010B3666,
        0x9E4E832485709116,
        0x9F3EC9BCFB80B357,
        0xA02E7469C7A5DF5F,
        0xA11D83F4C3554B38,
        0xA20BF926410B256D,
        0xA2F9D4C51039C526,
        0xA3E71796812C371A,
        0xA4D3C25E68DC57F2,
        0xA5BFD5DF24BCABE2,
        0xA6AB52D99E762253,
        0xA7963A0D4F99F3BA,
        0xA8808C384547C6EF,
        0xA96A4A1723C84A77,
        0xAA5374652A1C6D8D,
        0xAB3C0BDC358163D8,
        0xAC241134C4E99E1C,
        0xAD0B8525FC6AE15E,
        0xADF26865A8A1A557,
        0xAED8BBA84209E24E,
        0xAFBE7FA0F04D75C6,
        0xB0A3B5018D8844CE,
        0xB1885C7AA9824203,
        0xB26C76BB8CDF7CC3,
        0xB35004723C465E69,
        0xB433064B7B7C39BC,
        0xB5157CF2D0785040,
        0xB5F76912866D705C,
        0xB6D8CB53B0CA4ECB,
        0xB7B9A45E2E30BD34,
        0xB899F4D8AB63DF28,
        0xB979BD68A62D7E5C,
        0xBA58FEB2703A9E37,
        0xBB37B95931EF6E6F,
        0xBC15EDFEED32BBDE,
        0xBCF39D448030FE26,
        0xBDD0C7C9A817204F,
        0xBEAD6E2D03C52202,
        0xBF89910C1678AE89,
        0xC06531034A6FC648,
        0xC1404EADF38396DE,
        0xC21AEAA651BB9DBE,
        0xC2F5058593D93084,
        0xC3CE9FE3D9DB85F3,
        0xC4A7BA58377C5A03,
        0xC5805578B6A5470B,
        0xC65871DA59DDED9B,
        0xC73010111EB3043F,
        0xC80730B0001667F2,
        0xC8DDD448F8B845A5,
        0xC9B3FB6D055974E6,
        0xCA89A6AC27171B3A,
        0xCB5ED69565AFAF7F,
        0xCC338BB6D1C1742A,
        0xCD07C69D87027EF4,
        0xCDDB87D5AE726421,
        0xCEAECFEA80859B33,
        0xCF819F66474AB28E,
        0xD053F6D260896731,
        0xD125D6B73FDBB557,
        0xD1F73F9C70C0F683,
        0xD2C8320898AB212F,
        0xD398AE8179063DEA,
        0xD468B58BF13A2595,
        0xD53847AC00A69BE6,
        0xD6076564C899D941,
        0xD6D60F388E41968A,
        0xD7A445A8BC96AD5B,
        0xD8720935E6435EBD,
        0xD93F5A5FC7845236,
        0xDA0C39A548045ECB,
        0xDAD8A7847CB32F34,
        0xDBA4A47AA996D25A,
        0xDC703104439848E9,
        0xDD3B4D9CF24B2077,
        0xDE05FABF91B02C8E,
        0xDED038E633F36DA8,
        0xDF9A088A232535E6,
        0xE0636A23E2EE9B16,
        0xE12C5E2B3241455B,
        0xE1F4E5170D02A99B,
        0xE2BCFF5DADB2BE89,
        0xE384AD748F0E3B05,
        0xE44BEFD06DAC6C39,
        0xE512C6E54998B1AF,
        0xE5D9332667E7AD76,
        0xE69F350654483612,
        0xE764CCF6E29017EB,
        0xE829FB693044B398,
        0xE8EEC0CDA61F8649,
        0xE9B31D93F98EA94B,
        0xEA77122B2E315596,
        0xEB3A9F01975077F1,
        0xEBFDC484D953623C,
        0xECC08321EB30A61E,
        0xED82DB4517DB253D,
        0xEE44CD59FFAB62F3,
        0xEF0659CB99C52334,
        0xEFC781043579625F,
        0xF088436D7BA4B14F,
        0xF148A170700A00FD,
        0xF2089B7572A9E8D7,
        0xF2C831E4411672B0,
        0xF3876523F7C37731,
        0xF446359B13539551,
        0xF504A3AF71E1CF71,
        0xF5C2AFC65447D86A,
        0xF6805A445F611AB8,
        0xF73DA38D9D4A83EB,
        0xF7FA8C057E9F1E32,
        0xF8B7140EDBB181D9,
        0xF9733C0BF5C22856,
        0xFA2F045E7832AA72,
        0xFAEA6D6779B5F2E8,
        0xFBA577877D7D6EBD,
        0xFC60231E74634472,
        0xFD1A708BBE119B14,
        0xFDD4602E2A26F9FD,
        0xFE8DF263F957CA15,
        0xFF47278ADE8D012E
    }};

  //  uint64_t s = 0x172B83C7D517ADCE; // (2**( 1/8)-1)*2**64
  //  uint64_t c = 0xEAC0C6E7DD24392F; //  2**(-1/8)   *2**64
    uint64_t l = log2(n);
    uint64_t r = n<<(64-l);
    uint64_t m = 0u;
    uint64_t o = (1ull<<32)-1ull;
    // look at next bit add c[bit] to m
    // multiply r by s[bit]
    //
    for (size_t i=0;(i<64)&&r;++i){
      //cout << r << " " << s << " " << c << endl;
      if(r>=s){
        m|=uint64_t(1)<<(63-i);
        r=c+(c>>32)*(r>>32)+(((o&c)*(r>>32))>>32)+(((o&r)*(c>>32))>>32);
      }
      r=(r<<1)+(r>>32)*(r>>32)+(((o&r)*(r>>32))>>31);
    }
  }
*/
  /*
  // triangle indices are hashed in a way to allow
  // fast traversal if i is fixed and j ++ or --
  // and max(i) is not known in general
  // i<j or else undefined.
  // {{0 1 2 3 4 ... }= i
  // 0|- - - - - ...
  // 1|0 - - - - ...
  // 2|1 2 - - - ...
  // 3|3 4 5 - - ...
  // 4|6 7 8 9 - ...
  // ...
  // } = j
  struct triangle_index{
    size_t i;
    size_t j;
    triangle_index(const size_t& i, const size_t& j):i(i),j(j){
      // if(i>j) std::swap(i,j);
    }
    inline friend bool operator==(triangle_index a, triangle_index b){
      return ((a.i==b.i)&&(a.j==b.j))||((a.i==b.j)&&(a.j==b.i));
    }
  };
  template<typename T>
  struct triangle_matrix{
    size_t N;
    vector<T> data;
    triangle_matrix(const size_t& N):N(N){
      data.resize(N*(N-1)/2);
      //data.resize(N*N);
    }
    inline T& at(const size_t& i, const size_t& j){
      return data[j*(j-1)/2+i];
      //return data[i*N+j];
    }
  };
  */


  /*
  constexpr uint16_t morton_encode(const uint8_t& x, const uint8_t& y) {
    return (
            (x*0x0101010101010101ULL&0x8040201008040201ULL)
            *0x0102040810204081ULL>>49
           )&0x5555|
           (
            (y*0x0101010101010101ULL&0x8040201008040201ULL)
            *0x0102040810204081ULL>>48
           )&0xAAAA;
  }
  
  const uint32_t morton_encode(uint32_t x, uint32_t y, uint32_t z) {
    x = (x | (x << 16)) & 0x030000FF;
    x = (x | (x <<  8)) & 0x0300F00F;
    x = (x | (x <<  4)) & 0x030C30C3;
    x = (x | (x <<  2)) & 0x09249249;

    y = (y | (y << 16)) & 0x030000FF;
    y = (y | (y <<  8)) & 0x0300F00F;
    y = (y | (y <<  4)) & 0x030C30C3;
    y = (y | (y <<  2)) & 0x09249249;

    z = (z | (z << 16)) & 0x030000FF;
    z = (z | (z <<  8)) & 0x0300F00F;
    z = (z | (z <<  4)) & 0x030C30C3;
    z = (z | (z <<  2)) & 0x09249249;

    return x^(x<<1)^(z<<2);
  }
  
  const uint32_t morton_encode(uint32_t x, uint32_t y){
    x &= 0x0000ffff;                 
    x = (x ^ (x <<  8)) & 0x00ff00ff; 
    x = (x ^ (x <<  4)) & 0x0f0f0f0f; 
    x = (x ^ (x <<  2)) & 0x33333333; 
    x = (x ^ (x <<  1)) & 0x55555555;

    y &= 0x0000ffff;                  
    y = (y ^ (y <<  8)) & 0x00ff00ff;   
    y = (y ^ (y <<  4)) & 0x0f0f0f0f;     
    y = (y ^ (y <<  2)) & 0x33333333;       
    y = (y ^ (y <<  1)) & 0x55555555;

    return (x<<1)^y;
  }
  */
  /* TODO general morton encode that is fast
  template <typename T,size_t masksize,size_t blocksize=1>
  typename std::enable_if<std::is_integral<T>::value,T>::type
  constexpr morton_mask(){
    T x=0;
    for (size_t i=0; i<digits<T>(); i=i+masksize+blocksize){
      x^=((T(1)<<blocksize)-1)<<i;
    }
    return x;
    
    if ((d+masksize+blocksize)>(digits<T>())) return (T(1)<<blocksize-T(1));
    return ((T(1)<<blocksize-T(1))
           ^(morton_mask<T,masksize,blocksize>(d+masksize+blocksize)
           <<(blocksize+masksize)));
    
    return (blocksize>digits<T>())?(~T(0)):
           (d+masksize+blocksize>=(digits<T>()))?((T(1)<<blocksize)-T(1)):
             (((T(1)<<blocksize)-T(1))
              ^(morton_mask<T>(masksize,blocksize,d+masksize+blocksize)
               <<(blocksize+masksize)));
  }

  template<typename T,size_t n>
  typename std::enable_if<std::is_integral<T>::value,T>::type
  constexpr morton_shift(T x){
    x&=(T(1)<<(digits<T>()/(n+1)));
    for(size_t i=n*digits<T>()/4; i>=n; i=i/2){
      x=(x^(x<<i))&morton_mask<T,i*n,i>();
    }
    return x;
  }

  template<size_t N, typename T>
  typename std::enable_if<std::is_integral<T>::value,T>::type
  const morton_encode(const T x){
    return morton_shift<T,N-1>(x);
  }

  template<size_t N, typename T, typename... Ts>
  typename std::enable_if<std::is_integral<T>::value,T>::type
  const morton_encode(const T x, const Ts ... xs){
    return (morton_shift<T,N-1>(T(x))<<sizeof...(xs))
          ^(morton_encode<N,Ts...>(xs...));
  }
  */
}

namespace std {
  template <>
  struct hash<array<int32_t,2>>
  {
    //array<slow,fast> dimensinos
    size_t operator()(const array<uint32_t,2>& x) const {
      return (size_t(x[0])<<32)^(size_t(x[1]));
    }
  };
  template<>
  struct hash<array<uint32_t,2>> {
    size_t operator()(const array<uint32_t,2>& x) const {
      return (size_t(x[0])<<32)^(size_t(x[1]));
    }
  };
  /*
  template <>
  struct hash<wmath::triangle_index>
  {
    size_t operator()(const wmath::triangle_index& x) const{
      return x.i*(x.i-1)/2+x.j;
    }
  };
  */

}
#endif // WMATH_H

