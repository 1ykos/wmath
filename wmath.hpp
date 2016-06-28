#ifndef WMATH_H
#define WMATH_H

#include <array>
#include <algorithm>
#include <cstdint>
#include <functional>
#include <iostream>
#include <iterator>
#include <limits>
#include <iterator>
//#include <bitset>

using std::array;
using std::abs;
using std::vector;
using std::enable_if;
using std::min_element;
using std::remove_if;
using std::is_integral;
using std::is_unsigned;
using std::is_floating_point;
using std::min_element;
using std::swap;
//using std::cout;
//using std::endl;
using std::accumulate;
using std::unique;
using std::iterator;
using std::iterator_traits;
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

  template <typename T>
  constexpr size_t digits(){
    return size_t(std::numeric_limits<T>::digits);
  }

  struct popcount{
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
    }
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
  };
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
    const T mask = (std::numeric_limits<T>::digits-1);
    const T c = i&mask;
    return (n<<c)|(n>>((-c)&mask ));
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
  constexpr log2(const T& x,const T& lower=0,const T& upper=digits<T>()){
    return (upper-lower==T(1))?lower:(x&(T(0)-T(1)<<((upper+lower)/2))?
           log2(x,(upper+lower)/2,upper):
           log2(x,lower,(upper+lower)/2));
  }

  // levensthein coding in the range [0,2^51]
  const inline uint64_t universal_encode_uint64(uint64_t n){
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

  const inline uint64_t universal_encode_uint64(uint64_t n,uint64_t& l){
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

/*
  const inline uint64_t elias_delta_encode_uint64(const uint64_t n,uint64_t& len,uint64_t& lol){
    len = log2(n+1);
    lol = log2(len+1);
    return ((n+1)<<(64-2*lol-len))|(lol<<(64-lol));
  }
*/
  
  const inline uint64_t universal_decode_uint64(uint64_t i){
    uint64_t log2star = 0;
    while (i&(uint64_t(1)<<63-log2star)) ++log2star; // could be done with count leading zeroes
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

  const inline uint64_t universal_decode_uint64(uint64_t i,uint64_t& l){
    uint64_t log2star = 0;
    while (i&(uint64_t(1)<<63-log2star)) ++log2star; // could be done with count leading zeroes
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
    return (n<<1)^(n>>std::numeric_limits<T>::digits-1);
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
  class mean_variance_calculator{
    public:
      T M2=0;
      T sumw=0;
    public:
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

  /*
TODO: numerically stable implementation
calculate all moments with compensated summation

  template<typename T, size_t N>
  inline const void mean_moments(
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

