#ifndef WMATH_MATH_H
#define WMATH_MATH_H
#include "wmath_forward.hpp"
#include "wmath_bits.hpp"

namespace wmath{
  template <typename T> int signum(T val) {
    return (T(0) < val) - (val < T(0));
  }

  template <typename T0,typename T1>
  typename std::enable_if<std::is_unsigned<T1>::value,T0>::type
  constexpr pow(const T0& x,const T1& n){
    T1 m = n;
    T0 y = 1;
    T0 z = x;
    while (m) {
      if (m&1u) y*=z;
      z*=z;
      m>>=1;
    }
    return y;
  }
  
  template<typename T>
  double inline pow(const T& x,const double& e){
    return std::pow(x,e);
  }

  template<typename T0,typename T1,typename T2>
  constexpr T0 clip(const T0& n,const T1& l,const T2& h){
    return n<l?l:n>h?h:n;
  }
  
  template<typename T>
  typename std::enable_if<std::is_floating_point<T>::value,T>::type
  constexpr no_nan(const T& v){
    return isnan(v)?0:v;
  }

  template <typename T0,typename T1>
  typename std::enable_if<std::is_signed<T1>::value,double>::type
  constexpr pow(const T0& x,const T1& n){
    typedef typename make_unsigned<T1>::type U;
    return n<0?1.0/pow(x,U(-n)):pow(x,U(n));
  }
  
  template <typename T>
  typename std::enable_if<std::is_integral<T>::value,T>::type
  constexpr roundup_2_3(const T& n){
    if (n<0) return -roundup_2_3(-n);
    if (n==0) return 1;
    if (n<4) return n;
    if (n%2==0) return 2*roundup_2_3(n/2);
    if (n%3==0) return 3*roundup_2_3(n/3);
    return roundup_2_3(n+1);
  }
  
  template <typename T>
  typename std::enable_if<std::is_integral<T>::value,T>::type
  constexpr roundup_2(const T& n){
    if (n<0) return -roundup_2(-n);
    if (n==0) return 1;
    if (n<3) return n;
    if (n%2==0) return 2*roundup_2(n/2);
    return roundup_2(n+1);
  }

  template <typename T>
  typename std::enable_if<std::is_unsigned<T>::value,T>::type
  constexpr circadd(const T& a,const T& b){
    return a+b+(T(a+b)<a);
  }
  
  template <typename T>
  typename std::enable_if<std::is_unsigned<T>::value,T>::type
  constexpr circdff(const T& a,const T& b){
    return a-b-(a-b>a);
  }
  
  template <typename T>
  typename std::enable_if<std::is_unsigned<T>::value,T>::type
  constexpr circadd(const T& a,const T& b,const T& m){
    return (a+b)%m+((a+b)%m<a);
  }

  template<typename T>
  //typename std::enable_if<std::is_unsigned<T>::value,T>::type
  T
  fixponit_integer_inverse(const T& d);

  template<typename T>
  //typename std::enable_if<std::is_unsigned<T>::value,T>::type
  T
  constexpr fixpoint_integer_inverse(const T& d) {
    uint8_t lut[256] = {
     255u,254u,253u,252u,251u,250u,249u,248u,247u,246u,245u,244u,243u,242u,241u,
240u,240u,239u,238u,237u,236u,235u,234u,234u,233u,232u,231u,230u,229u,229u,228u,
227u,226u,225u,225u,224u,223u,222u,222u,221u,220u,219u,219u,218u,217u,217u,216u,
215u,214u,214u,213u,212u,212u,211u,210u,210u,209u,208u,208u,207u,206u,206u,205u,
204u,204u,203u,202u,202u,201u,201u,200u,199u,199u,198u,197u,197u,196u,196u,195u,
195u,194u,193u,193u,192u,192u,191u,191u,190u,189u,189u,188u,188u,187u,187u,186u,
186u,185u,185u,184u,184u,183u,183u,182u,182u,181u,181u,180u,180u,179u,179u,178u,
178u,177u,177u,176u,176u,175u,175u,174u,174u,173u,173u,172u,172u,172u,171u,171u,
170u,170u,169u,169u,168u,168u,168u,167u,167u,166u,166u,165u,165u,165u,164u,164u,
163u,163u,163u,162u,162u,161u,161u,161u,160u,160u,159u,159u,159u,158u,158u,157u,
157u,157u,156u,156u,156u,155u,155u,154u,154u,154u,153u,153u,153u,152u,152u,152u,
151u,151u,151u,150u,150u,149u,149u,149u,148u,148u,148u,147u,147u,147u,146u,146u,
146u,145u,145u,145u,144u,144u,144u,144u,143u,143u,143u,142u,142u,142u,141u,141u,
141u,140u,140u,140u,140u,139u,139u,139u,138u,138u,138u,137u,137u,137u,137u,136u,
136u,136u,135u,135u,135u,135u,134u,134u,134u,134u,133u,133u,133u,132u,132u,132u,
132u,131u,131u,131u,131u,130u,130u,130u,130u,129u,129u,129u,129u,128u,128u,128u,
127u
    };
    const auto l = log2(d);
    //T x = T(1)<<(digits(d)-1-l);
    //cout << std::hex << ((d-(T(1)<<l))>>(l-8)) << endl;
    //cout << (~T(0))/d << endl;
    T x;
    if (l<8) {
      x = T(1)<<(digits(d)-1-l);
    } else {
      if (digits(d)>(l+8)) x = T(lut[(d>>(l-8))-256])<<(digits(d)-l-8);
      else x = T(lut[(d>>(l-8))-256])>>(l+8-digits(d));
      //cout << std::dec << ((d>>(l-8))-256) << endl;
    }
    if (x==0) x=1;
    //T x = (~T(0))>>(log2(d)+1);
    // x = n/d  <->  d*x = n  <->  d = n/x
    // f(x) = (1<<digits(d))/x - d 
    // f'(x) = -(1<<digits(d))/(x*x)
    // x(i+1) = x - f(x)/f'(x) = x + ( x*((1<<digits(d)) - x*d) ) >> 64
    //x+=(x*(-x*d)>>digits(d);
    //cout << endl;
    while(true) {
      //cout << std::hex << "x = " << x << " -x*d = " << T(0)-x*d << endl;
      const auto lm = long_mul(x+1,T(0)-x*d);
      const T i = get<0>(lm);
      //cout << " i = " << i << endl;
      if (i) x+=i;
      else return x;
    }
    return x;
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
  constexpr floor_square_root(const T& n){
    if (n==0) return 0;
    T l = log2(n)/T(2);
    T t0 = ((n>>l)+(T(1)<<l))/T(2);
    T t1(0),t2(0); 
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

  /* this unfortunately does not work as std::linear_congruential_engine
   * has compile time parameters a c m ... TODO: runtime version of
   * linear_congruential_engine
  template <typename T,typename RNG>
  typename std::enable_if<
    std::is_integral<T>::value,
    linear_congruential_engine>::type
  >
  const inline get_linear_congruential_engine(const T& m,RNG rng){
    if (m<3) return linear_congruential_engine<T,1,1,m>;
    if (is_prime(m) return linear_congruential_engine<T,
    auto factors = prime_factorization(m);
    std::uniform_int_distribution<T> test_c(1,n-1);
    std::geometric_distribution<T> gdistr(1.0/3.0);
    T c;
    while(true){
      c = test_c(rng);
      for (auto it=factors.begin();it!=factors.end();++it) if(c%(*it)) continue;
      break; // c and m share no divisors
    }
    T a = 1;
    for (auto it=factors.begin();it!=factors.end();++it){
      a*=*it;
    }
    if (m%4==0) a*=2;
    T b = a;
    while(b<m-1){
      a = b;
      b*=2+gdistr(rng);
    }
    a+=1;
    return linear_congruential_engine<T,a,c,m>;
  }
  */
  

/* numerically stable and incremental mean and variance
 * start with T sumw=0.0, mean=0.0, M2=0.0;
 * or with T sumw=w_1, mean=x_1, M2=0;
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
    if (w<numeric_limits<T>::min()) return;
    const T temp = w+sumw;
    const T delta = x-mean;
    const T R = delta*w/temp;
    mean += R;
    M2   += sumw*delta*R; // sumw*(x-mean)*(x-mean)*w/(sumw+w)
    sumw  = temp;
  }

/* numerically stable and incremental mean
 * start with T sumw=0.0, mean=0.0
 * or with T sumw=w_1, mean=x_1;
 * and then call
 * mean_variance(x,w,sumw,mean)
 * for each pair of x and w */  
  template<typename T>
  typename enable_if<is_floating_point<T>::value>::type
  const inline mean_variance(
      const T &x,
      const T &w,
      T &sumw,
      T &mean
      ){
    T M2 = 0; // mean_variance does not depend on previous values of M2
    mean_variance(x,w,sumw,mean,M2);
    //mean += (x-mean)*w/(sumw+w);
    //sumw += w;
  }
  
  template<typename T>
  typename enable_if<is_floating_point<T>::value,T>::type
  const inline mean_variance_sumw_dx(
      const T &x,
      const T &w,
      const T &sumw,
      const T &mean,
      const T &M2
      ){
    return 0;
  }
  
  template<typename T>
  typename enable_if<is_floating_point<T>::value,T>::type
  const inline mean_variance_sumw_dw(
      const T &x,
      const T &w,
      const T &sumw,
      const T &mean,
      const T &M2
      ){
    return 1;
  }
  
  template<typename T>
  typename enable_if<is_floating_point<T>::value,T>::type
  const inline mean_variance_mean_dx(
      const T &x,
      const T &w,
      const T &sumw,
      const T &mean,
      const T &M2
      ){
    return w/sumw;
  }

  template<typename T>
  typename enable_if<is_floating_point<T>::value,T>::type
  const inline mean_variance_mean_dw(
      const T &x,
      const T &w,
      const T &sumw,
      const T &mean,
      const T &M2
      ){
    return (x-mean)/sumw;
  }

  template<typename T>
  typename enable_if<is_floating_point<T>::value,T>::type
  const inline mean_variance_M2_dx_ini( // initial
      const T &x,
      const T &w,
      const T &sumw,
      const T &mean,
      const T &M2
      ){
    return 2*w*(1-w/sumw)*(x-mean);
  }
  
  template<typename T>
  typename enable_if<is_floating_point<T>::value,T>::type
  const inline mean_variance_M2_dx_acc( // needs to be accumulated
      const T &x,
      const T &w,
      const T &w0,
      const T &sumw,
      const T &mean,
      const T &M2
      ){
    return -2*w0*w*(x-mean)/sumw;
  }

  template<typename T>
  typename enable_if<is_floating_point<T>::value,T>::type
  const inline mean_variance_M2_dw_ini(
      const T &x,
      const T &w,
      const T &sumw,
      const T &mean,
      const T &M2
      ){
    return (x-mean)*(x-mean);
  }


//       (x0-mean)^2
// +2*w3*(mean-x0)*(x3-mean)/sumw
// +2*w2*(mean-x0)*(x2-mean)/sumw
// +2*w1*(mean-x0)*(x1-mean)/sumw
// +2*w0*(mean-x0)*(x0-mean)/sumw
// + ...
  template<typename T>
  typename enable_if<is_floating_point<T>::value,T>::type
  const inline mean_variance_M2_dw_acc(
      const T &x,
      const T &w,
      const T &x0,
      const T &sumw,
      const T &mean,
      const T &M2
      ){
    return 2*w/sumw*(mean-x0)*(x-mean);
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
  void inline transpose(T m[], size_t h, size_t w)
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
  inline double golden_section_search(double  a, double  va,
                                      double  b, double  vb,
                                      const target& function,
                                      const double epsilon = 1e-8){
    //cerr << "entering golden section search" << endl;
    //cerr << a << " " << b << endl;
    //cerr << va << " " << vb << endl;
    const double  phi = (sqrt(5.0) + 1.0) * 0.5;
    const double iphi = 1.0/phi;
    while (true){
      const double c = b-(b-a)*iphi;
      const double d = a+(b-a)*iphi;
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
  
  
  template<class summand>
  const auto sum(
      const summand& s,
      const size_t& n0,
      const size_t& N
      );

  template<class summand>
  const auto sum(
      const summand& s,
      const size_t& n0,
      const size_t& N
      ){
    if (N==n0) return s(N);
    return sum(s,n0,N-1)+s(N);
    return sum(s,n0+1,N)+s(n0);
  }
  
  template<class factor>
  const auto product(
      const factor& f,
      const size_t& n0,
      const size_t& N
      ){
    if (N==n0) return f(N);
    return product(f,n0,N-1)*f(N);
  }

  template<typename T>
  constexpr T binom(
      const size_t& n,
      const size_t& k){
    if (k>n)   return 0;
    if (n-k<k) return binom<T>(n,n-k);
    T v(1);
    for (size_t i=1;i<=k;++i){
      v*=(n+1-i);
      v/=i;
    }
    return v;
  }

  template<class series_summand>
  inline auto aitken_delta_squared(
      const series_summand& s,
      const size_t& N
      ){
    auto s0  = s(0);
    auto s1  = s(1)+s0;
    auto s2  = s(2)+s1;
    auto axn = s2;
    for (size_t n=1;;++n){
      if (abs(s2-s1)==0)           break;
      if (abs((s2-s1)-(s1-s0))<numeric_limits<double>::epsilon()) break;
      axn = s2-pow(s2-s1,2)/((s2-s1)-(s1-s0));
      //cerr << n << endl;
      if (n>=N) return axn;
      s0=s1;
      s1=s2;
      s2+=s(n+2);
    }
    return axn;
  }
 
  template<typename T>
  constexpr complex<T> riemann_zeta_subsummand(
      const complex<T>& s,
      const size_t& n,
      const size_t& k
      ){
    return binom<T>(n,k)*T(k%2==0?1:-1)*pow(complex<T>(k+1),-s);
  }

  template<typename T>
  constexpr complex<T> riemann_zeta_summand(
      const complex<T>& s,
      const size_t& n
      ){
    return sum(bind(riemann_zeta_subsummand<T>,s,n,_1),0,n)/T(pow(T(2),T(n+1)));
  }

  template<typename T>
  constexpr complex<T> riemann_zeta(
      const complex<T>& s,
      const size_t& N = 256
      ){
    return sum(bind(riemann_zeta_summand<T>,s,_1),0,N)
          /(complex<T>(1)-pow(T(2),complex<T>(1)-s));
  }

  template<class target>
  inline double linesearch(double c1,double v1,target function){
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

  double inline inverf( // approximation to 4.5e-4
      const double& t
      ){
    return t-((0.010328*t+0.802853)*t+2.515517)
           /(((0.001308*t+0.189269)*t+1.432788)*t+1.0);
  }

  uint64_t inline morton_encode(uint64_t x,uint64_t y){
    x = (x|(x<<16))&0x0000FFFF0000FFFF;
    x = (x|(x<< 8))&0x00FF00FF00FF00FF;
    x = (x|(x<< 4))&0x0F0F0F0F0F0F0F0F;
    x = (x|(x<< 2))&0x3333333333333333;
    x = (x|(x<< 1))&0x5555555555555555;
    y = (y|(y<<16))&0x0000FFFF0000FFFF;
    y = (y|(y<< 8))&0x00FF00FF00FF00FF;
    y = (y|(y<< 4))&0x0F0F0F0F0F0F0F0F;
    y = (y|(y<< 2))&0x3333333333333333;
    y = (y|(y<< 1))&0x5555555555555555;
    return x|(y<<1);
  }

  uint32_t inline morton_1(uint64_t x)
  {
    x = x & 0x5555555555555555;
    x = (x|(x>> 1))&0x3333333333333333;
    x = (x|(x>> 2))&0x0F0F0F0F0F0F0F0F;
    x = (x|(x>> 4))&0x00FF00FF00FF00FF;
    x = (x|(x>> 8))&0x0000FFFF0000FFFF;
    x = (x|(x>>16))&0xFFFFFFFFFFFFFFFF;
    return x;
  }

  tuple<uint64_t,uint64_t> const inline morton_decode(uint32_t c){
    return tuple<uint64_t,uint64_t>{morton_1(c),morton_1(c<<1)};
  }


  template<typename T>
  typename std::enable_if<std::is_integral<T>::value,T>::type
  constexpr factorial(const T n){
    return n<2?1:n*factorial(n-1);
  }

  const void inline quaternion_multiply(
      const double& a,
      const double& b,
      const double& c,
      const double& d,
      const double& e,
      const double& f,
      const double& g,
      const double& h,
      double& i,
      double& j,
      double& k,
      double& l){
    i = a*e-b*f-c*g-d*h;
    j = a*f+b*e+c*h-d*g;
    k = a*g-b*h+c*e+d*f;
    l = a*h+b*g-c*f+d*e;
  }

  const void inline quaternion_multiply_vector(
      const double& a,
      const double& b,
      const double& c,
      const double& d,
      const double& f,
      const double& g,
      const double& h,
      double& j,
      double& k,
      double& l){
    j = a*f+c*h-d*g;
    k = a*g-b*h+d*f;
    l = a*h+b*g-c*f;
  }

  const void inline vector_multiply_quaternion(
      const double& b,
      const double& c,
      const double& d,
      const double& e,
      const double& f,
      const double& g,
      const double& h,
      double& j,
      double& k,
      double& l){
    j =  b*e+c*h-d*g;
    k = -b*h+c*e+d*f;
    l = +b*g-c*f+d*e;
  }

  const void inline quaternion_rotate_vector(
      const double& a,
      const double& b,
      const double& c,
      const double& d,
      const double& f,
      const double& g,
      const double& h,
      double& j,
      double& k,
      double& l){
    quaternion_multiply_vector(a,b,c,d,f,g,h,j,k,l);
    double _j=j,_k=k,_l=l;
    vector_multiply_quaternion(_j,_k,_l,0.5*a,-0.5*b,-0.5*c,-0.5*d,j,k,l);
  }

  struct quaternion{
    double a,b,c,d;
    const quaternion operator*(const quaternion& o) const {
      quaternion r;
      quaternion_multiply(a,b,c,d,o.a,o.b,o.c,o.d,r.a,r.b,r.c,r.d);
      return r;
    }
  };

  array<double,9> constexpr quaternion_to_matrix(
      const double& q0,
      const double& q1,
      const double& q2,
      const double& q3){
    return 
    {
      q0*q0+q1*q1-q2*q2-q3*q3,         2*(q1*q2-q0*q3),         2*(q1*q3+q0*q2),
              2*(q1*q2+q0*q3), q0*q0-q1*q1+q2*q2-q3*q3,         2*(q2*q3-q0*q1),
              2*(q1*q3-q0*q2),         2*(q2*q3+q0*q1), q0*q0-q1*q1-q2*q2+q3*q3
    };
  }
  
  // triangle indices are hashed in a way to allow fast traversal in i
  // {{ 0 1 2 3 4 5 ... }= i
  // 0| - - - - - - ...
  // 1| 0 - - - - - ...
  // 2| 1 2 - - - - ...
  // 3| 3 4 5 - - - ...
  // 4| 6 7 8 9 - - ...
  // 5| a b c d e - ...
  // } = j
  
  constexpr size_t triangle_index_pack(const size_t& i,const size_t& j){
    return (j*(j-1))/2+i;
  }
  
  constexpr tuple<size_t,size_t> triangle_index_unpack(const size_t& n){
    const size_t j = (floor_square_root(8*n+1)+1)/2;
    const size_t i = n-(j*(j-1))/2;
    return {i,j};
  }
}

#endif // WMATH_MATH_H
