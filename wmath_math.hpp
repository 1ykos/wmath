#ifndef WMATH_MATH_H
#define WMATH_MATH_H
#include "wmath_forward.hpp"
#include "wmath_bits.hpp"

namespace wmath{
  template <typename T> int signum(T val) {
    return (T(0) < val) - (val < T(0));
  }

  template <typename T0,typename T1>
  typename std::enable_if<
      std::is_unsigned<T1>::value
    &&std::is_integral<T1>::value
    ,T0>::type
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
  typename std::enable_if<
      std::is_signed<T1>::value
    &&std::is_integral<T1>::value,
    double>::type
  constexpr pow(const T0& x,const T1& n){
    typedef typename make_unsigned<T1>::type U;
    return n<0?1.0/pow(x,U(-n)):pow(x,U(n));
  }

  // second order approximation to tetration
  constexpr double tetration(
      const double a,
      const double x
      )
  {
    if (x> 0) return pow(a,tetration(a,x-1));
    if (x>-1) return 1+(2*log(a)/(1+log(a)))*x+(1-log(a))/(1+log(a))*pow(x,2);
    else      return log(tetration(a,x+1))/log(a);
  }
  
  // use this to find good fourier transform dimensions
  template <typename T>
  typename std::enable_if<std::is_unsigned<T>::value,T>::type
  constexpr roundup_2_3(const T& n){
    if (n==0) return 1;
    if (n<4) return n;
    if (n%2==0) return 2*roundup_2_3(n/2);
    if (n%3==0) return 3*roundup_2_3(n/3);
    return roundup_2_3(n+1);
  }
  
  // round up to next power of two
  template <typename T>
  typename std::enable_if<std::is_unsigned<T>::value,T>::type
  constexpr roundup_2(const T& n){
    return (T(1)<<(log2(n-1)+1));
  }

  // these three really are a strange bunch, do they even make sense?
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
  typename std::enable_if<std::is_signed<T>::value,tuple<T,T>>::type
  constexpr extended_euclidean(const T a, const T b)
  {
    T r0 = a;
    T r1 = b;
    T s0 = 1;
    T s1 = 0;
    T t0 = 0;
    T t1 = 1;
    while (r1) {
      T q = r0 / r1;
      r0 = r0 - q*r1; swap(r0,r1);
      s0 = s0 - q*s1; swap(s0,s1);
      t0 = t0 - q*t1; swap(t0,t1);
    }
    // Bézout coefficients = (old_s,old_t)
    // gcd = old_r
    return tuple<T,T>({s0,t0});
  }

  template<typename T>
  typename std::enable_if<std::is_unsigned<T>::value,tuple<T,T>>::type
  constexpr extended_euclidean(const T a, const T b)
  {
    T r0 = a;
    T r1 = b;
    T s0 = 1;
    T s1 = 0;
    T t0 = 0;
    T t1 = 1;
    size_t n = 0;
    while (r1) {
      T q = r0 / r1;
      r0 = r0>q*r1?r0-q*r1:q*r1-r0; swap(r0,r1);
      s0 = s0+q*s1; swap(s0,s1);
      t0 = t0+q*t1; swap(t0,t1);
      ++n;
    }
    // Bézout coefficients = (old_s,old_t)
    // gcd = old_r
    if (n%2) s0=b-s0;
    else     t0=a-t0;
    return tuple<T,T>({s0,t0});
  }

  template<typename T>
  typename std::enable_if<std::is_unsigned<T>::value,T>::type
  modular_inverse(const T& n, const T& m) {
    if (popcount(m)<=1) return modular_inverse_power2(n,m);
    return get<0>(extended_euclidean(n,m));
  }

  // is faster when a<b ...
  // this works for negative integers, nice
  template <typename T>
  typename std::enable_if<std::is_integral<T>::value,T>::type
  //constexpr
  gcd(const T& a,const T&b){
    return a?gcd(b%a,a):b;
  }

  template <class It>
  typename iterator_traits<It>::value_type
  const gcd(const It begin,const It end){
    auto min = min_element(
        begin,end,
        [](const auto& a,const auto& b)
        {return a&&(abs(a)<abs(b));}
        );
    if (min==end) return 0;
    auto t = *min;
    for (auto it=begin;it!=end;++it) t = gcd(t,*it);
    return t;
  }

  template <typename T>
  typename std::enable_if<std::is_integral<T>::value,T>::type
  constexpr lcm(const T& a, const T& b){
    return (a&&b)?(a/gcd(a,b)*b):0;
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
  
/* power of x to e via exponentiation by squaring
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
  pow_mod(T x,T e,const T& m){
    T result(1);
    tuple<T,T> lm{m,0};
    const auto i = fixpoint_integer_inverse(lm);
    while (e){
      if (e&1ull) result = mul_mod(x,result,m,i);
      x=mul_mod(x,x,m,i);
      e>>=1;
    }
    return result;
  }
  
  /*
  template<typename T>
  typename std::enable_if<std::is_integral<T>::value,T>::type
  power_mod_inplace(T& x,T& e,const T& m) {
    T result(1);
    while (e){
      if (e&1ull) result=(x*result)%m;
      x=(x*x)%m;
      e>>=1;
    }
    return result;
  }
  */

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
 * start with T sumw=0.0, mean=0.0, var=0.0;
 * or with T sumw=w_1, mean=x_1, var=0;
 * and then call
 * mean_variance(x,w,sumw,mean,var)
 * for each pair of x and w */  
  template<typename T>
  typename enable_if<is_floating_point<T>::value>::type
  const inline mean_variance(
      const T &x,
      const T &w,
      T &sumw,
      T &mean,
      T &var
      ){
    sumw += w;
    const T delta = x-mean;
    mean += delta*w/sumw;
    var  += ((x-mean)*delta-var)*w/sumw;
  }
  
  template<typename T>
  typename enable_if<is_floating_point<T>::value>::type
  const inline mean_variance_undo(
      const T &x,
      const T &w,
      T &sumw,
      T &mean,
      T &var
      ){
    sumw -= w;
    const T delta = x-mean;
    if (abs(sumw)<numeric_limits<double>::min()) {
      mean  = 0;
      var   = 0;
    } else {
      mean -= delta*w/sumw;
      var  -= ((x-mean)*delta-var)*w/sumw;
    }
    var = var<0?0:var;
  }

/* numerically stable and incremental mean and variance
 * start with T sumw=0.0, mean=0.0, var=0.0;
 * or with T sumw=w_1, mean=x_1, var=0;
 * and then call
 * mean_variance(x,w,sumw,mean,var)
 * for each pair of x and w */  
  template<typename T>
  typename enable_if<is_floating_point<T>::value>::type
  const inline mean_variance(
      const T &x,
      const T &w,
      T &sumw,
      T &mean,
      T &var,
      T &sumw2
      ){
    sumw2+=pow(w,2u);
    mean_variance(x,w,sumw,mean,var);
  }

  /* unbiased estimator for variance using reliability weigths
   */
  template<typename T>
  inline const T variance(
      const T &var,  // biased estimator
      const T &sumw,
      const T &sumw2
      ) {
    return var*sumw/(sumw-sumw2/sumw);
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
      ) {
    T var = 0; // mean_variance does not depend on previous values of var
    mean_variance(x,w,sumw,mean,var);
  }
  
  template<typename T>
  typename enable_if<is_floating_point<T>::value,T>::type
  const inline mean_variance_sumw_dx(
      const T &x,
      const T &w,
      const T &sumw,
      const T &mean,
      const T &var
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
      const T &var
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
      const T &var
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
      const T &var
      ){
    return (x-mean)/sumw;
  }

  // TODO: this is wrong now
  template<typename T>
  typename enable_if<is_floating_point<T>::value,T>::type
  const inline mean_variance_var_dx_ini( // initial
      const T &x,
      const T &w,
      const T &sumw,
      const T &mean,
      const T &var
      ){
    return 2*w*(1-w/sumw)*(x-mean);
  }
  
  // TODO: this is wrong now
  template<typename T>
  typename enable_if<is_floating_point<T>::value,T>::type
  const inline mean_variance_var_dx_acc( // needs to be accumulated
      const T &x,
      const T &w,
      const T &w0,
      const T &sumw,
      const T &mean,
      const T &var
      ){
    return -2*w0*w*(x-mean)/sumw;
  }

  // TODO: this is wrong now
  template<typename T>
  typename enable_if<is_floating_point<T>::value,T>::type
  const inline mean_variance_var_dw_ini(
      const T &x,
      const T &w,
      const T &sumw,
      const T &mean,
      const T &var
      ){
    return (x-mean)*(x-mean);
  }


  // TODO: this is wrong now
//       (x0-mean)^2
// +2*w3*(mean-x0)*(x3-mean)/sumw
// +2*w2*(mean-x0)*(x2-mean)/sumw
// +2*w1*(mean-x0)*(x1-mean)/sumw
// +2*w0*(mean-x0)*(x0-mean)/sumw
// + ...
  template<typename T>
  typename enable_if<is_floating_point<T>::value,T>::type
  const inline mean_variance_var_dw_acc(
      const T &x,
      const T &w,
      const T &x0,
      const T &sumw,
      const T &mean,
      const T &var
      ){
    return 2*w/sumw*(mean-x0)*(x-mean);
  }

/* non corrected sample variance, correction factor is n*(n-1) 
 * to be used in conjunction with mean_variance */
  template<typename T>
  inline const T sample_variance(
      const T& var,
      const T& sumw
      ){
    return var;
  }

  template<typename T>
  struct mean_variance_calculator{
    T var =0;
    T sumw=0;
    size_t n=0;
    T mean=0;
    void push(const T& x, const T& w){
      ++n;
      mean_variance(x,w,sumw,mean,var);
    } 
    T variance() const {
      return wmath::variance(var,sumw,n);
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
 
  template<typename T>
  typename enable_if<is_floating_point<T>::value>::type
  const inline mean_covariance(
      const T &x,
      const T &y,
      const T &w,
      T &sw,
      T &mx,
      T &my,
      T &vx,
      T &vy,
      T &cv
      ){
    if (w<numeric_limits<T>::min()) return;
    sw+=w;
    const T dx = x-mx;
    const T dy = y-my;
    mx += dx*w/sw;
    vx += (dx*(x-mx)-vx)*w/sw;
    my += dy*w/sw;
    vy += (dy*(y-my)-vy)*w/sw;
    cv += (w*dx*(y-my)-cv)*w/sw;
  }
  
  template<typename T>
  typename enable_if<is_floating_point<T>::value>::type
  const inline mean_covariance_undo(
      const T &x,
      const T &y,
      const T &w,
      T &sw,
      T &mx,
      T &my,
      T &vx,
      T &vy,
      T &cv
      ){
    if (w<numeric_limits<T>::min()) return;
    sw-=w;
    const T dx = x-mx;
    const T dy = y-my;
    mx -= dx*w/sw;
    vx -= (dx*(x-mx)-vx)*w/sw;
    my -= dy*w/sw;
    vy -= (dy*(y-my)-vy)*w/sw;
    cv -= (w*dx*(y-my)-cv)*w/sw;
  }

/*
  template<typename T>
  struct cc_pearson{
    size_t n = 0; 
    T M2_x   = 0;
    T M2_y   = 0;
    T cov    = 0;
    T mean_x = 0;
    T mean_y = 0;
    T sumw   = 0;
    T sumw2  = 0;
    inline void push(const T& x,const T& y,const T& w) {
      sumw      += w;
      sumw2     += w*w;
      const T dx = x-mean_x;
      const T dy = y-mean_y;
      mean_x    += (w/sumw) * dx;
      mean_y    += (w/sumw) * dy;
      cov       += w*dy*dy;
      // TODO
    }
    T operator()() const {
      return cov/sqrt(M2_x*M2_y);
    }
  }

  // Pearsons product moment correlation coefficient
  template<typename T>
  class cc_pearson_{
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
*/

  template<typename T>
  class cc_fisher{
    public:
      size_t n = 0;
      T sumw   = 0;
      T var_x  = 0;
      T var_y  = 0;
      T covar  = 0;
      cc_fisher& push(const T& x, const T& y, const T& w=T(1)){
        // use this approximation: w = 1/(vx*vy + vx*y + vy*x)
        var_x+=w*x*x;
        var_y+=w*y*y;
        covar+=w*x*y;
        sumw +=w;
        n    +=1;
      }
      T operator()(){
        return covar/sqrt(var_x*var_y);
      }
      inline void clear(){
        var_x = 0;
        var_y = 0;
        covar = 0;
        n     = 0;
        sumw  = 0;
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
  
  template<class iterator>
  auto const destructive_median(
      iterator begin,
      iterator end
      )
  {
    nth_element(begin,begin+(end-begin  )/2,end);
    nth_element(begin,begin+(end-begin+1)/2,end);
    const auto median = (*(begin+(end-begin  )/2)
                        +*(begin+(end-begin+1)/2))/2;
    return median;
  }
  
  template<class iterator>
  auto const destructive_median_mad(
      iterator begin,
      iterator end
      )
  {
    nth_element(begin,begin+(end-begin  )/2,end);
    nth_element(begin,begin+(end-begin+1)/2,end);
    const auto median = (*(begin+(end-begin  )/2)
                        +*(begin+(end-begin+1)/2))/2;
    transform(begin,end,begin,[&median](const auto v){return abs(v-median);});
    nth_element(begin,begin+(end-begin  )/2,end);
    nth_element(begin,begin+(end-begin+1)/2,end);
    const auto mad    = (*(begin+(end-begin  )/2)
                        +*(begin+(end-begin+1)/2))/2;
    return make_tuple(median,mad);
  }

  template<class iterator>
  auto const destructive_weighted_median_variance(
      iterator begin,
      iterator end
      )
  {
    auto median = get<0>(*begin); median = 0;
    if (begin==end) {
      return make_tuple(median,0.0);
    }
    sort(begin,end);
    auto lower_sum = get<1>(*begin);lower_sum = 0;
    auto upper_sum = get<1>(*begin);upper_sum = 0;
    auto sumw2     = get<1>(*begin);sumw2 = 0;
    auto it0 = begin;
    auto it1 = end-1;
    while (it0!=it1) {
      if (lower_sum<upper_sum) {
        lower_sum += get<1>(*it0);
        sumw2 += pow(get<1>(*it0),2);
        ++it0;
        continue;
      }
      if (upper_sum<lower_sum) {
        upper_sum += get<1>(*it1);
        sumw2 += pow(get<1>(*it1),2);
        --it1;
        continue;
      }
      lower_sum += get<1>(*it0);
      ++it0;
    }
    if (upper_sum<lower_sum) {
      median = get<0>(*it0);
    } else {
      if (lower_sum<upper_sum) {
        median = get<0>(*it1);
      } else {
        median = (get<0>(*it0)+get<1>(*it1))/2;
      }
    }
    const auto sumw = lower_sum+upper_sum;
    const double is2 = 4.0*pow(sumw,2)/sumw2;
    double variance = 0;
    double normalization = 0;
    for (it0=begin;it0!=end;++it0) {
      const double w = exp(-0.5*pow(lower_sum/sumw-0.5,2)*is2)
                      *get<1>(*it0);
                      //*pow(get<1>(*it0),3);
      variance += w*pow(get<0>(*it0)-median,2);
      normalization+=w;
    }
    variance/=normalization;
    return make_tuple(median,variance);
  }

  // copied from c++ reference to test
  template<class ForwardIt, class UnaryPredicate>
  ForwardIt _partition(ForwardIt first, ForwardIt last, UnaryPredicate p)
  {
    first = std::find_if_not(first, last, p);
    if (first == last) return first;
    for (auto i = std::next(first); i != last; ++i)
    {
      if (p(*i))
      {
        //std::iter_swap(*i, first);
        swap(*i,*first);
        ++first;
      }
    }
    return first;
  }
  
  // Dijkstra's three way partitioning
  template<class iterator,class threewaycomparator>
  tuple<iterator,iterator> three_way_partition(
      iterator begin,
      iterator end,
      threewaycomparator comp
      ) {
    auto next = begin;
    while (next!=end) {
      if (comp(*next)<0) {
        std::iter_swap(begin++,next++);
        continue;
      }
      if (comp(*next)>0) {
        std::iter_swap(next,--end);
      } else {
        ++next;
      }
    }
    return {begin,next};
  }
  
  template<class iterator,typename T,bool lower>
  auto lower_upper_weighted_quartile_element(
      iterator begin,
      T quartile,
      iterator end
      )
  {
    auto sum = T(0);
    const auto _begin = begin;
    while (begin!=end) {
      auto pivot = get<0>(*(begin+(end-begin)/2));
      auto [middle0,middle1] =
          three_way_partition(begin,end,
              [pivot](const auto& v){
              return (pivot<get<0>(v))-(get<0>(v)<pivot);
              });
      T lower_sum = 0;
      for (auto it=begin;it!=middle0;++it) {
        lower_sum+=get<1>(*it);
        if constexpr (lower) {
          if ((sum+lower_sum)>=quartile) {
            end=middle0;
            break;
          }
        } else {
          if ((sum+lower_sum)> quartile) {
            end=middle0;
            break;
          }
        }
      }
      if constexpr (lower) {
        if ((sum+lower_sum)>=quartile) continue;
      } else {
        if ((sum+lower_sum)> quartile) continue;
      }
      sum+=lower_sum;
      begin = middle0;
      for (auto it=begin;it!=middle1;++it) {
        sum+=get<1>(*it);
        if constexpr (lower) {
          if (sum>=quartile) return it;
        } else {
          if (sum> quartile) return it;
        }
      }
      if (middle1==end) return begin;
      begin = middle1;
    }
    return begin;
  }
  
  template<class iterator,typename T>
  auto lower_weighted_quartile_element(
      iterator begin,
      T quartile,
      iterator end
      )
  {
    return lower_upper_weighted_quartile_element<iterator,T,true>(
        begin,quartile,end);
  } 
  
  template<class iterator,typename T>
  auto upper_weighted_quartile_element(
      iterator begin,
      T quartile,
      iterator end
      )
  {
    return lower_upper_weighted_quartile_element<iterator,T,false>(
        begin,quartile,end);
  } 
  
  template<class iterator,typename T>
  auto weighted_quartile_value(
      iterator begin,
      T q,
      iterator end
      )
  {
    if (begin==end) {
      typename remove_reference<decltype(get<0>(*begin))>::type d = 0;
      return d;
    }
    auto element_hi = upper_weighted_quartile_element(begin,q,end);
    auto upper_value = get<0>(*element_hi);
    auto element_lo = lower_weighted_quartile_element(begin,q,element_hi+1);
    auto lower_value = get<0>(*element_lo);
    return lower_value+(upper_value-lower_value)/2;
  }
  
  template<class iterator,typename T>
  auto const primitive_weighted_quartile_value(
      iterator begin,
      T q,
      iterator end
      )
  {
    sort(begin,end,
        [](const auto& a,const auto &b){return get<0>(a)<get<0>(b);});
    T sum = 0;
    auto it = begin;
    for (;it!=end;++it) {
      sum+=get<1>(*it);
      if (sum>q) return get<0>(*it);
      if (sum<q) continue;
      const auto lower_value = get<0>(*it);
      for (++it;it!=end;++it) {
        sum+=get<1>(*it);
        if (sum>q) {
          const auto upper_value = get<0>(*it);
          return lower_value+(upper_value-lower_value)/2;
        }
      }
      return lower_value;
    }
    if (begin==end) {
      typename remove_reference<decltype(get<0>(*begin))>::type d = 0;
      return d;
    } else {
      --it;
      return get<0>(*it);
    }
  }
  
  template<class iterator,typename T>
  auto const primitive_lower_weighted_quartile_element(
      iterator begin,
      T q,
      iterator end
      )
  {
    if (begin==end) return begin;
    sort(begin,end,
        [](const auto& a,const auto &b){return get<0>(a)<get<0>(b);});
    T sum = 0;
    auto it = begin;
    for (;it!=end;++it) {
      sum+=get<1>(*it);
      if (sum>=q) return it;
    }
    return --it;
  }
  
  template<class iterator,typename T>
  auto const primitive_upper_weighted_quartile_element(
      iterator begin,
      T q,
      iterator end
      )
  {
    if (begin==end) return begin;
    sort(begin,end,
        [](const auto& a,const auto &b){return get<0>(a)<get<0>(b);});
    T sum = 0;
    auto it = begin;
    for (;it!=end;++it) {
      sum+=get<1>(*it);
      if (sum>q) return it;
    }
    return --it;
  }

  /* TODO
  template<typename T>
  typename std::enable_if<std::is_unsigned<T>::value,
           tuple<tuple<T,T>,tuple<T,T>>>::type
  long_division(
      const tuple<T,T>& n,
      const tuple<T,T>& d,
      )
  {
    tuple<T,T> x{0,0};
    auto d1=d;
    auto tmp=get<1>(d1);
    --get<1>(d1);
    get<0>(d1)-=get<1>(d1)>tmp;
    if (get<0>(d1)) {
      get<1>(x)=T(1)<<clz(get<0>(d1));
    } else {
      get<0>(x)=T(1)<<clz(get<1>(d1));
    }
    for (size_t j=0;j!=log2(digits<T>())*2;++j) {
      auto x1=x;
      auto tmp=get<1>(x1);
      ++get<1>(x1);
      get<0>(x1)+=get<1>(d1)<tmp;
      auto xd=long_mul(x,d);
      get<0>(x)=~get<0>(xd);
      get<1>(x)=~get<1>(xd);
      tmp = get<1>(xd);
      --get<1>(xd);
      get<1>(xd)-=get<1>(xd)>tmp;
      x += super_long_mul(x1,xd);
    }
    tuple<T,T> q = get<0>(super_long_mul(n,xd));
    if (n-long_mul(q,d)>=d) ++q;
    tuple<T,T> r = n-long_mul(q,d);
    return {q,r};
  }
  */
}

#endif // WMATH_MATH_H
