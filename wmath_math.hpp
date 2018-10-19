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
    return min(max(n,l),h);
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
    if (n<0) return -roundup_2_3(n);
    if (n==0) return 1;
    if (n<4) return n;
    if (n%2==0) return 2*roundup_2_3(n/2);
    if (n%3==0) return 3*roundup_2_3(n/3);
    return roundup_2_3(n+1);
  }
  
  template <typename T>
  typename std::enable_if<std::is_integral<T>::value,T>::type
  constexpr roundup_2(const T& n){
    if (n<0) return -roundup_2(n);
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
