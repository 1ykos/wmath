#ifndef WMATH_DERIVATIVES_H
#define WMATH_DERIVATIVES_H

namespace wmath{

  /*
  template<typename R>
  struct add{
    constexpr R operator(const R& a,const R& b) {
      return a+b;
    }
    constexpr tuple<R,R> J(const R& a,const R& b) {
      return {1,1};
    }
  };
  
  template<typename R>
  struct mul{
    constexpr R operator(const R& a,const R& b) {
      return a*b;
    }
    constexpr tuple<R,R> J(const R& a,const R& b) {
      return {b,a};
    }
  };

  template<typename R>
  struct inv{
    constexpr R operator(const R& a) {
      return 1/a;
    }
    constexpr tuple<R> J(const R& a,const R& b) {
      return {-1/(a*a)};
    }
  };
  
  template<typename R>
  struct exp{
    constexpr R operator(const R& a) {
      return exp(a);
    }
    constexpr tuple<R> J(const R& a,const R& b) {
      return {exp(a)};
    }
  };

  TODO dammit this is hard
  */
}
