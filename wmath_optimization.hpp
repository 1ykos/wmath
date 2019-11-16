#ifndef WMATH_OPTIMIZATION_H
#define WMATH_OPTIMIZATION_H
#include "wmath_forward.hpp"
#include "wmath_bits.hpp"

namespace wmath{ 
  using dlib::abs;
  using dlib::identity_matrix;
  using dlib::is_finite;
  using dlib::length;
  using dlib::length_squared;
  using dlib::matrix;
  using dlib::normalize;
  using dlib::ones_matrix;
  using dlib::squared;
  using dlib::sum;
  using dlib::tmp;
  using dlib::trans;
  using dlib::zeros_matrix;
  //using dlib::is_nan;
  
  struct rosenbrock{
    rosenbrock(const double a = 1,const double b = 100): a(a),b(b) {}
    const double a;
    const double b;
    const double operator()(
        const dlib::matrix<double,0,1>& x
        ) const {
      return pow((a-x(0)),2)+b*pow(x(1)-pow(x(0),2),2);
    }
    const double operator()(
        const dlib::matrix<double,0,1>& x,
              dlib::matrix<double,0,1>& J
        ) const {
      cout << x(0) << " " << x(1) << " " << pow((a-x(0)),2)+b*pow(x(1)-pow(x(0),2),2) << endl;
      J(0) = -4*b*x(0)*(x(1)-pow(x(0),2))-2*(a-x(0));
      J(1) = 2*b*(x(1)-pow(x(0),2));
      return pow((a-x(0)),2)+b*pow(x(1)-pow(x(0),2),2);
    }
  };

  struct count_stop_strategy{
    size_t n = 0;
    size_t m = 0;
    count_stop_strategy(const size_t n,const size_t m): n(n),m(m) {};
    bool operator()(
        const size_t i,
        const size_t j,
        const double v0,
        const double v1) const {
      if (j>=m) return true;
      if (i>=n) return true;
      else      return false;
    }
  };

  template<class target>
  struct linesearch_target_wrapper{
    const matrix<double,0,1>& x0;
    const matrix<double,0,1>& dir;
    const target& f;
    linesearch_target_wrapper(
        const matrix<double,0,1>& x0,
        const matrix<double,0,1>& dir,
        const target& f
        ): x0(x0),dir(dir),f(f){}
    matrix<double,0,1> x(const double a) const {
      matrix<double,0,1> xa = x0+dir*a;   
      return xa;
    }
    const double operator()(const double a) const {
      return f(x(a));
    }
  };
  
  template<class target>
  struct linesearch_target_wrapper2{
    const matrix<double,0,1>& x0;
    const matrix<double,0,1>& dir1;
    const matrix<double,0,1>& dir2;
    const target& f;
    linesearch_target_wrapper2(
        const matrix<double,0,1>& x0,
        const matrix<double,0,1>& dir1,
        const matrix<double,0,1>& dir2,
        const target& f
        ): x0(x0),dir1(dir1),dir2(dir2),f(f){}
    const matrix<double,0,1> x(const double a) const {
      return tmp(x0+dir1*a*(1-a)+dir2*a*a);
    }
    const double operator()(const double a) const {
      //cerr << a*(1-a) << " " << 1-a*(1-a) << " " << f(x(a)) << endl;
      return f(x(a));
    }
  };

  template<class target,class stop_strategy>
  void newton_minimize(
      dlib::matrix<double,0,1>& x,
      const target& f,
      const stop_strategy& stop){
    double v0=f(x);
    double v1=v0;
    for (size_t i=0;;++i){
      v0 = v1;
      //cerr << v0 << " " << trans(x);
      const matrix<double,0,1> dx = -f.J(x)*pinv(f.H(x));
      //cerr << trans(dx);
      //linesearch_target_wrapper ltw_c(x,dx,f);
      //const matrix<double,0,1> dir_d = normalize(f.J(x))*length(dx);
      //linesearch_target_wrapper ltw_d(x,dir_d,f);
      //const double c = golden_section_search(0,ltw_c(0),1,ltw_c(1),ltw_c);
      //const double d = golden_section_search(0,ltw_d(0),1,ltw_d(1),ltw_d);
      //cerr << c               << " " << ltw(c) << " "
      //     << pow(0.75,min_j) << " " << ltw(pow(0.75,min_j)) << endl;
      //if (ltw_c(c)<ltw_d(d)) x+= c*ltw_c.dir;
      //else                   x+= d*ltw_d.dir;
      const matrix<double,0,1> dir1 = -normalize(f.J(x))*length(dx);
      const matrix<double,0,1> dir2 = dx;
      linesearch_target_wrapper2 ltw(x,dir1,dir2,f);
      //linesearch_target_wrapper ltw(x,dir2,f);
      const double a = golden_section_search(-1,ltw(-1),1,ltw(1),ltw);
      x=ltw.x(a);
      /*cerr << a << " " << ltw(a) << " "
           << ltw(-1) << " "
           << ltw( 0) << " "
           << ltw( 1) << endl;*/
      v1 = f(x);
      if (stop(i,v0,v1)) break;
    }
  }

  template<long n,class target,class stop_strategy>
  void find_min(
      matrix<double,n,1>& x0,
      const target& f,
      const stop_strategy& stop,
      const double epsilon = 1e-8
      ){
    const double  phi = (sqrt(5.0) + 1.0) * 0.5;
    const double iphi = 1.0/phi;
    matrix<double,n,n> H  = zeros_matrix<double>(x0.nr(),x0.nr());
    matrix<double,n,1> J0 = zeros_matrix<double>(x0.nr(),x0.nc());
    matrix<double,n,1> J1 = zeros_matrix<double>(x0.nr(),x0.nc());
    const matrix<double,0,0> I = zeros_matrix<double>(x0.nr(),x0.nr());
    double v0 = f(x0,J0);
    matrix<double,n,1> dx = -epsilon*normalize(J0);
    matrix<double,n,1> x1 = x0+dx;
    double alpha = 1;
    double last_min_length = epsilon;
    size_t j=0;
    for (size_t i=0;;++i){
      double v1 = f(x1,J1);
      const matrix<double,n,1> yk = J1-J0;
      const matrix<double,n,1> sk = dx;
      // BFGS update:
      //cerr << "BFGS update" << endl;
      matrix<double,n,n> U = yk*trans(yk)/(trans(yk)*sk);
      if (trans(sk)*H*sk>numeric_limits<double>::epsilon())
        U-=(H*sk*trans(sk)*trans(H))/(trans(sk)*H*sk);
      U = H+U;
      if (is_finite(U)){
        dlib::cholesky_decomposition<matrix<double,n,n>> chol(U);
        if (sum(squared(U-chol.get_l()*trans(chol.get_l())))/(U.nr()*U.nc())
            < numeric_limits<double>::epsilon()){
          if (v1<=v0){
            H=U;
            ++j; // number of proper, successfull BFGS updates
          } else {
            const double s = last_min_length/length(dx);
            H=(s*U+H)/(1+s);
            //cerr << "meh" << endl;
          }
        }
      }
      if (v1<v0) {
        last_min_length=length(dx);
        // bounded by 1 as this is the solution for a quadratic function and it
        // can only get better by accident
        // growing fast when alpha is small, growing slower when alpha is close
        // to 1
        // This procedure was tested to be better than alpha*=c for all factors c
        alpha=sqrt(alpha);
        swap(x0,x1);
        swap(v0,v1);
        swap(J0,J1);
      } else { 
        alpha*=exp(-1);
      }
      dx = -alpha*pinv(H)*J0;
      if (!(length(dx)>epsilon)) {
        dx = -epsilon*normalize(J0);
      }
      x1 = x0+dx;
      if (stop(i,j,v0,v1)) break;
    }
  }
};

#endif // WMATH_OPTIMIZATION_H
