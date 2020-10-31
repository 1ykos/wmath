#ifndef WMATH_OPTIMIZATION_H
#define WMATH_OPTIMIZATION_H
#include "dlib/matrix.h"
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
  using dlib::randm;
  using dlib::squared;
  using dlib::sum;
  using dlib::pointwise_multiply;
  using dlib::dot;
  using dlib::tmp;
  using dlib::trans;
  using dlib::zeros_matrix;
  //using dlib::is_nan;


  template<class vector>
  struct null_projection{
    const void operator()(const vector& x) const {return;}
  };
  
  struct box_projection{
    matrix<double,0,2> box;
    const void operator()(matrix<double,0,1>& x)
    {
      for (size_t i=0;i!=box.nr();++i) x(i) = clip(x(i),box(i,0),box(i,1));
    }
  };

  struct rosenbrock{
    rosenbrock(const double a = 1,const double b = 100): a(a),b(b) {}
    const double a;
    const double b;
    const double operator()(
        const dlib::matrix<double,2,1>& x
        ) const {
      cout << pow((a-x(0)),2)+b*pow(x(1)-pow(x(0),2),2) << " "
           << x(0) << " " << x(1) << endl;
      return pow((a-x(0)),2)+b*pow(x(1)-pow(x(0),2),2);
    }
    const double operator()(
        const dlib::matrix<double,2,1>& x,
              dlib::matrix<double,2,1>& J
        ) const {
      // cout << x(0) << " "
      //      << x(1) << " "
      //      << pow((a-x(0)),2)+b*pow(x(1)-pow(x(0),2),2) << endl;
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
      const stop_strategy& stop
      ){
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

  template<long n,class function>
  const inline matrix<double,n,1>
  numerical_derivative(
      matrix<double,n,1>& x0,
      const matrix<double,n,1>& epsilon,
      const function& f
      )
  {
    matrix<double,n,1> J;
    matrix<double,n,1> x = x0;
    for (size_t i=0;i!=x.nr();++i) {
      x(i) = x0(i) + epsilon(i);
      J(i) = function(x);
      x(i) = x0(i) - epsilon(i);
      J(i)-= function(x);
      J(i)/= 2*epsilon(i);
    }
    return J;
  }

  template<long n,class target,class stop_strategy,class projection>
  void inline find_min_gradient(
      matrix<double,n,1>& x0,
      const target& f,
      const projection& p = null_projection<matrix<double,n,1>>{},
      const stop_strategy& stop = count_stop_strategy{4096,4096},
      const double epsilon = 1e-8
      ){
    const double  phi = (sqrt(5.0) + 1.0) * 0.5;
    const double iphi = 1.0/phi;
    const auto v0_J0 = f(x0);
    double v0 = get<0>(v0_J0);
    matrix<double,n,1> J0 = get<1>(v0_J0);
    matrix<double,n,1> dx;
    matrix<double,n,1> x1;
    double beta = epsilon;
    size_t j=0;
    for (size_t i=0;;++i) {
      while (length(J0)<numeric_limits<double>::min()) {
        x1 = x0+epsilon*normalize(randm(x0.nr(),1));
        const auto v1_J1 = f(x1);
        double v1 = get<0>(v1_J1);
        matrix<double,n,1> J1 = get<1>(v1_J1); 
        if ((v1<v0)&&(!isnan(v1))) {
          swap(x0,x1);
          swap(v0,v1);
          swap(J0,J1);
          beta = epsilon;
        }
        if (stop(i,j,v0,v1)) return;
      }
      dx = -beta*normalize(J0);
      x1 =  x0+dx;
      const auto v1_J1 = f(x1);
      double v1 = get<0>(v1_J1);
      matrix<double,n,1> J1 = get<1>(v1_J1); 
      if ((v1<v0)&&(!isnan(v1))) {
        beta*=2;
        swap(x0,x1);
        swap(v0,v1);
        swap(J0,J1);
      } else {
        if (beta<0.5*epsilon) break;
        beta*=iphi;
      }
      dx = -beta*normalize(J0);
      x1 = x0+dx;
      p(x1);
      if (stop(i,j,v0,v1)) return;
    }
  }
  
  template<long n,class target,class stop_strategy>
  void inline find_min_gradient(
      matrix<double,n,1>& x0,
      const target& f,
      const stop_strategy& stop = count_stop_strategy{4096,4096},
      const double epsilon = 1e-8
  ){find_min_gradient(x0,f,null_projection<matrix<double,n,1>>{},stop,epsilon);}

  template<long n,class target,class stop_strategy,class projection>
  void inline find_min(
      matrix<double,n,1>& x0,
      const target& f,
      const projection& p,
      const stop_strategy& stop = count_stop_strategy{4096,4096},
      const double epsilon = 1e-8
      ){
    constexpr size_t verbosity = 0;
    const double  phi = (sqrt(5.0) + 1.0) * 0.5;
    const double iphi = 1.0/phi;
    matrix<double,n,n> H  = zeros_matrix<double>(x0.nr(),x0.nr());
    const matrix<double,n,n> I = zeros_matrix<double>(x0.nr(),x0.nr());
    const auto fx0 = f(x0);
    auto v0 = get<0>(fx0);
    auto J0 = get<1>(fx0);
    matrix<double,n,1> dx = -epsilon*normalize(J0);
    matrix<double,n,1> x1 = x0+dx;
    double alpha = 1;
    double beta = epsilon;
    double last_min_length = epsilon;
    size_t j=0;
    size_t fallback = 0;
    for (size_t i=0;;++i) {
      const auto fx1 = f(x1);
      auto v1 = get<0>(fx1);
      auto J1 = get<1>(fx1);
      const matrix<double,n,1> yk = J1-J0;
      const matrix<double,n,1> sk = dx;
      // BFGS update:
      matrix<double,n,n> U = yk*trans(yk)/(trans(yk)*sk);
      if (trans(sk)*H*sk>numeric_limits<double>::epsilon())
        U-=(H*sk*trans(sk)*trans(H))/(trans(sk)*H*sk);
      U+=H;
      if (is_finite(U)) {
        dlib::cholesky_decomposition<matrix<double,n,n>> chol(U);
        if (sum(squared(U-chol.get_l()*trans(chol.get_l())))/(U.nr()*U.nc())
            < numeric_limits<double>::epsilon()){
          if (v1<=v0) {
            //cerr << "simple BFGS update" << endl;
            H=U;
            ++j; // number of proper, successfull BFGS updates
          } else {
            //cerr << "modified BFGS update" << endl;
            const double s = last_min_length/length(dx);
            H=(s*U+H)/(1+s);
            //cerr << "meh" << endl;
          }
        }
      }
      if ((v1<v0)&&(!isnan(v1))) {
        last_min_length=length(dx);
        // bounded by 1 as this is the solution for a quadratic function and it
        // can only get better by accident
        // growing fast when alpha is small, growing slower when alpha is close
        // to 1
        // This procedure was tested to be better than alpha*=c for all factors c
        // beta = abs(trans(x1-x0)*(J1-J0))/length_squared(J1-J0);
        if constexpr (verbosity>0) {
          if (!fallback)
            cerr << "successfull BFGS step with alpha = " << alpha << endl;
          else
            cerr << "successfull gradient descent step with beta = " << beta << endl;
        }
        if (fallback) {
          beta*=2;
          //alpha=1;//(1+15*alpha)/16;
          //beta = abs(length(J0)/length(J1-J0));
        } else {
          alpha=sqrt(alpha);
        }
        swap(x0,x1);
        swap(v0,v1);
        swap(J0,J1);
      } else {
        if (beta<0.5*epsilon) break;
        if (fallback) {
          beta*=iphi;
          //beta = abs(length(J0)/length(J1-J0));
        } else {
          alpha*=exp(-1);
        }
      }
      if ((v1>=v0)&&(!fallback)) {
        dx = -beta*normalize(J0);
        ++fallback;
        if constexpr (verbosity>0) {
          cerr << "alternate with gradient descent using beta = " << beta
               << " alpha = " << alpha << endl;
        }
      } else {
        //cerr << "computing BFGS step, alpha = " << alpha << endl;
        //dx = -alpha*pinv(H)*J0;
        //dlib::cholesky_decomposition<matrix<double,n,n>> chol(H*trans(H));
        //dlib::qr_decomposition<matrix<double,n,n>> QR(H*trans(H));
        //dx = dlib::qr_decomposition<matrix<double,n,n>>(H).solve(-alpha*J0);
        if (fallback>1) {
          //dx = trans(H)
          //  *dlib::qr_decomposition<matrix<double,n,n>>(H*trans(H)).solve(J0);
          dx =
             dlib::qr_decomposition<matrix<double,n,n>>(trans(H)*H)
                                  .solve(trans(H)*J0);
          if (is_finite(dx)) alpha = beta/length(dx);
          dx*=-alpha;
        } else {
          dx = -alpha*trans(H)
             *dlib::qr_decomposition<matrix<double,n,n>>(H*trans(H)).solve(J0);
        }
        //dx = -alpha*trans(H)
        //     *dlib::cholesky_decomposition<matrix<double,n,n>>(H*trans(H))
        //                                                      .solve(J0);
        //cerr << "BFGS step computed" << endl;
        //cerr << trans(dx) << endl;
        //cerr << trans(dx) << trans(-alpha*pinv(H)*J0);
        if ((!(length(dx)>min(epsilon,beta*exp(-1))))|(!is_finite(dx))) {
          if (length(dx)<=epsilon*exp(1)) alpha*=exp(1);
          dx = -beta*normalize(J0);
          ++fallback;
          if constexpr (verbosity>0) {
            cerr << "fallback to gradient descent with beta = " << beta << endl;
          }
        } else {
          fallback = 0;
        }
      }
      x1 = x0+dx;
      p(x1);
      dx = x1-x0;
      if (stop(i,j,v0,v1)) break;
    }
  }
  
  template<long n,class target,class stop_strategy>
  void inline find_min(
      matrix<double,n,1>& x0,
      const target& f,
      const stop_strategy& stop = count_stop_strategy{4096,4096},
      const double epsilon = 1e-8
  ){find_min(x0,f,null_projection<matrix<double,n,1>>{},stop,epsilon);}
  
  template<long n>
  class second_order_model{
    private:
      const tuple<double,matrix<double,n,1>,matrix<double,n,n>>
      project_a_J_H
      (
        const matrix<double,n,1> dx,
        const double dv
      ) const
      {
        const double l2 = length_squared(dx);
        const double a = dv/(1+l2+pow(l2,2)/2);
        const matrix<double,n,1>& dJ = a*dx;
        const matrix<double,n,n>& dH = a*dx*trans(dx);
        return {a,dJ,dH};
      }
      const tuple<matrix<double,n,1>,matrix<double,n,n>>
      project_J_H
      (
        const matrix<double,n,1> dx,
        const double dv
      ) const
      {
        const double l2 = length_squared(dx);
        const double a = dv/(l2+pow(l2,2)/2);
        const matrix<double,n,1>& dJ = a*dx;
        const matrix<double,n,n>& dH = a*dx*trans(dx);
        return {dJ,dH};
      }
    public:
    double lag = exp(1)*n;
    double v0 = 0;
    matrix<double,n,1> x0 = zeros_matrix<double>(n,1); // point of development
    // covariance matrix of trust region
    matrix<double,n,n> S  = zeros_matrix<double>(n,n);
    matrix<double,n,1> J  = zeros_matrix<double>(n,1);
    matrix<double,n,n> H  = zeros_matrix<double>(n,n);
    second_order_model(const long m)
    {
      x0 = zeros_matrix<double>(m,1);
      S  = zeros_matrix<double>(m,m);
      J  = zeros_matrix<double>(m,1);
      H  = zeros_matrix<double>(m,m);
      lag= exp(1)*m;
    }
    double operator() (
        const matrix<double,n,1>& dx,
        const double& da,
        const matrix<double,n,1>& dJ,
        const matrix<double,n,n>& dH
        )const {
      return da+dot(dx,dJ)+0.5*dot(dx,dH*dx);
    }
    double operator()(const matrix<double,n,1>& x)const {
      return (*this)(x-x0,v0,J,H);
    }
    const
    tuple<double,double,matrix<double,n,1>,matrix<double,n,n>>
    inline update_length_squared(
        const double& da,
        const matrix<double,n,1>& dJ,
        const matrix<double,n,n>& dH
    ) const
    {
      double v = 0;
      for (size_t i=0;i!=dH.nr();++i)
        for (size_t j=0;j!=dH.nc();++j)
          v+=(S*dH)(i,j)*(S*dH)(i,j);
      v*=3*0.25;
      for (size_t i=0;i!=dH.nr();++i)
        for (size_t j=0;j!=dH.nc();++j)
          v+=da*(S(i,j)*(dH*S)(i,j));
      v+=dot(dJ,S*dJ)+pow(da,2);
      double diff_da = 2*da;
      for (size_t i=0;i!=dH.nr();++i) // dot(S,dH*S)
        for (size_t j=0;j!=dH.nc();++j)
          diff_da += S(i,j)*(dH*S)(i,j);
      matrix<double,n,1> diff_dJ = 2*S*dJ;
      matrix<double,n,n> diff_dH = 3*0.5*S*S*dH+da*S*S;
      return {v,diff_da,diff_dJ,diff_dH};
    }
    const
    tuple<double,double,matrix<double,n,1>,matrix<double,n,n>>
    inline weak_target(
        const matrix<double,n,1>& dx,
        const double& dv,
        const double& da,
        const matrix<double,n,1>& dJ,
        const matrix<double,n,n>& dH
    ) const
    {
      auto [v,diff_da,diff_dJ,diff_dH] = update_length_squared(da,dJ,dH);
      //return {v,diff_da,diff_dJ,diff_dH};
      const double ddv = dv-(*this)(dx,da,dJ,dH);
      /*v*=lag;
      diff_da*=lag;
      diff_dJ*=lag;
      diff_dH*=lag;*/
      v+=pow(ddv,2);
      diff_da-=2*da*ddv;
      diff_dJ-=2*dx*ddv;
      diff_dH-=dx*trans(dx)*ddv;
      return {v,diff_da,diff_dJ,diff_dH};
    }
    void test_weak_target() const {
      std::minstd_rand rng{std::random_device{}()};
      std::normal_distribution<double> nd;
      matrix<double,n,1> J = zeros_matrix<double>(x0.nr(),1);
      matrix<double,n,3> S = zeros_matrix<double>(x0.nr(),x0.nr());
      matrix<double,n,3> H = zeros_matrix<double>(x0.nr(),x0.nr());
      matrix<double,n,1> dx;
      for (size_t i=0;i!=x0.nr();++i) {
        dx(i)=nd(rng);
        J(i) =nd(rng);
      }
      for (size_t i=0;i!=x0.nr();++i) {
        for (size_t j=0;j!=x0.nr();++j) {
          H(i,j)=nd(rng);
          S(i,j)=nd(rng);
        }
      }
      //H=trans(H)*H;
      S=trans(S)*S;
      double a = 3;
      double dv = 5;
      const auto [v,diff_da,diff_dJ,diff_dH] = weak_target(dx,dv,a,J,H);
      const double eps = 1e-8;
      cout << a << endl;
      cout << diff_da << " " <<
        (get<0>(weak_target(dx,dv,a+eps,J,H))
        -get<0>(weak_target(dx,dv,a-eps,J,H)))/(2*eps) << endl;
      for (size_t i=0;i!=3;++i) {
        J(i)+=eps;
        const double p = get<0>(weak_target(dx,dv,a,J,H));
        J(i)-=2*eps;
        const double m = get<0>(weak_target(dx,dv,a,J,H));
        cout << diff_dJ(i) << " " << (p-m)/(2*eps) << endl;
        J(i)+=eps;
      }
      for (size_t i=0;i!=3;++i) for (size_t j=0;j!=3;++j) {
        H(i,j)+=eps;
        const double p = get<0>(weak_target(dx,dv,a,J,H));
        H(i,j)-=2*eps;
        const double m = get<0>(weak_target(dx,dv,a,J,H));
        cout << diff_dH(i,j) << " " << (p-m)/(2*eps) << endl;
        H(i,j)+=eps;
      }
    }
    void weak_update(const matrix<double,n,1>& x, const double v) {
      //min(
      // least change of previous model
      // lag*(3/4*dot(S*dH,S*dH)+da*dot(S,dH*S)+dot(dJ,S*dJ)+pow(da,2))
      // new function value
      // pow(dv-(da+dot(dJ,dx)+dot(dx,H*dx)/2),2) )
      const double  phi = (sqrt(5.0) + 1.0) * 0.5;
      const double iphi = 1.0/phi;
      const matrix<double,n,1> dx = x-x0;
      //cout << (*this)(x) << endl;
      const double dv0 = v-(*this)(x);
      auto [_dJ,_dH] = project_J_H(dx,dv0);
      double stepsize = length_squared(_dJ);
      for (size_t i=0;i!=_dH.nr();++i)
        for (size_t j=0;j!=_dH.nc();++j)
          stepsize+=_dH(i,j)*_dH(i,j);
      stepsize = sqrt(stepsize);
      matrix<double,n,1> dJ   = zeros_matrix<double>(S.nr(),1);
      matrix<double,n,n> dH   = zeros_matrix<double>(S.nr(),S.nc());
      const double da = 0; // don't change v0
      double best_value = get<0>(weak_target(dx,dv0,da,dJ,dH));
      //cerr << "begin optimisation" << endl;
      //cerr << pow(dv0,2) << " " << best_value << " " << dv0 << endl;
      for (size_t i=0;i!=256;) {
        const auto [ignore,diff_da,diff_dJ,diff_dH] =
          weak_target(dx,dv0,da,dJ,dH);
        double difflength = length_squared(diff_dJ);
        for (size_t i=0;i!=dH.nr();++i)
          for (size_t j=0;j!=dH.nc();++j)
            difflength+=diff_dH(i,j)*diff_dH(i,j);
        difflength = 1.0/difflength;
        while (i!=256) {
          const double value = get<0>(weak_target(
                dx,
                dv0,
                da,
                dJ-difflength*stepsize*diff_dJ,
                dlib::make_symmetric(dH-difflength*stepsize*diff_dH)
                //dH-difflength*stepsize*diff_dH
                ));
          //cerr << value << " " << stepsize << " " << i << endl;
          if (value<best_value) {
            best_value = value;
            //cerr << "update" << endl;
            dJ-=difflength*stepsize*diff_dJ;
            dH=make_symmetric(dH-difflength*stepsize*diff_dH);
            //dH-=difflength*stepsize*diff_dH;
            stepsize*=2;
            break;
          } else {
            stepsize*=iphi;
          }
          ++i;
        }
      }
      //cerr << "supposed to be:" << v << endl;
      //cerr << "before optimisation: " << (*this)(x) << endl;
      J+=dJ;
      H+=dH;
      //cerr << "after optimisation: " << (*this)(x) << endl;
      S=S*lag/(lag+1)+dx*trans(dx)/(lag+1);
    }
    void strong_update(const matrix<double,n,1>& x, const double v) {
      const double  phi = (sqrt(5.0) + 1.0) * 0.5;
      const double iphi = 1.0/phi;
      const matrix<double,n,1> dx = x-x0;
      const double dv0 = v-(*this)(x);
      auto [dJ,dH] = project_J_H(dx,dv0);
      double stepsize = length_squared(dJ);
      for (size_t i=0;i!=dH.nr();++i)
        for (size_t j=0;j!=dH.nc();++j)
          stepsize+=dH(i,j)*dH(i,j);
      stepsize = sqrt(stepsize);
      const double da = 0; // don't change v0
      double best_value = get<0>(update_length_squared(da,dJ,dH));
      //cerr << "begin optimisation" << endl;
      //cerr << pow(dv0,2) << " " << best_value << " " << dv0 << endl;
      for (size_t i=0;i!=256;) {
        auto [ignore,diff_da,diff_dJ,diff_dH] =
          update_length_squared(da,dJ,dH);
        diff_dJ=-diff_dJ-dot(-diff_dJ,dx)*dx/length_squared(dx);
        diff_dH=-
          diff_dH - sum(pointwise_multiply(-diff_dH,dx*trans(dx)))*x*trans(x)
                     * pow(length_squared(dx),-2);
        double difflength = length_squared(diff_dJ);
        for (size_t i=0;i!=dH.nr();++i)
          for (size_t j=0;j!=dH.nc();++j)
            difflength+=diff_dH(i,j)*diff_dH(i,j);
        difflength = 1.0/difflength;
        while (i!=256) {
          auto [ddJ,ddH] = project_J_H(dx,v-(*this)(dx)-(*this)(
                dx,
                da,
                dJ+stepsize*difflength*diff_dJ,
                dH+stepsize*difflength*diff_dH
                ));
          const double value = get<0>(update_length_squared(da,dJ+ddJ,dH+ddH));
          //cerr << value << " " << stepsize << " " << i << endl;
          if (value<best_value) {
            best_value = value;
            //cerr << "update" << endl;
            dJ+=ddJ;
            dH=make_symmetric(dH+ddH);
            stepsize*=2;
            break;
          } else {
            stepsize*=iphi;
          }
          ++i;
        }
      }
      //cerr << "supposed to be:" << v << endl;
      //cerr << "before optimisation: " << (*this)(x) << endl;
      J+=dJ;
      H+=dH;
      //cerr << "after optimisation: " << (*this)(x) << endl;
      S=S*lag/(lag+1)+dx*trans(dx)/(lag+1);
    }
    void set_zero(const matrix<double,n,1>& x,const double& v) {
      v0 = v;
      matrix<double,n,1> dx = x-x0;
      S=S*lag/(lag+1)+dx*trans(dx)/(lag+1);
      x0 = x;
    }
    matrix<double,n,1> minimise(double alpha,double beta,size_t& flag) {
      //matrix<double,n,1> dx = trans(H)
      //      *dlib::qr_decomposition<matrix<double,n,n>>(H*trans(H)).solve(J);
      matrix<double,n,1> dx =
      dlib::qr_decomposition<matrix<double,n,n>>(trans(H)*H).solve(trans(H)*J);
      if (is_finite(dx)) {
        flag = 1;
        return -alpha*normalize(dx);
      }
      //cerr << "# second derivative failed " << endl;
      dx = -beta*normalize(J);
      flag = 2;
      return dx;
    }
  };

  template<
    long n,
    class target,
    class proposer,
    class projection,
    class stop_strategy
  >
  void inline find_min_numerical(
      matrix<double,n,1>& x0,
      const target& f,
            proposer&   propose,
      const projection& project,
      const stop_strategy& stop = count_stop_strategy{4096,4096}
      ){
    const double  phi = (sqrt(5.0) + 1.0) * 0.5;
    const double iphi = 1.0/phi;
    double v0 = f(x0);
    second_order_model<n> model(x0.nr());
    model.set_zero(x0,v0);
    matrix<double,n,1> x1 = propose(x0);
    double v1 = f(x1);
    model.strong_update(x1,v1);
    size_t flag = 0;
    double alpha = 1;
    double beta = 1e-20;
    for (size_t i=0;;++i)
    {
      if (i&1u) {
        flag = 0;
        x1 = model.minimise(alpha,beta,flag);
        if (isnan(length(x1))) {
          flag = 0;
          x1 = propose(x0);
          v1 = f(x1);
        } else {
          x1+= x0;
          v1 = f(x1);
          if (isnan(v1)) {
            flag = 0;
            x1 = propose(x0);
            v1 = f(x1);
          }
        }
      } else {
        flag = 0;
        x1 = propose(x0);
        v1 = f(x1);
      }
      if (!isnan(v1)) {
        if (v1<v0) {
          swap(v0,v1);
          switch (flag) {
            case 0:
              alpha = beta = length(x0-x1);
              break;
            case 1:
              alpha=sqrt(alpha);
              break;
            case 2:
              beta*=2;
              break;
          }
          x0=x1;
          model.set_zero(x0,v0);
        } else {
          switch (flag) {
            case 0:
              model.strong_update(x1,v1);
              break;
            case 1:
              model.weak_update(x1,v1);
              if (alpha>1e-20) alpha*=exp(-1);
              break;
            case 2:
              model.weak_update(x1,v1);
              if (beta>1e-20) beta*=iphi;
              break;
          }
        }
      }
      if (stop(i,i,v0,v1)) break;
    }
  }
  
  template<long n,class target,class proposer,class stop_strategy>
  void inline find_min_numerical(
      matrix<double,n,1>& x0,
      const target& f,
      proposer& propose,
      const stop_strategy& stop = count_stop_strategy{4096,4096}
  )
  {
    find_min_numerical(x0,f,propose,null_projection<matrix<double,n,1>>{},stop);
  }
};

#endif // WMATH_OPTIMIZATION_H
