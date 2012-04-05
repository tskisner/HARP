

#ifndef _BOOST_UBLAS_LANCZOS_
#define _BOOST_UBLAS_LANCZOS_

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/exception.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/bindings/lapack.hpp>

#include <moat.hpp>

namespace boost { namespace numeric { namespace ublas {

  // This is based directly on the IETL code.  Will be replaced with
  // something more suitable once testcases work.

  template <class T>
  class FortranMatrix {
    private:
      T* p;
      std::size_t n_;
      std::size_t m_;
    public:
      typedef std::size_t size_type;
      FortranMatrix(size_type n, size_type m) : n_(n), m_(m) { p = new T[m*n]; };
      ~FortranMatrix() { delete[] p; };
      T* data() { return p; };
      const T* data() const { return p; };
      T operator()(size_type i, size_type j) const { return p[i+j*n_]; };
      T& operator()(size_type i, size_type j) { return p[i+j*n_]; };
      void resize(size_type n, size_type m) { m_=m; n_=n; delete[] p; p = new T[m*n]; };
      size_type nrows() { return n_; };
      size_type ncols() { return m_; };
      size_type minor() { return n_; };
  };

  template <class T> 
  struct number_traits {
    typedef T magnitude_type;
  };
    
  template <class T>
  struct number_traits<std::complex<T> > {
    typedef T magnitude_type;
  };

  template <class VS>
  struct vectorspace_traits {
    typedef typename VS::vector_type vector_type;
    typedef typename VS::size_type size_type;
    typedef typename VS::scalar_type scalar_type;
    typedef typename number_traits<scalar_type>::magnitude_type magnitude_type;
  };

  template<class V>
  class vectorspace {
    public:
      typedef V vector_type;
      typedef typename V::value_type scalar_type;
      typedef typename V::size_type size_type;
    
      vectorspace(size_type n):n_(n){}
    
      inline size_type vec_dimension() const {
        return n_;
      }
      vector_type new_vector() const {
        return vector_type(n_);
      }
    
      void project(vector_type&) const {
      }
    
    private:
      size_type n_;
  };
  
  template <class V, class S> class scaled_vector_wrapper;
  
  template <class V>
  class vector_wrapper : public boost::shared_ptr<V> {
    typedef boost::shared_ptr<V> super_type;
    public:
      vector_wrapper(V* p) : boost::shared_ptr<V>(p) {}
      operator V& () { return *super_type::get();}
      operator const V& () const { return *super_type::get();}
      const vector_wrapper operator += (const vector_wrapper& x) { *super_type::get() += *x.get(); return *this;}
      const vector_wrapper operator -= (const vector_wrapper& x) { *super_type::get() -= *x.get(); return *this;}
      template <class T> const vector_wrapper& operator *= (T x) { *super_type::get() *= x; return *this;}
      template <class T> const vector_wrapper& operator /= (T x) { *super_type::get() /= x; return *this;}
      template <class S>
      const vector_wrapper& operator += (const scaled_vector_wrapper<V,S>& x) { *super_type::get() += x.scalar()*x.vector(); return *this;}
      template <class S>
      const vector_wrapper& operator -= (const scaled_vector_wrapper<V,S>& x) { *super_type::get() -= x.scalar()*x.vector(); return *this;}
      template <class S>
      const vector_wrapper& operator = (const scaled_vector_wrapper<V,S>& x) { *super_type::get() = x.scalar()*x.vector(); return *this;}
  };
  
  template<class VS>
  void project(typename vectorspace_traits<VS>::vector_type& v, const VS& vs) {
    vs.project(v);
  }
  
  template<class V>
  class wrapper_vectorspace {
    public:
      typedef vector_wrapper<V> vector_type;
      typedef typename V::value_type scalar_type;
      typedef typename V::size_type size_type;
    
      wrapper_vectorspace(size_type n):n_(n){}
    
      inline size_type vec_dimension() const{
        return n_;
      }
      vector_type new_vector() const {
        return vector_wrapper<V>(new V(n_));
      }
    
      void project(vector_type& src) const {
      }  
    private:
      size_type n_;
  };
  
  template <class VS>
  typename vectorspace_traits<VS>::vector_type new_vector(const VS& vs) {
    return vs.new_vector();
  }

  template <class VS>
  typename vectorspace_traits<VS>::size_type vec_dimension(const VS& vs) {
    return vs.vec_dimension();
  } 

  template <class V, class VS>
  void project(vector_wrapper<V>& v, const VS& vs) {
    vs.project(v);
  }
 
  template <class V, class S>
  class scaled_vector_wrapper {
    public:
      scaled_vector_wrapper(const vector_wrapper<V>& v, S s) : v_(v), s_(s) {}
   
      const V& vector() const { return *v_.get();}
      S scalar() const { return s_;}
    private:
      const vector_wrapper<V>& v_;
      S s_;
  }; 

  template <class VS>
  class Tmatrix {
    public:
      typedef typename vectorspace_traits<VS>::scalar_type scalar_type; 
      typedef typename vectorspace_traits<VS>::vector_type vector_type;
      typedef typename vectorspace_traits<VS>::magnitude_type magnitude_type;
      typedef typename vectorspace_traits<VS>::size_type size_type;
    
      Tmatrix() : error_tol(pow(std::numeric_limits<magnitude_type>::epsilon(),2./3.)) {}
    
      inline void push_back(magnitude_type a, magnitude_type b);
      inline void push_back(std::pair<magnitude_type,magnitude_type> a_and_b) { push_back(a_and_b.first,a_and_b.second);}
    
      const std::vector<magnitude_type>& eigenvalues(bool discard_ghosts=true) const {
        if(!computed) compute();
        if (discard_ghosts)
          return eigval_distinct_noghost;
        else
          return eigval_distinct;
      }
    
      const std::vector<magnitude_type>& errors(bool discard_ghosts=true) const {
        if(!computed) compute();
        if (discard_ghosts)
          return err_noghost;
        else
          return err;
      }
    
      const std::vector<int>& multiplicities(bool discard_ghosts=true) const {
        if(!computed) compute();
        if (discard_ghosts)
          return multiplicty_noghost;
        else
          return multiplicty;
      }
    
    protected:
      std::vector<magnitude_type> alpha;
      std::vector<magnitude_type> beta;
      magnitude_type error_tol;
      mutable magnitude_type thold; 
    
    private:
      mutable bool computed;
      void compute() const;
      mutable magnitude_type multol; 
      mutable std::vector<magnitude_type> err; 
      mutable std::vector<magnitude_type> err_noghost;
      mutable std::vector<magnitude_type> eigval_distinct; // distinct eigen values.
      mutable std::vector<magnitude_type> eigval_distinct_noghost; // distinct eigen values.  
      mutable std::vector<int> multiplicty; 
      mutable std::vector<int> multiplicty_noghost;
      magnitude_type alpha_max;
      magnitude_type beta_max;
      magnitude_type beta_min;     
    }; 

  template <class VS>
  void Tmatrix<VS>::push_back(magnitude_type a, magnitude_type b) {
    computed = false;
    alpha.push_back(a);
    beta.push_back(b);
    if(alpha.size() == 1) {
      alpha_max = a;
      beta_min = beta_max = b;
    } else {
      if(a > alpha_max)
        alpha_max = a;
      if(b > beta_max) beta_max = b; 
      if(b < beta_min) beta_min = b;
    } 
  }
  
  template <class VS>
  void Tmatrix<VS>::compute() const {
    err.resize(0,0);
    eigval_distinct.resize(0,0);
    multiplicty.resize(0,0);
    
    err_noghost.resize(0,0);
    eigval_distinct_noghost.resize(0,0);
    multiplicty_noghost.resize(0,0);
    
    computed = true;
    int info,n;
    std::vector<magnitude_type> eval(alpha.size());
    std::vector<magnitude_type> etemp(alpha.size());

    // on return from stev function, eval contains the eigen values.
    n = alpha.size();

    eval = alpha;
    etemp = beta;
    FortranMatrix<magnitude_type> z2(n,n);

    info = boost::numeric::bindings::lapack::stev ( 'V', n, &(eval[0]), &(etemp[0]), z2.data() );

    if (info > 0)
      throw std::runtime_error("LAPACK error, stev function failed.");
    
    // tolerance values:
    multol = std::max(alpha_max,beta_max) * 2 * std::numeric_limits<magnitude_type>::epsilon() * (1000 + n); 
    thold = std::max(eval[0],eval[n-1]);
    thold = std::max(error_tol * thold, 5 * multol);
    
    // error estimates of eigen values starts:    
    // the unique eigen values selection, their multiplicities and corresponding errors calculation follows:
    
    magnitude_type temp = eval[0];
    eigval_distinct.push_back(eval[0]);
    int multiple = 1;
    
    for(int i = 1; i < n ; i++) {
      if((eval[i]- temp) > thold) {
        eigval_distinct.push_back(eval[i]);
        temp = eval[i];
        multiplicty.push_back(multiple);
        if(multiple > 1) {
          err.push_back(0.);
        } else {
          err.push_back(fabs(*beta.rbegin() * z2(n-1,i-1))); // *beta.rbegin() = betaMplusOne.
          multiple = 1;
        }
      } else {
        multiple++;
      }
    }
    
    // for last eigen value.
    multiplicty.push_back(multiple);
    if(multiple > 1) {
      err.push_back(0); 
    } else {
      err.push_back(fabs(*beta.rbegin() * z2(n-1,n-1))); // *beta.rbegin() = betaMplusOne.
    }

    // the unique eigen values selection, their multiplicities and corresponding errors calculation ends.
    
    // ghosts calculations starts:
    std::vector<magnitude_type> beta_g(alpha.size() - 1);
    std::vector<magnitude_type> alpha_g(alpha.size() - 1);
    
    std::copy(alpha.begin() + 1, alpha.end(), alpha_g.begin());
    std::copy(beta.begin() + 1, beta.end(), beta_g.begin());
    
    std::vector<magnitude_type> eval_g(alpha_g.size()); 

    std::vector<magnitude_type> etemp_g(alpha_g.size());
    eval_g = alpha_g;
    etemp_g = beta_g;

    info = boost::numeric::bindings::lapack::stev ( 'N', n-1, &(eval_g[0]), &(etemp_g[0]), NULL );

    if (info > 0)
      throw std::runtime_error("LAPACK error, stev function failed.");
    
    typename std::vector<magnitude_type>::iterator k;
    int i = 0, t2 = 0;
    for(k = eigval_distinct.begin(); k != eigval_distinct.end(); k++,i++) { 
      if(multiplicty[i] == 1) { // test of spuriousness for the eigenvalues whose multiplicity is one.
        for(int j = t2; j < n-1; j++,t2++) { // since size of reduced matrix is n-1
          if((eval_g[j] - *k) >= multol) break;
    
          if(fabs(*k - eval_g[j]) < multol) {
            multiplicty[i] = 0;
            err[i] = 0; // if eigen value is a ghost => error calculation not required, 0=> ignore error.
            t2++;
            break;
          }   
        }
      }
    } // end of outer for.
    
    i = 0;
    for(k = eigval_distinct.begin(); k != eigval_distinct.end(); k++,i++) {
      if(multiplicty[i] != 0) {
        eigval_distinct_noghost.push_back(*k);
        multiplicty_noghost.push_back(multiplicty[i]);
        err_noghost.push_back(err[i]);
      }
    }
  } // end of compute.

  // class Info starts, contains info about the eigen vector calculation.

  template < class magnitude_type=double >
    class Info {
    public:
      enum errorinfo {ok = 0, no_eigenvalue, not_calculated};      
      Info() {}    
      Info(std::vector<int> M1, std::vector<int> M2, std::vector<int> Ma,
        std::vector<magnitude_type> Eigenvalue,std::vector<magnitude_type> Residuum,
        std::vector<errorinfo> Status):
        m1_(M1),
        m2_(M2),
        ma_(Ma),
        eigenvalue_(Eigenvalue),
        residuum_(Residuum),
        status_(Status){}
    
      int m1(int i) const {return m1_[i];}
      int m2(int i) const {return m2_[i];}
      int ma(int i) const {return ma_[i];}
      int size() {return m1_.size();}
      magnitude_type eigenvalue(int i) const {return eigenvalue_[i];}
      magnitude_type residual(int i) const {return residuum_[i];}
      errorinfo error_info(int i) const {return status_[i];}    
    private:
      std::vector<int> m1_;
      std::vector<int> m2_;
      std::vector<int> ma_;
      std::vector<magnitude_type> eigenvalue_;
      std::vector<magnitude_type> residuum_;
      std::vector<errorinfo> status_;
  };


  // iteration classes

  template <class T>
  class basic_iteration {
    public:       
      basic_iteration(unsigned int max_iter, T reltol = 0., T abstol = 0.) : error(0), i(0), max_iter_(max_iter), rtol_(reltol), atol_(abstol) { }        
      bool finished(T r,T lambda) {
        if (converged(r,lambda))
          return true;
        else if (i < max_iter_)
          return false;
        else {
          fail(1,"maximum number of iterations exceeded");
          return true;
        }
      }
    
      inline bool converged(T r, T lambda) {
        return (r <= rtol_ * std::fabs(lambda) || r < atol_); // relative or absolute tolerance.
      }    
      
      inline void operator++() { ++i; }  
      inline bool first() { return i == 0; }    
      inline int error_code() { return error; }    
      inline unsigned int iterations() { return i; }    
      inline T relative_tolerance() { return rtol_; }
      inline T absolute_tolerance() { return atol_; }
      inline unsigned int max_iterations() { return max_iter_; }    
      inline void fail(int err_code) { error = err_code; }  
      inline void fail(int err_code, const std::string& msg) { error = err_code; err_msg = msg; }
    
    protected:
      int error;
      unsigned int i;    
      unsigned int max_iter_;
      T rtol_;
      T atol_;
      std::string err_msg;
  };
  
  template <class T, class Derived>
  class basic_lanczos_iteration {
    public:         
      basic_lanczos_iteration(unsigned int max_iter, T r = 0., T a = 0.) : error(0), i(0), 
      max_iter_(max_iter), rtol_(r), atol_(a) { }   
    
      template <class Tmatrix>
      bool finished(const Tmatrix& tmatrix) {
        if (static_cast<const Derived&>(*this).converged(tmatrix))
          return true;
        else if (i < max_iter_)
          return false;
        else {
          fail (1, "maximum number of iterations exceeded");
          return true;
        }
      }
    
      bool converged() const { return false;}
      void operator++() { ++i; }    
      bool first() const { return i == 0; }
      int error_code() const { return error; }  
      unsigned int iterations() const { return i; }  
      inline unsigned int max_iterations() { return max_iter_; }    
      T relative_tolerance() const { return rtol_; }
      T absolute_tolerance() const { return atol_; }   
      inline void fail(int err_code){ error = err_code; }  
      inline void fail(int err_code, const std::string& msg) { error = err_code; err_msg = msg; }
    
    protected:
      int error;
      unsigned int i;    
      unsigned int max_iter_;
      T rtol_;
      T atol_;
      std::string err_msg;
  };
  
  template <class T>
  class lanczos_iteration_nlowest : public basic_lanczos_iteration<T,lanczos_iteration_nlowest<T> > {
    typedef basic_lanczos_iteration<T,lanczos_iteration_nlowest<T> > super_type;
    public:         
      lanczos_iteration_nlowest(unsigned int max_iter, unsigned int n= 1, 
      T r = 100.*std::numeric_limits<T>::epsilon(), 
      T a = 100.*std::numeric_limits<T>::epsilon())
      : basic_lanczos_iteration<T,lanczos_iteration_nlowest<T> >(max_iter,r,a), n_(n) { }   
        
      template <class Tmatrix>
      bool converged(const Tmatrix& tmatrix) const { 
        if(super_type::iterations()>1) {
          const std::vector<T>& errs = tmatrix.errors();
          const std::vector<T>& vals = tmatrix.eigenvalues();      
          if(vals.size()<n_)
            return false;
          else { 
            for(unsigned int i = 0; i < n_; i++)
              if (errs[i] > std::max(super_type::absolute_tolerance(),super_type::relative_tolerance()*std::abs(vals[i])))
                return false;
              return true;
          }
        }
        return false;
      }    
    
    private:
      unsigned int n_; 
  };
  

  template <class T>
  class lanczos_iteration_nhighest : public basic_lanczos_iteration<T,lanczos_iteration_nhighest<T> > {
    typedef basic_lanczos_iteration<T,lanczos_iteration_nhighest<T> > super_type;
    public:     
    
      lanczos_iteration_nhighest(unsigned int max_iter, unsigned int n= 1,
        T r = 100.*std::numeric_limits<T>::epsilon(), 
        T a = 100.*std::numeric_limits<T>::epsilon())
        : basic_lanczos_iteration<T,lanczos_iteration_nhighest<T> >(max_iter,r,a), n_(n){}
    
      template <class Tmatrix>
      bool converged(const Tmatrix& tmatrix) const {
        if(super_type::iterations()>1) { 
          const std::vector<T>& errs = tmatrix.errors();
          const std::vector<T>& vals = tmatrix.eigenvalues();
  
          if(errs.size()<n_)
            return false;
          else { 
            for(int i = 0; i < n_; i++)
              if (errs[errs.size()-i - 1] > std::max(super_type::absolute_tolerance(),
                 super_type::relative_tolerance()*std::abs(vals[vals.size()-i-1])))
                return false;
              return true;
          }
        } 
        return false;
      }   
        
    private:
      unsigned int n_; 
  };
  
  
  template <class T>
  class fixed_lanczos_iteration : public basic_lanczos_iteration<T,fixed_lanczos_iteration<T> > {
    public:         
      fixed_lanczos_iteration(unsigned int max_iter)
      : basic_lanczos_iteration<T,fixed_lanczos_iteration<T> >(max_iter,0.,0.) { }   
        
      template <class Tmatrix>
      bool converged(const Tmatrix& ) const { return false;}    
  };


  template <class T>
  class bandlanczos_iteration_nlowest {
    public:
      bandlanczos_iteration_nlowest(unsigned int max_iter,T def_tol, 
        T dep_tol,T ghost_tol,
        bool ghost_discarding,unsigned int evs)
        : max_iter_(max_iter), def_tol_(def_tol),
        dep_tol_(dep_tol), ghost_tol_(ghost_tol),
        ghost_discarding_(ghost_discarding), evs_(evs) {
        i=0;
      };

      bool finished() const {
        if ( i < max_iter_ )
          return false;
        else
          return true;
      }
      inline void operator++() { ++i; };
      inline void operator--() { --i; };
      inline bool first() { return i == 0; };
      inline unsigned int iterations() { return i; };
      inline unsigned int evs() { return evs_; };
      inline unsigned int max_iter() { return max_iter_; };
      inline T def_tol() { return def_tol_; };
      inline T dep_tol() { return dep_tol_; };
      inline T ghost_tol() { return ghost_tol_; };
      inline bool ghost_discarding() { return ghost_discarding_; };
      inline bool low() { return true; };
    private:
      unsigned int i;
      unsigned int max_iter_;
      unsigned int evs_;
      T def_tol_;
      T dep_tol_;
      T ghost_tol_;
      bool ghost_discarding_;
  };

  template <class T>
  class bandlanczos_iteration_nhighest {
    public:
      bandlanczos_iteration_nhighest(unsigned int max_iter,T def_tol,
        T dep_tol,T ghost_tol,
        bool ghost_discarding, unsigned int evs)
        : max_iter_(max_iter), def_tol_(def_tol),
        dep_tol_(dep_tol), ghost_tol_(ghost_tol),
        ghost_discarding_(ghost_discarding), evs_(evs) {
        i=0;
      };
      bool finished() const {
        if ( i < max_iter_ )
          return false;
        else
          return true;
      }
      inline void operator++() { ++i; };
      inline void operator--() { --i; };
      inline bool first() { return i == 0; };
      inline unsigned int iterations() { return i; };
      inline unsigned int evs() { return evs_; };
      inline unsigned int max_iter() { return max_iter_; };
      inline T def_tol() { return def_tol_; };
      inline T dep_tol() { return dep_tol_; };
      inline T ghost_tol() { return ghost_tol_; };
      inline bool ghost_discarding() { return ghost_discarding_; };
      inline bool low() { return false; };
    private:
      unsigned int i;
      unsigned int max_iter_;
      unsigned int evs_;
      T def_tol_;
      T dep_tol_;
      T ghost_tol_;
      bool ghost_discarding_;
  };
  
  // class lanczos starts:
  template <class MATRIX, class VS>
  class lanczos : public Tmatrix<VS> { 
    typedef Tmatrix<VS> super_type;
    
    public:
      typedef typename vectorspace_traits<VS>::vector_type vector_type;
      typedef typename vectorspace_traits<VS>::scalar_type scalar_type;
      typedef typename vectorspace_traits<VS>::magnitude_type magnitude_type;
    
      lanczos(const MATRIX& matrix, const VS& vec); // constructor.
    
      template <class IT, class GEN>
      void calculate_eigenvalues(IT& iter, GEN gen) {generate_tmatrix(iter,gen);}
    
      template <class IT>
      void more_eigenvalues(IT& iter) { generate_tmatrix(iter); }
    
      template <class IN, class OUT, class GEN>
      void eigenvectors(IN in_eigvals_start, IN in_eigvals_end , OUT eig_vectors, Info<magnitude_type>& inf, GEN gen, int maxiter=0);    
  
    private:
      template <class IN> void find_m1m2(IN in_eigvals_start, IN in_eigvals_end);
      // m1 m2 finder for eigen vector calculation.
    
      template <class GEN> std::pair<magnitude_type,magnitude_type> make_first_step(GEN gen);
      std::pair<magnitude_type,magnitude_type> make_step(int j, vector_type& vec3);
  
      template <class IT, class GEN> void generate_tmatrix(IT& iter, GEN gen); 
      template <class IT> void generate_tmatrix(IT& iter);     
      // T matrix generator, Used in eigenvalues & eigenvectors calculation.

      const MATRIX& matrix_; 
      const VS vecspace_;
      vector_type startvector;
      vector_type vec2;
      unsigned int n; // index of vec2
      std::vector<int> M1, M2, Ma;
    
  }; // end of class lanczos.
  
  //-----------------------------------------------------------------------  
  // implementation of member functions start:  
  template <class MATRIX, class VS> //constructor:
  lanczos<MATRIX, VS>::lanczos(const MATRIX& matrix, const VS& vec):
    matrix_(matrix),
    vecspace_(vec),
    startvector(new_vector(vec)),
    vec2(new_vector(vec)),
    n(0) {}

    
  //-----------------------------------------------------------------------  
  // eigen vectors calculation starts:
  template <class MATRIX, class VS>
  template <class IN, class OUT, class GEN>
  void lanczos<MATRIX, VS>::eigenvectors(IN in_eigvals_start, IN in_eigvals_end , OUT eig_vectors, Info<magnitude_type>& inf, GEN gen_, int maxiter) {
    vector_type vec3 = new_vector(vecspace_); // a temporary vector.
    std::vector<vector_type> eigvectors; // contains ritz vectors.
    std::vector<std::vector<magnitude_type> > Tvectors; // contains
                                             // eigenvectors of T matrix.
    // calculation of eigen vectors of T matrix(consists of alphas & betas):    
    int n1 =  super_type::alpha.size();
    magnitude_type mamax, error, lambda;
    std::pair<magnitude_type,magnitude_type> a_and_b;
    unsigned int ma = 0, deltam;
    int nth, maMax = 0, count;
    std::vector<magnitude_type> eigenval_a, residuum;
    std::vector<typename Info<magnitude_type>::errorinfo> status;
    
    find_m1m2(in_eigvals_start, in_eigvals_end);    
    std::vector<int>::iterator M1_itr = M1.begin();
    std::vector<int>::iterator M2_itr = M2.begin();
    
    while(in_eigvals_start !=  in_eigvals_end) {
    
      int maxcount = 10;
      lambda = 0; count = 0;
      typename Info<magnitude_type>::errorinfo errInf = Info<magnitude_type>::ok;
      
      // calculation of ma starts:
      if(*M1_itr != 0 && *M2_itr != 0) {
        ma = (3 * (*M1_itr) + *M2_itr)/4 + 1;
        deltam = ((3 * (*M1_itr) + 5 * (*M2_itr))/8 + 1 - ma)/10 + 1;
      } else {         
        if(*M1_itr != 0 && *M2_itr == 0) {
          ma = (5 * (*M1_itr))/4 + 1;
          mamax = std::min((11 * n1)/8 + 12, (13 * (*M1_itr))/8 + 1);
          deltam = int((mamax - ma)/10) + 1;
          if (maxiter > 0)
            maxcount = maxiter/deltam;
        }      
        else {
          errInf = Info<magnitude_type>::no_eigenvalue;
          ma = 0;
        }
      } // calculation of ma ends.

      eigvectors.push_back(new_vector(vecspace_)); 
      // new ritz vector is being added in eigvectors.
      
      std::vector<magnitude_type> Tvector;        
      Tvectors.push_back(Tvector); 
      // new T matrix vector is being added in Tvectors.
            
      if (ma == 0)
        eigvectors.back()*=0.;
      
      if(ma != 0) {
        std::vector<magnitude_type> eval(ma);
        FortranMatrix<magnitude_type> z(ma,ma);

        std::vector<magnitude_type> etemp(ma);
    
        // on return, z contains all orthonormal eigen vectors of T matrix.
        do {
          if(ma > super_type::alpha.size()) { // size of T matrix is to be increased.
            fixed_lanczos_iteration<magnitude_type> iter_temp(ma);
            generate_tmatrix(iter_temp,gen_);
          }
          count++;

          eval = super_type::alpha;
          etemp = super_type::beta;

          int lapack_info = boost::numeric::bindings::lapack::stev ( 'V', ma, &(eval[0]), &(etemp[0]), z.data() );

          if (lapack_info > 0)
            throw std::runtime_error("LAPACK error, stev function failed.");
          
          // search for the value of nth starts, where nth is the nth eigen vector in z.            
          for(nth = ma -1; nth >= 0; nth--)
            if(fabs(eval[nth]- *in_eigvals_start) <= super_type::thold) 
              break;
          
          // search for the value of ith ends, where ith is the ith eigen vector in z.            
          if(nth == -1) {
            error = 0;  ma = 0;
            eigvectors.back()*=0.;
            errInf = Info<magnitude_type>::no_eigenvalue;
          } else { 
            error = fabs(super_type::beta[ma-1] * z(ma - 1, nth)); // beta[ma - 1] = betaMplusOne.
            if(error > super_type::error_tol) {
              ma += deltam; 
              eval.resize(ma);
              etemp.resize(ma);
              z.resize(ma,ma);
            }
          } // end of else
        } while(error > super_type::error_tol && count < maxcount);
        
        if(error > super_type::error_tol) {
          eigvectors.back()*=0.;
          errInf = Info<magnitude_type>::not_calculated;
        } else { // if error is small enough.
          if(ma != 0) {
            for(int i = 0; i < ma; i++)
              (Tvectors.back()).push_back(z(i,nth));
            if(ma > maMax) maMax = ma;
            lambda = eval[nth];
          } // end of if(ma != 0), inner.
        } // end of else{//if error is small enough.
      } // end of if(ma != 0).
      
      eigenval_a.push_back(lambda); // for Info object.
      Ma.push_back(ma); // for Info object.
      status.push_back(errInf);
      in_eigvals_start++; 
      M1_itr++;  
      M2_itr++; 
    } // end of while(in_eigvals_start !=  in_eigvals_end)
    
    // basis transformation of eigen vectors of T. These vectors are good 
    // approximation of eigen vectors of actual matrix.  
    typename std::vector<vector_type>::iterator eigenvectors_itr;
    typename std::vector<std::vector<magnitude_type> >::iterator Tvectors_itr;
    
    a_and_b = make_first_step(gen_);
    if(a_and_b.first != super_type::alpha[0]|| a_and_b.second != super_type::beta[0])
      throw std::runtime_error("T-matrix problem at first step");
      
    eigenvectors_itr = eigvectors.begin();
    Tvectors_itr = Tvectors.begin(); 
    
    while(eigenvectors_itr !=  eigvectors.end()) {
      if(!Tvectors_itr->empty()) {
        *eigenvectors_itr = (*Tvectors_itr)[0]*startvector;
        *eigenvectors_itr += (*Tvectors_itr)[1]*vec2; 
      }
      eigenvectors_itr++;
      Tvectors_itr++;
    }
    n=2;
    for(int j = 2; j < maMax; j++) {
      a_and_b = make_step(j-1,vec3);
      if(a_and_b.first != super_type::alpha[j-1]|| a_and_b.second != super_type::beta[j-1])
        throw std::runtime_error("T-matrix problem");
      
      ++n;
      eigenvectors_itr = eigvectors.begin();
      Tvectors_itr = Tvectors.begin();       
      while(eigenvectors_itr !=  eigvectors.end()){
        if(Tvectors_itr->size() > j)
          *eigenvectors_itr+=(*Tvectors_itr)[j]*vec2; 
        // vec2 is being added in one vector of eigvectors.
        eigenvectors_itr++;
        Tvectors_itr++;
      } // end of while loop.      
    } // end of for(int j = 2; j < maMax; j++).
    // end of basis transformation.  
    
    // copying to the output iterator & residuum calculation starts:    
    int i = 0;
    for(eigenvectors_itr = eigvectors.begin(); 
        eigenvectors_itr != eigvectors.end(); eigenvectors_itr++) {
      *eig_vectors = *eigenvectors_itr;
      eig_vectors++;
      moat::la::prod ( matrix_, *eigenvectors_itr, vec3 );
      vec3-=eigenval_a[i++]*(*eigenvectors_itr);
      
      // now vec3 is (A*v - eigenval_a*v); *eigenvectors_itr) is being added in vec3.
      residuum.push_back(norm_2(vec3));
    } // copying to the output iterator ends.    
    inf =  Info<magnitude_type>(M1, M2, Ma, eigenval_a,residuum, status);
  } // end of void eigenvector(....).
  
  //------------------------------------------------------  
  template <class MATRIX, class VS> template <class IN>
  void lanczos<MATRIX, VS>::find_m1m2(IN in_eigvals_start,  IN in_eigvals_end) {
    int info,  m2counter = 0;
    unsigned int n = 1;
    IN in_eigvals = in_eigvals_start;
    M1.resize(in_eigvals_end - in_eigvals_start,0);
    M2.resize(in_eigvals_end - in_eigvals_start,0);
    
    while(m2counter < (in_eigvals_end - in_eigvals_start)  && (n < super_type::alpha.size()) ) { 
      n++; // n++ == 2, at first time in this loop.
      std::vector<magnitude_type> eval(n);

      std::vector<magnitude_type> etemp(n);
      eval = super_type::alpha;
      etemp = super_type::beta;

      info = boost::numeric::bindings::lapack::stev ( 'N', n, &(eval[0]), &(etemp[0]), NULL );

      if (info > 0)
        throw std::runtime_error("LAPACK error, stev function failed.");
      
      std::vector<int>::iterator M1_itr = M1.begin();        
      std::vector<int>::iterator M2_itr = M2.begin();
      in_eigvals = in_eigvals_start;
      
      while(in_eigvals != in_eigvals_end) {        
        if(*M1_itr == 0 || *M2_itr == 0) {
          typename std::vector<magnitude_type> ::const_iterator lb, ub;
          ub = std::lower_bound(eval.begin(),eval.end(),*in_eigvals+super_type::thold);
          lb = std::upper_bound(eval.begin(),eval.end(),*in_eigvals-super_type::thold);
          if (*M1_itr == 0 && ub-lb >= 1)
            *M1_itr = n;          
          if (*M2_itr == 0 && ub-lb >= 2) {
            *M2_itr = n;
            m2counter++;
          }          
        } // end of "if(*M1_itr == 0 || *M2_itr ...".        
        in_eigvals++;
        M1_itr++;
        M2_itr++;
      } // end of inner while loop.
    } // end of outer while loop.
  } // end of function find_m1m2.

  
  //------------------------------------------------------  
  // generation of alpha, beta of T matrix starts:  
  template <class MATRIX, class VS> template <class IT, class GEN>
  void lanczos<MATRIX, VS>::generate_tmatrix(IT& iter, GEN gen_) {
    vector_type vec3 = new_vector(vecspace_);
    std::pair<magnitude_type,magnitude_type> a_and_b;
    
    if(super_type::alpha.size() == 0) {
      a_and_b = make_first_step(gen_);
      push_back(a_and_b); // member of T-matrix class.
      n=1;
    }    
    
    generate_tmatrix(iter);
  } // generation of alpha, beta of T matrix ends.
  
  template <class MATRIX, class VS> template <class IT>
  void lanczos<MATRIX, VS>::generate_tmatrix(IT& iter) {
    vector_type vec3 = new_vector(vecspace_);
    std::pair<magnitude_type,magnitude_type> a_and_b;
    
    for(int j = 0; j < n; j++) 
      ++iter;    
    if(super_type::alpha.size() == 0) 
      throw std::runtime_error("T matrix error, size of T matrix is zero, more_eigenvalues() cannot be called.");        
    
    while(!iter.finished(*this)) {          
      a_and_b = make_step(n,vec3);
      if (n==super_type::alpha.size()) 
        push_back(a_and_b); // member of T-matrix class 
      ++iter;
      ++n;
    }
  } // generation of alpha, beta of T matrix ends.


  //------------------------------------------------------  
  // generation of one step of alpha, beta:  
  template <class MATRIX, class VS> template <class GEN>
  std::pair<typename lanczos<MATRIX, VS>::magnitude_type,typename lanczos<MATRIX, VS>::magnitude_type> 
  lanczos<MATRIX, VS>::make_first_step(GEN gen) {
    magnitude_type a, b;
    generate(startvector,gen);      
    project(startvector,vecspace_);
    startvector/=norm_2(startvector); // normalization of startvector.
    moat::la::prod ( matrix_, startvector, vec2 );
    a = static_cast < magnitude_type > ( moat::la::dot ( startvector, vec2 ) );
    vec2 -= a*startvector;
    b = norm_2(vec2);   
    vec2 /= b;
    return std::make_pair(a,b);
  }  
    
  template <class MATRIX, class VS>
  std::pair<typename lanczos<MATRIX, VS>::magnitude_type, typename lanczos<MATRIX, VS>::magnitude_type> 
  lanczos<MATRIX, VS>::make_step(int j,vector_type& vec3) {
    magnitude_type a, b;
    b = super_type::beta[j-1];
    moat::la::prod ( matrix_, vec2, vec3 );
    a = static_cast < magnitude_type > ( moat::la::dot ( vec2, vec3 ) );
    vec3-=a*vec2;
    vec3-=b*startvector;
    b = norm_2(vec3);
    vec3/=b;
    std::swap(vec2,startvector); 
    std::swap(vec3,vec2);
    return std::make_pair(a,b);
  }




}}}

#endif
