namespace wmath{
  using std::allocator_traits;
  using std::array;
  using std::cerr;
  using std::conditional;
  using std::cout;
  using std::enable_if;
  using std::endl;
  using std::false_type;
  using std::get;
  using std::index_sequence;
  using std::index_sequence_for;
  using std::initializer_list;
  using std::integral_constant;
  using std::is_const;
  using std::is_fundamental;
  using std::is_same;
  using std::is_signed;
  using std::is_trivially_copyable;
  using std::make_unique;
  using std::numeric_limits;
  using std::pair;
  using std::setw;
  using std::string;
  using std::swap;
  using std::true_type;
  using std::tuple;
  using std::tuple_size;
  using std::unique_ptr;
  using std::remove_const;

  template<typename test, template<typename...> class ref>
  struct is_specialization : std::false_type {};

  template<template<typename...> class ref, typename... args>
  struct is_specialization<ref<args...>, ref>: std::true_type {};

  template<typename T>
  typename std::enable_if<std::is_unsigned<T>::value,T>::type
  constexpr bijective_square(const T& n) {
    T i = n;
    T j = n+1;
    if (j&1) i>>=1;
    if (n&1) j>>=1;
    return i*j;
  }

  /* now found in "wmath_bits.hpp" as modular_inverse_power2 T
  template<typename T>
  typename std::enable_if<std::is_unsigned<T>::value,T>::type
  constexpr modular_inverse(const T& a) {
    T x(1u);
    for (size_t i(1);i!=digits<T>();++i) x*=T(2u)-a*x;
    return x;
  }
  */

  template<typename T,class enable = void>
  class hash;
  
  template<>
  class hash<uint8_t>{
    private:
      const uint8_t a = 97u;
      const uint8_t i = modular_inverse_power2(a);
      const uint8_t b = 111u;
    public:
      typedef typename true_type::type is_injective;
      typedef typename true_type::type unhash_defined;
      constexpr size_t digits() const {
        return wmath::digits<uint8_t>();
      }
      constexpr uint8_t operator()(const uint8_t v) const {
        return (v+b)*a;
      }
      constexpr uint8_t unhash(const uint8_t v) const {
        return v*i-b;
      }
  };
  
  template<>
  class hash<uint16_t>{
    private:
      const uint16_t a = 43581u;
      const uint16_t i = modular_inverse_power2(a);
      const uint16_t b = 36690u;
    public:
      typedef typename true_type::type is_injective;
      typedef typename true_type::type unhash_defined;
      constexpr size_t digits() const {
        return wmath::digits<uint16_t>();
      }
      constexpr uint16_t operator()(const uint16_t v) const {
        return (v+b)*a;
      }
      constexpr uint16_t unhash(const uint16_t v) const {
        return v*i-b;
      }
  };

  template<>
  class hash<uint32_t>{
    private:
      const uint32_t a = 3370923577ul;
      const uint32_t i = modular_inverse_power2(a);
    public:
      typedef typename true_type::type is_injective;
      typedef typename true_type::type unhash_defined;
      constexpr size_t digits() const {
        return wmath::digits<uint32_t>();
      }
      constexpr uint32_t operator()(uint32_t v) const {
        v^= v>>16;
        v*= a;
        v^= v>>16;
        return v;
      }
      constexpr uint32_t unhash(uint32_t v) const {
        v^= v>>16;
        v*= i;
        v^= v>>16;
        return v;
      }
  };

  template<>
  class hash<uint64_t>{
    private:
      const uint64_t  a = 15864664792644967873ull;
      const uint64_t  i = modular_inverse_power2(a);
    public:
      typedef typename true_type::type is_injective;
      typedef typename true_type::type unhash_defined;
      constexpr size_t digits() const {
        return wmath::digits<uint64_t>();
      }
      uint64_t constexpr operator()(uint64_t v) const {
        v^= v>>32;
        v*= a;
        v^= v>>32;
        return v;
      }  
      uint64_t constexpr unhash(uint64_t v) const {
        v^= v>>32;
        v*= i;
        v^= v>>32;
        return v;
      }
  };

  template<typename S>
  class hash<S,
        typename enable_if<is_integral<S>::value&&is_signed<S>::value>::type>
  {
    private:
      using U = typename make_unsigned<S>::type;
      const hash<typename remove_const<U>::type> hasher{};
    public:
      typedef typename true_type::type is_injective;
      typedef typename true_type::type unhash_defined;
      constexpr size_t digits() const {
        return hasher.digits();
      }
      constexpr U operator()(const S& v) const {
        return hasher(U(v));
      }
      constexpr S unhash(const U& v) const {
        return S(hasher.unhash(v));
      }
  };
  
  template<typename T>
  class hash<T,typename enable_if<
  (is_fundamental<T>::value)&&(!is_integral<T>::value)
  >::type>
  {
    public:
      typedef typename false_type::type is_injective;
      typedef typename false_type::type unhash_defined;
      constexpr size_t digits() const {
        return wmath::digits<size_t>();
      }
      const inline size_t operator()(const T& v) const {
        size_t h = 0;
        for (size_t i=0;i!=sizeof(T);++i) {
          h^=*(reinterpret_cast<const char*>(&v)+i);
          rol(h,wmath::digits<char>());
          if (((i%sizeof(size_t))==(sizeof(size_t)-1))||(i==(sizeof(T)-1))) {
            h = hash<size_t>{}(h);
          }
          //cerr << "# " << h << endl;
        }
        return h;
      }
  };


  template<>
  class hash<string>
  {
    public:
      typedef typename false_type::type is_injective;
      typedef typename false_type::type unhash_defined;
      constexpr size_t digits() {
        return wmath::digits<size_t>();
      }
      const inline size_t operator()(const string& s) const {
        return std::hash<string>{}(s);
      }
  };

  template<typename T>
  class hash<
    T,
    typename enable_if<
      is_specialization<
        typename remove_const<T>::type,
        tuple
      >::value
    >::type
  >
  {
    private:
      template<size_t i = 0>
      const size_t impl(const T& v,const size_t& h=0u)
      const
      {
        if constexpr (i==tuple_size<T>::value) {
          return h;
        } else {
          const auto e = get<i>(v);
          const size_t t = hash<decltype(e)>{}(e);
          return impl<i+1>(v,h^(size_t(2u*i+1u)*t));
        }
      }
    public:
      typedef typename false_type::type is_injective;
      typedef typename false_type::type unhash_defined;
      constexpr size_t digits() const {
        return wmath::digits<size_t>();
      }
      const inline size_t operator()(const T& v) const {
        return impl(v);
      }
  };
  
  template<typename T>
  class hash<
    T,
    typename enable_if<
      is_specialization<typename remove_const<T>::type,vector>::value
    >::type
  >
  {
    public:
      typedef typename false_type::type is_injective;
      typedef typename false_type::type unhash_defined;
      constexpr size_t digits() const {
        return wmath::digits<size_t>();
      }
      const inline size_t operator()(const T& v) const {
        size_t h = 0;
        for (size_t i=0;i!=v.size();++i) {
          auto e = v[i];
          h^= size_t(2u*i+1u)*hash<decltype(e)>{}(e);
        }
        return h;
      }
  };
}
