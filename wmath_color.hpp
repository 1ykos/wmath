namespace wmath{
  uint64_t constexpr apply_srgb_gamma(const uint64_t& v) {
    if (v<=57753066345969864ull) {
      return 323*v/25; // 12.92 * v
    } else {
      // TODO integer nth root and replace 2.4 by 12/5
      return uint64_t(1.8446744073709552e+19*
          (1.055*pow(v*5.421010862427522e-20,2.4)))
        -1014570924054025338ull;
    }
  }
  uint64_t constexpr revert_srgb_gamma(const uint64_t& v) {
    if (v<=746170797781551360ull) {
      return 25*v/323;
    } else {
      // TODO integer nth root and replace 1/2.4 by 5/12 
      return uint64_t(1.8446744073709552e+19*
      pow(5.421010862427522e-20*(v+1014570924054025338ull)/1.055,1.0/2.4));
    }
  }
  const inline void hsv_rgb(
      uint64_t& h,
      uint64_t& s,
      uint64_t& v
      ){
    const uint64_t c = 3074457345618258602ull;
    const uint64_t o = ~uint64_t(0);
    const uint64_t i = h/c;
    const uint64_t f = 6*(h%c);
    //cerr << i << " " << f/double(o) << " "
    //     << (h%c)/double(o) << endl;
    const uint64_t p = get<0>(long_mul(v,o-s));
    const uint64_t q = get<0>(long_mul(v,o-get<0>(long_mul(s,f))));
    const uint64_t t = get<0>(long_mul(v,o-get<0>(long_mul(s,o-f))));
    switch (i){
      case 0:
      case 6:
        h = v;
        s = t;
        v = p;
        return;
      case 1:
        h = q;
        s = v;
        v = p;
        return;
      case 2:
        h = p;
        s = v;
        v = t;
        return;
      case 3:
        h = p;
        s = q;
        return;
      case 4:
        h = t;
        s = p;
        return;
      case 5:
        h = v;
        s = p;
        v = q;
        return;
    };
  }
}
