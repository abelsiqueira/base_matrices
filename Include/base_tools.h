#include <iostream>

namespace base_matrices {

  template < typename T >
  T max (T a, T b) {
    if (a > b)
      return a;
    else
      return b;
  }


#ifdef LONG
  typedef long int bmInt;
#define CHOLMOD(name) cholmod_l_ ## name
#else
  typedef int bmInt;
#define CHOLMOD(name) cholmod_ ## name
#endif

}
