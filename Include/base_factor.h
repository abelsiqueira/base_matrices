#ifndef base_factor_h
#define base_factor_h

#include "base_sparse.h"

/* =================================
 * class base_factor
 * 
 * We define the class for the cholesky factor
 * When an object of this class is created, it must be associated with
 * an ***existing*** base_common object. Hence, we declare a default
 * constructor, but do not define him, preventing calls to it to be
 * compiled.
 * ==================================
 */

namespace base_matrices {
  class base_factor {
    public:
      base_factor (); //Not defined, because not allowed!
      base_factor (const base_common &, size_t = 0);
      base_factor (const base_factor &);
      void operator= (const base_factor &);
      virtual ~base_factor ();

      //CHOLMOD Check functions
      int check () const;
      int print (std::string) const;
      //End

      //CHOLMOD Cholesky functions
      void analyze (const base_sparse &);
      void factorize (const base_sparse &, double = 0);
      //End

      //CHOLMOD Core functions
      int reallocate (size_t);
      int change (int, int, int, int, int);
      int pack ();
      int reallocate_column (size_t, size_t);
      int change_xtype (int);
      //End

      //Other functions
      void error (int, std::string) const;

      //access to base_common &
      const base_common & get_common () const {
        return *common;
      }

      friend class base_dense;
      friend class base_sparse;
    protected:
      cholmod_factor *factor;
      const base_common *common;
      cholmod_common * get_cholmod_common () const {
        return common->common;
      }
      cholmod_factor * get_factor () const {
        return factor;
      }
  };
}

#endif
