#ifndef base_dense_h
#define base_dense_h

#include "base_factor.h"

/* =================================
 * class base_dense
 * 
 * We define the class for dense matrices
 * When an object of this class is created, it must be associated with
 * an ***existing*** base_common object. Hence, we declare a default
 * constructor, but do not define him, preventing calls to it to be
 * compiled.
 * ==================================
 */

namespace base_matrices {
  class base_dense {
    public:
      base_dense (); //Not defined, because not allowed!
      base_dense (const base_common &, size_t = 0, size_t = 0, size_t = 0, int = CHOLMOD_REAL);
      base_dense (const base_common &, std::string, size_t, size_t, int = CHOLMOD_REAL);
      base_dense (const base_dense &);
      void operator= (const base_dense &);
      virtual ~base_dense ();

      //CHOLMOD Check functions
      int check () const;
      int print (std::string) const;
      void read (FILE *);
      int write (FILE *) const;
      int write (FILE *, std::string) const;
      //End

      //CHOLMOD Cholesky functions
      void solve (int, const base_factor &, const base_dense &);
      //End

      //CHOLMOD Core functions
      void sparse_to_dense (const base_sparse &); 
      int xtype_change (int);
      //End

      //CHOLMOD MatrixOps functions
      double norm (int = 2) const;
      int sdmult (const base_sparse &, int, double [2], double [2], const base_dense &);
      //End

      //other functions
      size_t get_nrow () const { return dense->nrow; }
      size_t get_ncol () const { return dense->ncol; }
      double dot (const base_dense &) const;
      void print_more (std::ostream & = std::cout) const;
      void saxpy (const base_dense &, double);
      void scale (const base_dense &, double);
      void scale (double);

      double get (size_t, size_t) const;
      void set (size_t, size_t, double);

      base_dense operator+ (const base_dense &) const;
      base_dense operator- (const base_dense &) const;
      base_dense operator* (double) const;
      friend base_dense operator* (double, const base_dense &);
      void error (int, std::string) const;

      //access to base_common &
      const base_common & get_common () const {
        return *common;
      }

      friend class base_sparse;
    protected:
      cholmod_dense *dense;
      const base_common *common;
      cholmod_common * get_cholmod_common () const {
        return common->common;
      }
      cholmod_dense * get_dense () const {
        return dense;
      }
  };
}

#endif
