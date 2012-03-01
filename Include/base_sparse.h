#ifndef base_sparse_h
#define base_sparse_h

#include "base_common.h"

/* =================================
 * class base_sparse
 * 
 * We define the class for sparse matrices.
 * When an object of this class is created, it must be associated with
 * an ***existing*** base_common object. Hence, we declare a default
 * constructor, but do not define him, preventing calls to it to be
 * compiled.
 * ==================================
 */


namespace base_matrices {
  class base_dense;
  class base_factor;
  class base_triplet;
  class base_sparse {
    public:
      base_sparse (); //Not defined, because not allowed!
      base_sparse (const base_common &, size_t = 0, size_t = 0, size_t = 0, int = 1, int = 1, int = 0, int = CHOLMOD_REAL);
      base_sparse (const base_common &, std::string, size_t, size_t, size_t, int = CHOLMOD_REAL);
      base_sparse (const base_sparse &);
      void operator= (const base_sparse &);
      virtual ~base_sparse ();

      //CHOLMOD Check functions
      int check () const;
      int print (std::string) const;
      void read (FILE *);
      int write (FILE *) const;
      int write (FILE *, const base_sparse &) const;
      int write (FILE *, const std::string) const;
      int write (FILE *, const base_sparse &, std::string) const;
      //End

      //CHOLMOD Cholesky functions
      void spsolve (int, const base_factor &, const base_sparse &);
      //End

      //CHOLMOD Core functions
      int reallocate (size_t);
      UF_long nnz () const;
      void transpose (const base_sparse &, int = 2);
      int sort ();
      void band (const base_sparse &, UF_long, UF_long, int);
      int band_inplace (UF_long, UF_long, int);
      void aat (const base_sparse &, bmInt * = 0, size_t = 0, int = 1);
      void copy (const base_sparse &, int, int);
      void add (const base_sparse &, const base_sparse &, double = 1, double = 1, int = 1, int = 0);
      void add (const base_sparse &, const base_sparse &, double[2], double[2], int = 1, int = 0);
      int xtype_change (int);
      //From Dense:
      void dense_to_sparse (const base_dense &, int = 1);
      //From Factor:
      void factor_to_sparse (const base_factor &);
      //From Triplet:
      void triplet_to_sparse (const base_triplet &, size_t);
      //End

      //CHOLMOD MatrixOps functions
      int drop (double);
      double norm (int = 0) const;
      void horzcat (const base_sparse &, const base_sparse &, int = 1);
      int scale (base_dense &, int);
      void ssmult (const base_sparse &, const base_sparse &, int = 0, int = 1, int = 0);
      void submatrix (const base_sparse &, bmInt *, UF_long, bmInt *, UF_long, int = 1, int = 0);
      void vertcat (const base_sparse &, const base_sparse &, int = 1);
      int symmetry (int, bmInt *, bmInt *, bmInt *, bmInt *);
      //End

      //Other functions
      void print_more (std::ostream & = std::cout) const;
      size_t get_nrow () const { return sparse->nrow; }
      size_t get_ncol () const { return sparse->ncol; }
      void error (int, std::string) const;

      //access to base_common &
      const base_common & get_common () const {
        return *common;
      }

      friend class base_dense;
      friend class base_factor;
    protected:
      cholmod_sparse *sparse;
      const base_common *common;
      cholmod_common * get_cholmod_common () const {
        return common->common;
      }
      cholmod_sparse * get_sparse () const {
        return sparse;
      }
  };
}

#endif
