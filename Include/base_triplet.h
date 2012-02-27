#ifndef base_triplet_h
#define base_triplet_h

#include "base_dense.h"

/* =================================
 * class base_triplet
 * 
 * We define the class for triplet matrices
 * When an object of this class is created, it must be associated with
 * an ***existing*** base_common object. Hence, we declare a default
 * constructor, but do not define him, preventing calls to it to be
 * compiled.
 * ==================================
 */

namespace base_matrices {
  class base_triplet {
    public:
      base_triplet (); //Not defined, because not allowed!
      base_triplet (const base_common &, size_t = 0, size_t = 0, size_t = 0, int = 0, int = CHOLMOD_REAL);
      base_triplet (const base_triplet &);
      void operator= (const base_triplet &);
      virtual ~base_triplet ();

      //CHOLMOD Check functions
      int check () const;
      int print (std::string) const;
      void read (FILE *);
      //End

      //CHOLMOD Core functions
      int reallocate (size_t);
      void sparse_to_triplet (const base_sparse &); 
      int xtype_change (int);
      //End

      //other functions
      void print_more (std::ostream & = std::cout) const;
      void error (int, std::string) const;
      size_t get_nrow () const { return triplet->nrow; }
      size_t get_ncol () const { return triplet->ncol; }
      size_t get_nnz () const { return triplet->nnz; }
      size_t get_nzmax () const { return triplet->nzmax; }

      //access to common and base_common &
      const base_common & get_common () const {
        return *common;
      }

      friend class base_sparse;
      friend class base_dense;
    protected:
      cholmod_triplet *triplet;
      const base_common *common;
      cholmod_common * get_cholmod_common () const {
        return common->common;
      }
      cholmod_triplet * get_triplet () const {
        return triplet;
      }
  };
}

#endif
