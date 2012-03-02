#ifndef base_common_h
#define base_common_h

#include "cholmod.h"
#include <string>
#include <ostream>
#include <iostream>
#include "base_tools.h"


/* ==============================================
 * Class base_common
 *
 * This class in a wraper for cholmod_common.
 * Since cholmod_start(&c) must be the first call to solve
 * cholesky by CHOLMOD, then it is in the constructor.
 * Also, since cholmod_finish(&c) must be the last call, 
 * then it is in the destructor.
 * Note that base_common will be the first variable created
 * because the others will have to point to it.
 * ==============================================
 */


namespace base_matrices {

  class base_sparse;

  class base_common {
    public:
      base_common ();
      base_common (const base_common &); //Copy not allowed
      virtual ~base_common ();
      void operator= (const base_common &); //Copy not allowed

      //CHOLMOD Check functions
      int check () const;
      int print (std::string) const;
      int get_status () const { return common->status; };
      //End

      //CHOLMOD Core functions
      int defaults ();
      size_t maxrank (size_t);
      int allocate_work (size_t, size_t, size_t);
      int free_work ();
      UF_long clear_flag ();
      int error (int, const std::string, int, const std::string) const;
      //End
      
      //Other functions
      void set_error_handler (void (*eh)(int, const char*, int, const char*));
      void useSimplicial () {
        common->supernodal = CHOLMOD_SIMPLICIAL;
      }
      void useSupernodal () {
        common->supernodal = CHOLMOD_SUPERNODAL;
      }

      friend class base_sparse;
      friend class base_dense;
      friend class base_factor;
      friend class base_triplet;
      friend base_sparse add (const base_sparse &, const base_sparse &, double[2], double[2], int, int);
    protected:
      cholmod_common *common;
  };
}

#endif
