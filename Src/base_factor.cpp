/* Copyright 2012 Abel Soares Siqueira
 *
 * base_matrices 
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
#include "base_factor.h"

namespace base_matrices {

  base_factor::base_factor (const base_common & cm, size_t n) {
    common = &cm;
    if (n == 0)
      factor = 0;
    else
      factor = CHOLMOD(allocate_factor) (n, common->common);
  }

  base_factor::base_factor (const base_factor & A) {
    common = A.common;
    if (A.factor == 0) {
      error (16, "ERROR: copying uninitialized base_factor");
      factor = 0;
    } else {
      factor = CHOLMOD(copy_factor) (A.factor, common->common);
    }
  }
  
  void base_factor::operator= (const base_factor & A) {
    if (A.factor == 0) {
      error (25, "ERROR: cannot assign to uninitialized base_factor");
      return;
    }
    if (A.factor == factor)
      return;
    if (factor != 0)
      CHOLMOD(free_factor) (&factor, common->common);
    common = A.common;
    factor = CHOLMOD(copy_factor) (A.factor, common->common);
  }

  base_factor::~base_factor () {
    CHOLMOD(free_factor) (&factor, common->common);
  }

  //CHOLMOD Check functions
  int base_factor::check () const {
    return CHOLMOD(check_factor) (factor, common->common);
  }

  int base_factor::print (std::string name) const {
    return CHOLMOD(print_factor) (factor, name.c_str(), common->common);
  }
  //End

  //CHOLMOD Cholesky functions
  void base_factor::analyze (const base_sparse & A) {
    if (A.sparse == 0) {
      error (51, "ERROR: cannot analyze uninitialized base_sparse");
      return;
    }
    if (factor != 0)
      CHOLMOD(free_factor) (&factor, common->common);
    factor = CHOLMOD(analyze) (A.sparse, common->common);
  }

  void base_factor::factorize (const base_sparse & A, double beta) {
    if (A.sparse == 0) {
      error (59, "ERROR: cannot factorize uninitialized base_sparse");
      return;
    } else if (factor == 0) {
      error (62, "ERROR: cannot factorize. base_sparse must be analyzed first");
      return;
    }
    double beta2[2] = {beta, 0};
    if (beta == 0)
      CHOLMOD(factorize) (A.sparse, factor, common->common);
    else
      CHOLMOD(factorize_p) (A.sparse, beta2, 0, 0, factor, common->common);
  }

  void base_sparse::spsolve (int sys, const base_factor & L, const base_sparse & B) {
    if (L.common != B.common) {
      L.error (70, "ERROR: base_factor and base_sparse have different base_common");
      return;
    } else if (L.factor == 0) {
      L.error (73, "ERROR: base_factor uninitialized");
      return;
    } else if (B.sparse == 0) {
      L.error (76, "ERROR: base_sparse uninitialized");
      return;
    }
    if (sparse != 0)
      CHOLMOD(free_sparse) (&sparse, common->common);
    common = L.common;
    sparse = CHOLMOD(spsolve) (sys, L.factor, B.sparse, common->common);
  }
  //End
  
  //CHOLMOD Core functions
  int base_factor::reallocate (size_t nznew) {
    return CHOLMOD(reallocate_factor) (nznew, factor, common->common);
  }

  int base_factor::change (int to_xtype, int to_ll, int to_super, int to_packed, 
      int to_monotonic) {
    return CHOLMOD(change_factor) (to_xtype, to_ll, to_super, to_packed, to_monotonic, 
        factor, common->common);
  }

  int base_factor::pack () {
    return CHOLMOD(pack_factor) (factor, common->common);
  }

  int base_factor::reallocate_column (size_t j, size_t need) {
    return CHOLMOD(reallocate_column) (j, need, factor, common->common);
  }

  void base_sparse::factor_to_sparse (const base_factor & L) {
    if (L.factor == 0) {
      L.error (107, "ERROR: base_factor uninitialized");
      return;
    }
    if (sparse != 0)
      CHOLMOD(free_sparse) (&sparse, common->common);

    common = L.common;
    sparse = CHOLMOD(factor_to_sparse) (L.factor, common->common);
  }

  int base_factor::change_xtype (int to_xtype) {
    return CHOLMOD(factor_xtype) (to_xtype, factor, common->common);
  }
  //End

  //Other functions
  void base_factor::error (int line, std::string message) const {
    common->error (-10, "base_factor.cpp", line, message);
  }

}
