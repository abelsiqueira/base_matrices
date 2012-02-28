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
/* Copyright 2012 Abel Soares Siqueira
 *
 * base_matrices - version 0.6.2.1
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
/* Copyright 2012 Abel Soares Siqueira
 *
 * Turn on the lights - version 0.8
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
#include "base_triplet.h"

namespace base_matrices {

  base_triplet::base_triplet (const base_common & cm, size_t nrow, size_t ncol, size_t nzmax, int stype, int xtype) {
    common = &cm;
    if ( (nrow == 0) || (ncol == 0) ) {
      triplet = 0;
    } else 
      triplet = CHOLMOD(allocate_triplet) (nrow, ncol, nzmax, stype, xtype, common->common);
  }

  base_triplet::base_triplet (const base_triplet & A) {
    common = A.common;
    if (A.triplet == 0) {
      error (31, "ERROR: copying uninitialized base_triplet");
      triplet = 0;
    } else {
      triplet = CHOLMOD(copy_triplet) (A.triplet, common->common);
    }
  }

  void base_triplet::operator= (const base_triplet & A) {
    if (A.triplet == 0) {
      error (40, "ERROR: cannot assign to uninitialized base_triplet");
      return;
    }
    if (A.triplet == triplet)
      return;
    if (triplet != 0)
      CHOLMOD(free_triplet) (&triplet, common->common);
    common = A.common;
    triplet = CHOLMOD(copy_triplet) (A.triplet, common->common);
  }

  base_triplet::~base_triplet () {
    CHOLMOD(free_triplet) (&triplet, common->common);
  }

  //CHOLMOD Check functions
  int base_triplet::check () const {
    return CHOLMOD(check_triplet) (triplet, common->common);
  }

  int base_triplet::print (std::string name) const {
    return CHOLMOD(print_triplet) (triplet, name.c_str (), common->common);
  }

  void base_triplet::read (FILE *f) {
    if (f == 0) {
      error (64, "ERROR: file is not open");
      return;
    }
    triplet = CHOLMOD(read_triplet) (f, common->common);
  }

  //End

  //CHOLMOD Core functions
  int base_triplet::reallocate (size_t nznew) {
    return CHOLMOD(reallocate_triplet) (nznew, triplet, common->common);
  }

  int base_triplet::xtype_change (int to_xtype) {
    return CHOLMOD(triplet_xtype) (to_xtype, triplet, common->common);
  }

  void base_sparse::triplet_to_sparse (const base_triplet & L, size_t nzmax) {
    if (L.triplet == 0) {
      L.error (107, "ERROR: base_triplet uninitialized");
      return;
    }
    if (sparse != 0)
      CHOLMOD(free_sparse) (&sparse, common->common);

    common = L.common;
    sparse = CHOLMOD(triplet_to_sparse) (L.triplet, nzmax, common->common);
  }
  //End

  //Other functions
  void base_triplet::print_more (std::ostream & outstr) const {
    if (triplet == 0) {
      outstr << "base_triplet not initialized\n";
      return;
    }
    double * px = static_cast < double * > (triplet->x);
    UF_long * pi = static_cast < UF_long * > (triplet->i);
    UF_long * pj = static_cast < UF_long * > (triplet->j);
    size_t nnz = triplet->nnz;

    for (size_t i = 0; i < nnz; i++) {
      outstr << "(" << pi[i]
             << "," << pj[i]
             << ") = " << px[i]
             << "\n";
    }

  }

  void base_triplet::error (int line, std::string message) const {
    common->error (-10, "base_triplet.cpp", line, message);
  }

}
