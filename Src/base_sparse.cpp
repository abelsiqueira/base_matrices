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
#include "base_sparse.h"

namespace base_matrices {

  base_sparse::base_sparse (const base_common & cm, size_t nrow, size_t ncol, size_t nzmax, int sorted, int packed, int stype, int xtype) {
    common = &cm;
    if ( (nrow == 0) || (ncol == 0) ) {
      sparse = 0;
    } else 
      sparse = CHOLMOD(allocate_sparse) (nrow, ncol, nzmax, sorted, packed, stype, xtype, common->common);
  }

  base_sparse::base_sparse (const base_common & cm, std::string type, size_t nrow, size_t ncol, size_t nzmax, int xtype) {
    common = &cm;
    if ( (nrow == 0) || (ncol == 0) ) {
      sparse = 0;
    } else if (type == "speye") {
      sparse = CHOLMOD(speye) (nrow, ncol, xtype, common->common);
    } else if (type == "spzeros") {
      sparse = CHOLMOD(spzeros) (nrow, ncol, nzmax, xtype, common->common);
    } else {
      sparse = 0;
    }
  }

  base_sparse::base_sparse (const base_sparse & A) {
    common = A.common;
    if (A.sparse == 0) {
      error (31, "ERROR: copying uninitialized base_sparse");
      sparse = 0;
    } else {
      sparse = CHOLMOD(copy_sparse) (A.sparse, common->common);
    }
  }

  void base_sparse::operator= (const base_sparse & A) {
    if (A.sparse == 0) {
      error (40, "ERROR: cannot assign to uninitialized base_sparse");
      return;
    }
    if (A.sparse == sparse)
      return;
    if (sparse != 0)
      CHOLMOD(free_sparse) (&sparse, common->common);
    common = A.common;
    sparse = CHOLMOD(copy_sparse) (A.sparse, common->common);
  }

  base_sparse::~base_sparse () {
    CHOLMOD(free_sparse) (&sparse, common->common);
  }

  //CHOLMOD Check functions
  int base_sparse::check () const {
    return CHOLMOD(check_sparse) (sparse, common->common);
  }

  int base_sparse::print (std::string name) const {
    return CHOLMOD(print_sparse) (sparse, name.c_str (), common->common);
  }

  void base_sparse::read (FILE *f) {
    if (f == 0) {
      error (64, "ERROR: file is not open");
      return;
    }
    sparse = CHOLMOD(read_sparse) (f, common->common);
  }

  int base_sparse::write (FILE *f) const {
    if (f == 0) {
      error (72, "ERROR: file is not open");
      return 0;
    } else if (sparse == 0) {
      error (75, "ERROR: cannot write uninitialized base_sparse");
      return 0;
    }
    return CHOLMOD(write_sparse) (f, sparse, 0, 0, common->common);
  }

  int base_sparse::write (FILE *f, const base_sparse & Z) const {
    if (f == 0) {
      error (83, "ERROR: file is not open");
      return 0;
    } else if (sparse == 0) {
      error (86, "ERROR: cannot write uninitialized base_sparse");
      return 0;
    } else if (Z.sparse == 0) {
      error (89, "ERROR: pattern base_sparse is uninitialized");
      return 0;
    }
    return CHOLMOD(write_sparse) (f, sparse, Z.sparse, 0, common->common);
  }
  
  int base_sparse::write (FILE *f, std::string comments) const {
    if (f == 0) {
      error (97, "ERROR: file is not open");
      return 0;
    } else if (sparse == 0) {
      error (100, "ERROR: cannot write uninitialized base_sparse");
      return 0;
    }
    return CHOLMOD(write_sparse) (f, sparse, 0, comments.c_str(), common->common);
  }

  int base_sparse::write (FILE *f, const base_sparse & Z, std::string comments) const {
    if (f == 0) {
      error (108, "ERROR: file is not open");
      return 0;
    } else if (sparse == 0) {
      error (111, "ERROR: cannot write uninitialized base_sparse");
      return 0;
    } else if (Z.sparse == 0) {
      error (114, "ERROR: pattern base_sparse is uninitialized");
      return 0;
    }
    return CHOLMOD(write_sparse) (f, sparse, Z.sparse, comments.c_str(), common->common);
  }
  //End

  //CHOLMOD Cholesky functions
  //spsolve is in base_factor.cpp
  //End
  
  //CHOLMOD Core functions
  int base_sparse::reallocate (size_t nznew) {
    return CHOLMOD(reallocate_sparse) (nznew, sparse, common->common);
  }

  long base_sparse::nnz () const {
    if (sparse == 0) {
      error (132, "ERROR: base_sparse is uninitialized");
      return 0;
    }
    return CHOLMOD(nnz) (sparse, common->common);
  }

  void base_sparse::transpose (const base_sparse & A, int values) {
    if (A.sparse == 0) {
      error (140, "ERROR: cannot transpose uninitialized base_sparse");
      return;
    }
    if (sparse != 0)
      CHOLMOD(free_sparse) (&sparse, common->common);
    common = A.common;
    sparse = CHOLMOD(transpose) (A.sparse, values, common->common);
  }

  int base_sparse::sort () {
    if (sparse == 0) {
      error (151, "ERROR: cannot sort uninitialized base_sparse");
      return 0;
    }
    return CHOLMOD(sort) (sparse, common->common);
  }

  void base_sparse::band (const base_sparse & A, long k1, long k2, int mode) {
    if (A.sparse == 0) {
      error (159, "ERROR: cannot extract band from uninitialized base_sparse");
      return;
    }
    if (sparse != 0)
      CHOLMOD(free_sparse) (&sparse, common->common);
    common = A.common;
    sparse = CHOLMOD(band) (A.sparse, k1, k2, mode, common->common);
  }

  int base_sparse::band_inplace (long k1, long k2, int mode) {
    if (sparse == 0) {
      error (170, "ERROR: cannot extract band from uninitialized base_sparse");
      return 0;
    }
    return CHOLMOD(band_inplace) (k1, k2, mode, sparse, common->common);
  }

  void base_sparse::aat (const base_sparse & A, bmInt *fset, size_t size, int mode) {
    if (A.sparse == 0) {
      error (178, "ERROR: cannot calculate AAT of uninitialized base_sparse");
      return;
    }
    cholmod_sparse *tempSparse = sparse;
    common = A.common;
    sparse = CHOLMOD(aat) (A.sparse, fset, size, mode, common->common);
    if (tempSparse != 0)
      CHOLMOD(free_sparse) (&tempSparse, common->common);
  }

  void base_sparse::copy (const base_sparse & A, int stype, int mode) {
    if (A.sparse == 0) {
      error (189, "ERROR: cannot copy uninitialized base_sparse");
      return;
    }
    if (sparse != 0)
      CHOLMOD(free_sparse) (&sparse, common->common);
    common = A.common;
    sparse = CHOLMOD(copy) (A.sparse, stype, mode, common->common);
  }
 
  void base_sparse::add (const base_sparse & A, const base_sparse & B, double alpha, double beta, int values, int sorted) {
    if (A.common != B.common) {
      error (200, "ERROR: base_sparses have different base_common objects");
      return;
    } else if (A.sparse == 0) {
      error (203, "ERROR: left side operand is uninitialized");
      return;
    } else if (B.sparse == 0) {
      error (206, "ERROR: right side operand is uninitialized");
      return;
    }
    if (sparse != 0)
      CHOLMOD(free_sparse) (&sparse, common->common);
    common = A.common;
    double Alpha[2] = {alpha, 0};
    double Beta[2] = {beta, 0};
    sparse = CHOLMOD(add) (A.sparse, B.sparse, Alpha, Beta, values, sorted, common->common);
  }

  void base_sparse::add (const base_sparse & A, const base_sparse & B, double alpha[2], double beta[2], int values, int sorted) {
    if (A.common != B.common) {
      error (200, "ERROR: base_sparses have different base_common objects");
      return;
    } else if (A.sparse == 0) {
      error (203, "ERROR: left side operand is uninitialized");
      return;
    } else if (B.sparse == 0) {
      error (206, "ERROR: right side operand is uninitialized");
      return;
    }
    if (sparse != 0)
      CHOLMOD(free_sparse) (&sparse, common->common);
    common = A.common;
    sparse = CHOLMOD(add) (A.sparse, B.sparse, alpha, beta, values, sorted, common->common);
  }

  int base_sparse::xtype_change (int to_xtype) {
    return CHOLMOD(sparse_xtype) (to_xtype, sparse, common->common);
  }

  //dense_to_sparse is in base_dense.cpp
  //factor_to_sparse is in base_factor.cpp
  //triplet_to_sparse is in base_triplet.cpp
  //End
  
  //CHOLMOD MatrixOps functions
  int base_sparse::drop (double tol) {
    if (sparse == 0) {
      error (225, "ERROR: cannot drop from uninitialized base_sparse");
      return 0;
    }
    return CHOLMOD(drop) (tol, sparse, common->common);
  }

  double base_sparse::norm (int n) const {
    if (sparse == 0) {
      error (233, "ERROR, cannot calculate norm of uninitialized base_sparse");
      return 0;
    }
    return CHOLMOD(norm_sparse) (sparse, n, common->common);
  }

  void base_sparse::horzcat (const base_sparse & A, const base_sparse & B, int values) {
    if (A.common != B.common) {
      error (241, "ERROR: base_sparses have different base_common objects");
      return;
    } else if (A.sparse == 0) {
      error (244, "ERROR: first base_sparse is uninitialized");
      return;
    } else if (B.sparse == 0) {
      error (247, "ERROR: second base_sparse is uninitialized");
    }
    if (sparse != 0)
      CHOLMOD(free_sparse) (&sparse, common->common);
    common = A.common;
    sparse = CHOLMOD(horzcat) (A.sparse, B.sparse, values, common->common);
  }

  void base_sparse::ssmult (const base_sparse & A, const base_sparse & B, int stype, int values, int sorted) {
    if (A.common != B.common) {
      error (257, "ERROR: base_sparses have different base_common objects");
      return;
    } else if (A.sparse == 0) {
      error (260, "ERROR: left side operand is uninitialized");
      return;
    } else if (B.sparse == 0) {
      error (263, "ERROR: right side operand is uninitialized");
      return;
    }
    if (sparse != 0)
      CHOLMOD(free_sparse) (&sparse, common->common);
    common = A.common;
    sparse = CHOLMOD(ssmult) (A.sparse, B.sparse, stype, values, sorted, common->common);
  }

  void base_sparse::submatrix (const base_sparse & A, bmInt *rset, long rsize, bmInt *cset, long csize, int values, int sorted) {
    if (A.sparse == 0) {
      error (274, "ERROR: cannot extract submatrix from uninitialized base_sparse");
      return;
    }
    if (sparse != 0)
      CHOLMOD(free_sparse) (&sparse, common->common);
    common = A.common;
    sparse = CHOLMOD(submatrix) (A.sparse, rset, rsize, cset, csize, values, sorted, common->common);
  }

  void base_sparse::vertcat (const base_sparse & A, const base_sparse & B, int values) {
    if (A.common != B.common) {
      error (285, "ERROR: base_sparses have different base_common objects");
      return;
    } else if (A.sparse == 0) {
      error (288, "ERROR: first base_sparse is uninitialized");
      return;
    } else if (B.sparse == 0) {
      error (291, "ERROR: second base_sparse is uninitialized");
      return;
    }
    if (sparse != 0)
      CHOLMOD(free_sparse) (&sparse, common->common);
    common = A.common;
    sparse = CHOLMOD(vertcat) (A.sparse, B.sparse, values, common->common);
  }

  int base_sparse::symmetry (int option, bmInt *xmatched, bmInt *pmatched, bmInt *nzoffdiag,
      bmInt *nzdiag) {
    if (sparse == 0) {
      error (303, "ERROR: cannot verify symmetry of uninitialized base_sparse");
      return 0;
    }
    return CHOLMOD(symmetry) (sparse, option, xmatched, pmatched, nzoffdiag, nzdiag,
        common->common);
  }
  //End

  //Other functions
  void base_sparse::print_more (std::ostream & outstr) const {
    if (sparse == 0) {
      outstr << "base_sparse not initialized\n";
      return;
    }
    typedef double * ptrd;
    typedef int * ptri;
    ptrd px = ptrd(sparse->x);
    ptri pp = ptri(sparse->p);
    ptri pi = ptri(sparse->i);
    size_t n = sparse->ncol;
    size_t nz = pp[n];
    outstr << "x = ";
    for (size_t i = 0; i < nz; i++)
      outstr << px[i] << " ";
    outstr << "\n";
    outstr << "i = ";
    for (size_t i = 0; i < nz; i++)
      outstr << pi[i] << " ";
    outstr << "\n";
    outstr << "p = ";
    for (size_t i = 0; i <= n; i++)
      outstr << pp[i] << " ";
    outstr << "\n";

  }

  void base_sparse::print_matlab (std::ostream & outstr) const {
    if (sparse == 0) {
      outstr << "base_sparse not initialized\n";
      return;
    }
    double * px = static_cast < double * > (sparse->x);
    bmInt * pi = static_cast < bmInt * > (sparse->i);
    bmInt * pp = static_cast < bmInt * > (sparse->p);
    bmInt nrow = sparse->nrow, ncol = sparse->ncol;

    outstr << "A = sparse(" << nrow << ',' << ncol << ");" << std::endl;
    for (bmInt j = 0; j < ncol; j++) {
      bmInt pstart = pp[j], pend = pp[j+1];
      while (pstart != pend) {
        outstr << "A(" << pi[pstart]+1 << ',' << j+1 << ") = "
               << px[pstart] << ';' << std::endl;
        pstart++;
      }
    }
  }

  void base_sparse::error (int line, std::string message) const {
    common->error (-10, "base_sparse.cpp", line, message);
  }

}
