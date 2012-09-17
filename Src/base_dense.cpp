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
#include "base_dense.h"
#include <cassert>
#include <cmath>

namespace base_matrices {

  base_dense::base_dense (const base_common & cm, size_t nrow, size_t ncol, size_t d,
      int xtype) {
    common = &cm;
    if ( (nrow == 0) || (ncol == 0) ) {
      dense = 0;
    } else if (d < nrow) {
      dense = CHOLMOD(allocate_dense) (nrow, ncol, nrow, xtype, common->common);
    } else {
      dense = CHOLMOD(allocate_dense) (nrow, ncol, d, xtype, common->common);
    }
  }

  base_dense::base_dense (const base_common & cm, std::string type, size_t nrow, size_t ncol, int xtype) {
    common = &cm;
    if ( (nrow == 0) || (ncol == 0) )
      dense = 0;
    else if (type == "zeros")
      dense = CHOLMOD(zeros) (nrow, ncol, xtype, common->common);
    else if (type == "ones")
      dense = CHOLMOD(ones) (nrow, ncol, xtype, common->common);
    else if (type == "eye")
      dense = CHOLMOD(eye) (nrow, ncol, xtype, common->common);
    else
      dense = 0;
  }

  base_dense::base_dense (const base_dense & b) {
    common = b.common;
    if (b.dense == 0) {
      error (34, "ERROR: copying uninitialized base_dense");
      dense = 0;
    } else {
      dense = CHOLMOD(copy_dense) (b.dense, common->common);
    }
  }

  void base_dense::operator= (const base_dense & b) {
    if (b.dense == 0) {
      error (43, "ERROR: cannot assign to uninitialized base_dense");
      return;
    }
    if (b.dense == dense)
      return;
    if (dense == 0) {
      common = b.common;
      dense = CHOLMOD(copy_dense) (b.dense, common->common);
    } else {
      size_t nrow = dense->nrow, ncol = dense->ncol, d = dense->d;
      int xtype = dense->xtype;
      if ( (nrow != b.dense->nrow) || (ncol != b.dense->ncol) || 
          (d != b.dense->d) || (xtype != b.dense->xtype) ) {
        CHOLMOD(free_dense) (&dense, common->common);
        common = b.common;
        dense = CHOLMOD(copy_dense) (b.dense, common->common);
      } else {
        double * pthisx = static_cast < double * > (dense->x);
        double * pbx = static_cast < double * > (b.dense->x);

        for (size_t i = 0; i < nrow; i++)
          for (size_t j = 0; j < ncol; j++)
            pthisx[i + j*d] = pbx[i + j*d];
      }
    }
  }

  base_dense::~base_dense () {
    CHOLMOD(free_dense) (&dense, common->common);
  }

  //CHOLMOD Check functions
  int base_dense::check () const {
    return CHOLMOD(check_dense) (dense, common->common);
  }

  int base_dense::print (std::string name) const {
    return CHOLMOD(print_dense) (dense, name.c_str(), common->common);
  }

  void base_dense::read (FILE *f) {
    if (f == 0) {
      error (63, "ERROR: file not open");
      return;
    }
    if (dense != 0)
      CHOLMOD(free_dense) (&dense, common->common);
    dense = CHOLMOD(read_dense) (f, common->common);
  }

  int base_dense::write (FILE *f) const {
    if (f == 0) {
      error (73, "ERROR: file not open");
      return 0;
    } else if (dense == 0) {
      error (78, "ERROR: base_dense uninitialized");
      return 0;
    }
    return CHOLMOD(write_dense) (f, dense, 0, common->common);
  }

  int base_dense::write (FILE *f, std::string comments) const {
    if (f == 0) {
      error (81, "ERROR: file not open");
      return 0;
    } else if (dense == 0) {
      error (89, "ERROR: base_dense uninitialized");
      return 0;
    }
    return CHOLMOD(write_dense) (f, dense, comments.c_str(), common->common);
  }
  //End

  //CHOLMOD Cholkesy functions
  void base_dense::solve (int sys, const base_factor & L, const base_dense & B) {
    if (L.common != B.common) {
      error (91, "ERROR: base_factor and base_dense have different base_common objects");
      return;
    } else if (L.factor == 0) {
      error (94, "ERROR: base_factor uninitialized");
      return;
    } else if (B.dense == 0) {
      error (97, "ERROR: base_dense uninitialized");
      return;
    }
    if (dense != 0)
      CHOLMOD(free_dense) (&dense, common->common);
    common = L.common;
    dense = CHOLMOD(solve) (sys, L.factor, B.dense, common->common);
  }
  //End

  //CHOLMOD Core functions
  void base_dense::sparse_to_dense (const base_sparse & A) {
    if (A.sparse == 0) {
      error (110, "ERROR: base_sparse uninitialized");
      return;
    }
    if (dense != 0)
      CHOLMOD(free_dense) (&dense, common->common);
    common = A.common;
    dense = CHOLMOD(sparse_to_dense) (A.sparse, common->common);
  }

  int base_dense::xtype_change (int to_xtype) {
    return CHOLMOD(dense_xtype) (to_xtype, dense, common->common);
  }

  //From Sparse
  void base_sparse::dense_to_sparse (const base_dense & A, int values) {
    if (A.dense == 0) {
      A.error (126, "ERROR: base_dense uninitialized");
      return;
    }
    if (sparse != 0)
      CHOLMOD(free_sparse) (&sparse, common->common);
    common = A.common;
    sparse = CHOLMOD(dense_to_sparse) (A.dense, values, common->common);
  }
  //End

  //CHOLMOD MatrixOps functions
  double base_dense::norm (int n) const {
    if (dense == 0) {
      error (139, "ERROR: base_dense uninitialized");
      return 0;
    }

    return CHOLMOD(norm_dense) (dense, n, common->common);
  }

  int base_dense::sdmult (const base_sparse & A, int transpose, double alpha[2], double beta[2], const base_dense & X) {
    //If Y (this) is emtpy we create it.
    //If Y does not match the resulting multiplication, we erase it and create it.
    //If Y match, then we use it.
    if (A.sparse == 0) {
      error (150, "ERROR: base_sparse uninitialized");
      return 0;
    } else if (X.dense == 0) {
      error (153, "ERROR: base_dense uninitialized");
      return 0;
    }
    uint nx, ny;
    if (transpose) {
      nx = A.get_ncol();
      ny = A.get_nrow();
    } else {
      nx = A.get_nrow();
      ny = A.get_ncol();
    }

    if ( (ny != X.dense->nrow) || (A.common != X.common) )
      return 0;

    cholmod_dense *Xdense = X.dense;
    if (dense != 0) {
      if (this == &X) {
        dense = 0;
      } else if ( (dense->nrow != nx) || (dense->ncol != Xdense->ncol) || (dense->xtype != Xdense->xtype) ) {
        CHOLMOD(free_dense) (&dense, common->common);
        dense = 0;
      }
    }

    common = X.common;
    if (dense == 0) {
      dense = CHOLMOD(allocate_dense) (nx, Xdense->ncol, nx, Xdense->xtype, common->common);
      int vecsize = nx*(Xdense->ncol);
      typedef double * ptrd;
      ptrd px = ptrd(dense->x);
      for (int i = 0; i < vecsize; i++)
        px[i] = 0;
    }

    int exitflag = CHOLMOD(sdmult) (A.sparse, transpose, alpha, beta, Xdense, dense, common->common);
    if (this == &X)
      CHOLMOD(free_dense) (&Xdense, common->common);
    return exitflag;
  }
  //End
  
  //Other functions
  double base_dense::dot (const base_dense & v) const {
    //This must only work if the number of columns is 1, and
    //if the matrices are real
    if (v.dense == 0) {
      error (212, "ERROR: base_dense uninitialized");
      return 0;
    } else if ( (dense->ncol != 1) || (v.dense->ncol != 1) ) {
      error (215, "ERROR: dot only works for base_dense with 1 column");
      return 0;
    } else if ( (dense->xtype != CHOLMOD_REAL) || (v.dense->xtype != CHOLMOD_REAL) ) {
      error (218, "ERROR: dot only works with real matrices");
      return 0;
    } else if (common != v.common) {
      error (221, "ERROR: base_dense have different base_common objects");
      return 0;
    } else if (dense->nrow != v.dense->nrow) {
      error (224, "ERROR: number of rows must match");
      return 0;
    }
    double val = 0.0;
    size_t n = dense->nrow;
    typedef double * ptrd;
    ptrd pthis = ptrd(dense->x), pv = ptrd(v.dense->x);
    for (size_t i = 0; i < n; i++)
      val += pthis[i]*pv[i];
    return val;

  }
  void base_dense::print_more (std::ostream & outstr) const {
    if (dense == 0) {
      outstr << "base_dense not initialized\n";
      return;
    }
    typedef double * ptrd;
    ptrd p = ptrd(dense->x);
    int m = dense->nrow, n = dense->ncol;
    int d = dense->d;
    for (int i = 0; i < m; i++) {
      for (int j = 0; j < n; j++) {
        outstr.width(3);
        outstr << p[i + j*d] << " ";
      }
      outstr << "\n";
    }
  }

  void base_dense::saxpy (const base_dense & B, double alpha) {
    if (alpha == 0)
      return;

    if (dense == 0) {
      error (-1, "ERROR: -1");
      return;
    } else if (B.dense == 0) {
      error (-1, "ERROR: -1");
      return;
    }

    size_t row = B.dense->nrow, col = B.dense->ncol;
    size_t dB = B.dense->d, dthis = dense->d;
    int xt = B.dense->xtype;

    if (common != B.common) {
      error (-1, "ERROR, -1");
      return;
    } else if ( (dense->nrow != row) || (dense->ncol != col) ) {
      error (-1, "ERROR, -1");
      return;
    } else if ( (xt != CHOLMOD_REAL) || (dense->xtype != CHOLMOD_REAL) ) {
      error (-1, "ERROR, -1");
      return;
    }

    double * pthisx = static_cast < double * > (this->dense->x);
    double * pBx = static_cast < double * > (B.dense->x);

    for (size_t i = 0; i < row; i++)
      for (size_t j = 0; j < col; j++)
        pthisx[i + j*dthis] += alpha * pBx[i + j*dB];

  }

  void base_dense::scale (double alpha) {
    if (dense == 0) {
      error (__LINE__, "ERROR: uninitialized base_dense");
      return;
    } else if (alpha == 1.0)
      return;

    double * px = static_cast < double * > (dense->x);
    size_t nrow = dense->nrow, ncol = dense->ncol, d = dense->d;

    for (size_t i = 0; i < nrow; i++)
      for (size_t j = 0; j < ncol; j++)
        px[i + j*d] *= alpha;
  }

  void base_dense::scale (const base_dense & A, double alpha) {
    if (A.dense == 0) {
      error (__LINE__, "ERROR: uninitialized base_dense");
      return;
    }

    double * pAx = static_cast < double * > (A.dense->x);
    size_t nrow = A.dense->nrow, ncol = A.dense->ncol, d = A.dense->d;
    int xtype = A.dense->xtype;

    if (this->dense != 0) {
      if ( (nrow != dense->nrow) || (ncol != dense->ncol) || (d < dense->d) || (xtype != dense->xtype) || (common != A.common) ) {
        CHOLMOD(free_dense) (&dense, common->common);
        common = A.common;
      }
    }

    if (this->dense == 0)
      this->dense = CHOLMOD(allocate_dense) (nrow, ncol, d, xtype, common->common);
    double * pthisx = static_cast < double * > (this->dense->x);

    for (size_t i = 0; i < nrow; i++)
      for (size_t j = 0; j < ncol; j++)
        pthisx[i + j*d] = alpha * pAx[i + j*d];
  }

  double base_dense::get (size_t i, size_t j) const {
    if ( dense == 0 ) {
      error (__LINE__, "ERROR: Uninitialized base_dense");
      return 0;
    }
    if ( dense->xtype != CHOLMOD_REAL ) {
      error (__LINE__, "ERROR: not double");
      return 0;
    }

    size_t nrow = dense->nrow, ncol = dense->ncol, d = dense->d;
    double * pthisx = static_cast < double * > (this->dense->x);

    if ( (i == 0) || (i > nrow) || (j == 0) || (j > ncol) ) {
      error (__LINE__, "ERROR: Out of Range");
      return 0;
    }
    i--; j--;

    return pthisx[i + j*d];
  }

  void base_dense::set (size_t i, size_t j, double x) {
    if ( dense == 0 ) {
      error (__LINE__, "ERROR: Uninitialized base_dense");
      return;
    }
    if ( dense->xtype != CHOLMOD_REAL ) {
      error (__LINE__, "ERROR: not double");
      return;
    }

    size_t nrow = dense->nrow, ncol = dense->ncol, d = dense->d;
    double * pthisx = static_cast < double * > (this->dense->x);

    if ( (i == 0) || (i > nrow) || (j == 0) || (j > ncol) ) {
      error (__LINE__, "ERROR: Out of Range");
      return;
    }
    i--; j--;

    pthisx[i + j*d] = x;

    return;
  }

  base_dense base_dense::operator+ (const base_dense & A) const {
    if (dense == 0) {
      error (216, "ERROR: left side operand is uninitialized");
      return *this;
    } else if (A.dense == 0) {
      error (219, "ERROR: right side operand is uninitialized");
      return *this;
    }

    size_t row = A.dense->nrow, col = A.dense->ncol;
    size_t dA = A.dense->d, dthis = dense->d;
    int xt = A.dense->xtype;

    if (common != A.common) {
      error (225, "ERROR: operands have different base_common objects");
      return *this;
    } else if ( (dense->nrow != row) || (dense->ncol != col) ) {
      error (228, "ERROR: operands have different dimensions");
      return *this;
    } else if ( (xt != CHOLMOD_REAL) || (dense->xtype != CHOLMOD_REAL) ) {
      error (271, "ERROR: operator+ defined for CHOLMOD_REAL only");
      return *this;
    }

    size_t dR = max (dA, dthis);

    base_dense R(*common, row, col, dR, xt);
    double * pthisx = static_cast < double * > (dense->x);
    double * pAx = static_cast < double * > (A.dense->x);
    double * pRx = static_cast < double * > (R.dense->x);

    for (size_t i = 0; i < row; i++)
      for (size_t j = 0; j < col; j++)
        pRx[i + j*dR] = pthisx[i + j*dthis] + pAx[i + j*dA];

    return R;
  }

  base_dense base_dense::operator- (const base_dense & A) const {
    if (dense == 0) {
      error (244, "ERROR: left side operand is uninitialized");
      return *this;
    } else if (A.dense == 0) {
      error (247, "ERROR: right side operand is uninitialized");
      return *this;
    }

    size_t row = A.dense->nrow, col = A.dense->ncol;
    size_t dA = A.dense->d, dthis = dense->d;
    int xt = A.dense->xtype;

    if (common != A.common) {
      error (253, "ERROR: operands have different base_common objects");
      return *this;
    } else if ( (dense->nrow != row) || (dense->ncol != col) ) {
      error (256, "ERROR: operands have different dimensions");
      return *this;
    } else if ( (xt != CHOLMOD_REAL) || (A.dense->xtype != CHOLMOD_REAL) ) {
      error (312, "ERROR: operator- defined for CHOLMOD_REAL only");
      return *this;
    }
    
    size_t dR = max (dA, dthis);

    base_dense R(*common, row, col, dR, xt);
    double * pthisx = static_cast < double * > (dense->x);
    double * pAx = static_cast < double * > (A.dense->x);
    double * pRx = static_cast < double * > (R.dense->x);

    for (size_t i = 0; i < row; i++)
      for (size_t j = 0; j < col; j++)
        pRx[i + j*dR] = pthisx[i + j*dthis] - pAx[i + j*dA];

    return R;
  }

  base_dense base_dense::operator* (double alpha) const {
    if (dense == 0) {
      error (332, "ERROR: base_dense is uninitialized");
      return *this;
    }

    base_dense R(*this);
    size_t row = R.dense->nrow, col = R.dense->ncol;
    size_t d = R.dense->d;
    double *px = static_cast < double * > (R.dense->x);

    for (size_t i = 0; i < row; i++)
      for (size_t j = 0; j < col; j++)
        px[i + j*d] = alpha * px[i + j*d];

    return R;
  }

  base_dense operator* (double alpha, const base_dense & A) {
    if (A.dense == 0) {
      A.error (332, "ERROR: base_dense is uninitialized");
      return A;
    }

    base_dense R(A);
    size_t row = R.dense->nrow, col = R.dense->ncol;
    size_t d = R.dense->d;
    double *px = static_cast < double * > (R.dense->x);

    for (size_t i = 0; i < row; i++)
      for (size_t j = 0; j < col; j++)
        px[i + j*d] = alpha * px[i + j*d];

    return R;
  }

  void base_dense::error (int line, std::string message) const {
    common->error (-10, "base_dense.cpp", line, message);
  }


}
