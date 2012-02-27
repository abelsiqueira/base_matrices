

namespace base_matrices {

  //base_dense functions
  inline base_dense ones (base_common & cm, size_t nrow, size_t ncol, int xtype = CHOLMOD_REAL) {
    base_dense A(cm, "ones", nrow, ncol, xtype);
    return A;
  }

  inline base_dense zeros (base_common & cm, size_t nrow, size_t ncol, int xtype = CHOLMOD_REAL) {
    base_dense A(cm, "zeros", nrow, ncol, xtype);
    return A;
  }

  inline base_dense eye (base_common & cm, size_t nrow, size_t ncol, int xtype = CHOLMOD_REAL) {
    base_dense A(cm, "eye", nrow, ncol, xtype);
    return A;
  }

  inline base_dense read (base_common & cm, FILE *f) {
    base_dense A(cm);
    A.read(f);
    return A;
  }

  inline base_dense solve (int sys, const base_factor & L, const base_dense & B) {
    base_dense X( B.get_common () );
    X.solve (sys, L, B);
    return X;
  }

  inline base_dense sparse_to_dense (const base_sparse & A) {
    base_dense B( A.get_common () );
    B.sparse_to_dense (A);
    return B;
  }

  inline base_dense full (const base_sparse & A) {
    base_dense B( A.get_common () );
    B.sparse_to_dense (A);
    return B;
  }

  inline double norm (const base_dense & A, int n = 2) {
    return A.norm (n);
  }

  inline base_dense sdmult (const base_sparse & A, int transpose, double alpha[2], const base_dense & X) {
    base_dense Y( A.get_common() );
    double beta[2] = {0, 0};
    Y.sdmult (A, transpose, alpha, beta, X);
    return Y;
  }

  inline base_dense sdmult (const base_sparse & A, int transpose, double alpha[2], double beta[2], const base_dense & X, const base_dense & Y) {
    base_dense B(Y);
    B.sdmult (A, transpose, alpha, beta, X);
    return B;
  }

  inline base_dense sdmult (const base_sparse & A, int transpose, double alpha, const base_dense & X) {
    double Alpha[2] = {alpha, 0};
    return sdmult (A, transpose, Alpha, X);
  }
  
  inline base_dense sdmult (const base_sparse & A, int transpose, double alpha, double beta, const base_dense & X, const base_dense & Y) {
    double Alpha[2] = {alpha, 0};
    double Beta[2] = {beta, 0};
    return sdmult (A, transpose, Alpha, Beta, X, Y);
  }

  inline base_dense operator* (const base_sparse & A, const base_dense & X) {
    return sdmult(A, 0, 1, X);
  }

  inline double dot (const base_dense & v, const base_dense & w) {
    return v.dot(w);
  }
  //end of base_dense functions

  //base_factor functions
  inline base_factor analyze (const base_sparse & A) {
    base_factor L( A.get_common() );
    L.analyze (A);
    return L;
  }

  inline base_factor factorize (const base_sparse & A, const base_factor & L) {
    base_factor M(L);
    M.factorize (A);
    return M;
  }
  //end of base_factor functions

  //base_sparse functions
  inline base_sparse speye (const base_common & cm, size_t nrow, size_t ncol, int xtype = CHOLMOD_REAL) {
    base_sparse A(cm, "speye", nrow, ncol, 0, xtype);
    return A;
  }

  // A = speye (v, m) = speye (v, m, 0) => A : m x size(v,1), i.e. A  * v is well defined
  // A = speye (v, m, 1)                => A : size(v,1) x m, i.e. A' * v is well defined 
  inline base_sparse speye (const base_dense & v, size_t nrow, size_t transpose = 0) {
    size_t d1, d2;
    if (transpose == 0) {
      d1 = v.get_nrow ();
      d2 = nrow;
    } else {
      d1 = nrow;
      d2 = v.get_nrow ();
    }
    base_sparse A( v.get_common (), "speye", d1, d2, 0, CHOLMOD_REAL);
    return A;
  }

  //A = speye (v) => A : size(v,1) x size(v,1)
  inline base_sparse speye (const base_dense & v) {
    size_t d = v.get_nrow ();
    base_sparse A( v.get_common (), "speye", d, d, 0, CHOLMOD_REAL);
    return A;
  }

  inline base_sparse spzeros (const base_common & cm, size_t nrow, size_t ncol, size_t nzmax, int xtype = CHOLMOD_REAL) {
    base_sparse A(cm, "spzeros", nrow, ncol, nzmax, xtype);
    return A;
  }

  inline base_sparse read (const base_common & cm, FILE *f) {
    base_sparse A(cm);
    A.read(f);
    return A;
  }

  inline base_sparse spsolve (int sys, const base_factor & L, const base_sparse & B) {
    base_sparse X( L.get_common () );
    X.spsolve(sys, L, B);
    return X;
  }

  inline double nnz (const base_sparse & A) {
    return A.nnz();
  }

  inline base_sparse transpose (const base_sparse & A, int values = 2) {
    base_sparse B( A.get_common() );
    B.transpose (A, values);
    return B;
  }

  inline base_sparse band (const base_sparse & A, UF_long k1, UF_long k2, int mode) {
    base_sparse B( A.get_common() );
    B.band (A, k1, k2, mode);
    return B;
  }

  inline base_sparse aat (const base_sparse & A, bmInt *fset = 0, size_t size = 0, int mode = 1) {
    base_sparse B( A.get_common() );
    B.aat (A, fset, size, mode);
    return B;
  }

  inline base_sparse add (const base_sparse & A, const base_sparse & B, double alpha = 1, double beta = 1, int values = 1, int sorted = 0) {
    base_sparse C( A.get_common() );
    C.add (A, B, alpha, beta, values, sorted);
    return C;
  }

  inline base_sparse add (const base_sparse & A, const base_sparse & B, double alpha[2], double beta[2], int values = 1, int sorted = 0) {
    base_sparse C( A.get_common() );
    C.add (A, B, alpha, beta, values, sorted);
    return C;
  }

  inline base_sparse operator+ (const base_sparse & A, const base_sparse & B) {
    return add (A, B);
  }

  inline base_sparse operator- (const base_sparse & A, const base_sparse & B) {
    return add (A, B, 1, -1);
  }

  inline double norm (const base_sparse & A, int n = 0) {
    return A.norm(n);
  }

  inline base_sparse horzcat (const base_sparse & A, const base_sparse & B, int values = 1) {
    base_sparse C( A.get_common() );
    C.horzcat (A, B, values);
    return C;
  }

  inline base_sparse ssmult (const base_sparse & A, const base_sparse & B, int stype = 0, int values = 1, int sorted = 0) {
    base_sparse C( A.get_common() );
    C.ssmult (A, B, stype, values, sorted);
    return C;
  } 

  inline base_sparse operator* (const base_sparse & A, const base_sparse & X) {
    return ssmult (A, X);
  }

  inline base_sparse submatrix (const base_sparse & A, bmInt *rset, UF_long rsize, bmInt *cset, UF_long csize, int values = 1, int sorted = 0) {
    base_sparse B( A.get_common() );
    B.submatrix (A, rset, rsize, cset, csize, values, sorted);
    return B;
  }

  inline base_sparse vertcat (const base_sparse & A, const base_sparse & B, int values = 1) {
    base_sparse C( A.get_common() );
    C.vertcat (A, B, values);
    return C;
  }
  //end of base_sparse functions
}
