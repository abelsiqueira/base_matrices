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
/* base_matrices: test3
 * 
 * This test is made with the retangular matrices.
 * 
 */

#include "base_matrices.h"
#include <iostream>
#include <cstdio>

using namespace base_matrices;
using namespace std;

int main () {

  base_common c;
  base_sparse A(c);
  A.read(stdin); //A is m by n. We intend to solve AA'x = b
  int m = A.get_nrow();
  double alpha[2] = {1,0}, beta[2] = {0,0};
  base_dense e(c, "ones", m, 1);
  base_dense b(c), x(c), dif(c);
  base_factor L(c);

  L.analyze (A);
  L.factorize (A);
  b.sdmult (A, 1, alpha, beta, e); //temp <- A'*e
  b.sdmult (A, 0, alpha, beta, b); //b <- A*temp = A*A'*e
  x.solve (CHOLMOD_A, L, b);
  dif = x - e;
  cout << "|x - e| = " << dif.norm() << endl << endl;

  return (0) ;
}
