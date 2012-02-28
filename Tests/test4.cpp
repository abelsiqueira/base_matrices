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
/* base_matrices: test4
 * 
 * This test is made with the retangular matrices and the extra functions
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
  base_dense b(c), x(c);
  base_factor L(c);
  b.read(stdin);

  L = analyze (A);
  L = factorize (A,L);
  x = solve (CHOLMOD_A, L, b);
  cout << "|Ax - b| = " << ( A*sdmult(A, 1, 1, x) - b ).norm() << endl;
  cout << "|x|^2 - dot(x,x) = " << norm(x,2) * norm(x,2) - dot (x,x) << endl << endl;
//  full(A*speye(c, 11, 1)).print_more ();
//  add( A, speye (c, A.get_nrow(), A.get_ncol()), 1, -1, 1, 1).print_more();

  return (0) ;
}
