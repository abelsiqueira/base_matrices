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
/* base_matrices: test1
 * 
 * This test reads a matrix from stdin and solve the system Ax = b, with
 * b = [1 1 ... 1]'.
 *
 */

#include "base_matrices.h"

using namespace base_matrices;
using namespace std;

int main (void) {

  double one [2] = {1,0}, m1 [2] = {-1,0} ;
  base_common c;
  base_sparse A(c);
  base_dense x(c);
  base_factor L(c);
  base_dense r(c);

  A.read(stdin); //reads from stdin 
  A.print("A");

  //creates b = [1 1 ... 1]' compatible with A
  base_dense b(c, "ones", A.get_nrow(), 1); 

  L.analyze(A); //Analyze A to make factorization
  L.factorize(A); //Factorize A into L
  x.solve(CHOLMOD_A,L,b); //Solve the system Ax = b, using L
  r = b; 
  r.sdmult (A, 0, one, m1, x); //r <- one * A*x + m1 *r = Ax - b
  cout << "|Ax - b| = " << r.norm() << endl << endl;

  return 0 ;
}
