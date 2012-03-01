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
/* base_matrices: test_dense
 * 
 * This tests base_dense's routines
 * b = [1 1 ... 1]'.
 *
 */

#include "base_matrices.h"
#include <cassert>

using namespace base_matrices;
using namespace std;

int main (void) {

//  double one [2] = {1,0}, m1 [2] = {-1,0} ;
  size_t nrow = 10, ncol = 5;
  base_common c;
  base_dense x(c, nrow, ncol);
  base_dense e(c, "ones", ncol, 1);
  base_dense y(e), z(c);

  for (size_t i = 0; i < ncol; i++)
    y.set (i + 1, 1, static_cast <double> (i + 1));

  int I = ncol * (ncol + 1) / 2;
  int S = I * (2 * ncol + 1) / 3;
  z = e;
  assert (y.dot(z) == I);
  z = y;
  assert (y.dot(z) == S);
  z = -1 * z;
  assert (y.dot(z) == -S);
  z = e + y;
  assert (y.dot(z) == S + I);
  z = y - e;
  assert (y.dot(z) == S - I);

  return (0) ;

}
