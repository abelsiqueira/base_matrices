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
#include <fstream>
#include <iostream>
#include <cstdlib>
#include <cstdlib>
#include <cstdlib>
#include <string>

using namespace std;

int main (int argc, char **argv) {

  if (argc < 2)
    return 1;

  string filename("ret_mm" + string(argv[1]) + ".mtx");
  ofstream file(filename.c_str());
  int m = atoi(argv[1]);
  int n = m+1;
  
  file << m << " " << n << " " << 2*m << endl;

  for (int i = 0; i < m; i++) {
    file << i << " " << i << " 1\n";
    file << i << " " << i+1 << " -1\n";
  }

  return 0;

}
