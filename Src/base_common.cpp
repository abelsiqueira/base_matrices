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
#include "base_common.h"

namespace base_matrices {

  base_common::base_common () {
    common = new cholmod_common;
    CHOLMOD(start)(common);
  }

  base_common::~base_common () {
    CHOLMOD(finish)(common);
    delete common;
  }

  //CHOLMOD Check functions
  int base_common::check () const {
    return CHOLMOD(check_common) (common);
  }

  int base_common::print (std::string name) const {
    return CHOLMOD(print_common) (name.c_str(), common);
  }
  //End

  //CHOLMOD Core functions
  int base_common::defaults () {
    return CHOLMOD(defaults) (common);
  }

  size_t base_common::maxrank (size_t n) {
    return CHOLMOD(maxrank) (n, common);
  }

  int base_common::allocate_work (size_t nrow, size_t iworksize, size_t xworksize) {
    return CHOLMOD(allocate_work) (nrow, iworksize, xworksize, common);
  }

  int base_common::free_work () {
    return CHOLMOD(free_work) (common);
  }

  UF_long base_common::clear_flag () {
    return CHOLMOD(clear_flag) (common);
  }

  int base_common::error (int status, std::string file, int line, std::string message) const {
    return CHOLMOD(error) (status, file.c_str(), line, message.c_str(), common);
  }
  //End

  //Other functions
  void base_common::set_error_handler (void (*error_handler) (int status, const char *file,
      int line, const char *message) ) {
    common->error_handler = error_handler;
  }
  
}
