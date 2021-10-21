/*
 *    This file is part of CasADi.
 *
 *    CasADi -- A symbolic framework for dynamic optimization.
 *    Copyright (C) 2010-2014 Joel Andersson, Joris Gillis, Moritz Diehl,
 *                            K.U. Leuven. All rights reserved.
 *    Copyright (C) 2011-2014 Greg Horn
 *
 *    CasADi is free software; you can redistribute it and/or
 *    modify it under the terms of the GNU Lesser General Public
 *    License as published by the Free Software Foundation; either
 *    version 3 of the License, or (at your option) any later version.
 *
 *    CasADi is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *    Lesser General Public License for more details.
 *
 *    You should have received a copy of the GNU Lesser General Public
 *    License along with CasADi; if not, write to the Free Software
 *    Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 *
 */


/**
 *  Example program demonstrating the usage of C-code generated from CasADi
 *  Generated code is encapsulated in self-contained code with the entry points
 *  defined below.
 *  Note that other usage, e.g. accessing the internal data structures in the
 *  generated files is not recommended and subject to change.
 *
 *  We show how the generated code can be used from C (or C++), without requiring
 *  any CasADi classes as well as how to use it from C++ using CasADi's external.
 *
 *  Joel Andersson, 2013-2015
 */

#include <stdio.h>

// C++ (and CasADi) from here on
#include <casadi/casadi.hpp>
using namespace casadi;
using namespace std;

// void usage_cplusplus(){
//   cout << "---" << endl;
//   cout << "Usage from CasADi C++:" << endl;
//   cout << endl;

//   // Use CasADi's "external" to load the compiled function
//   Function f = external("f");

//   // Use like any other CasADi function
//   vector<double> x = {1, 2, 3, 4};
//   vector<DM> arg = {reshape(DM(x), 2, 2), 5};
//   vector<DM> res = f(arg);

//   cout << "result (0): " << res.at(0) << endl;
//   cout << "result (1): " << res.at(1) << endl;
// }


int main(){

  // Variables
  SX x = SX::sym("x", 2, 2);
  SX y = SX::sym("y");

  // Simple function
  Function f("f", {x, y}, {sqrt(y)-1, sin(x)-y});
  Function g("g", {x, y}, {sqrt(y)-1, cos(x)-y});

  auto C = CodeGenerator("gen.c");
  // Generate C-code
  C.add(f);
  C.add(g);

  C.generate();

  // Compile the C-code to a shared library
  string compile_command = "gcc -fPIC -shared -O3 gen.c -o gen.so";
  int flag = system(compile_command.c_str());
  casadi_assert(flag==0, "Compilation failed");

  auto f_loaded = external("f","gen.so");
  auto g_loaded = external("g","gen.so");

  vector<double> x_arg = {1, 2, 3, 4};
  vector<DM> arg = {reshape(DM(x_arg), 2, 2), 5};
  vector<DM> f_res = f_loaded(arg);
  vector<DM> g_res = g_loaded(arg);

  cout << "f_result (0): " << f_res.at(0) << endl;
  cout << "f_result (1): " << f_res.at(1) << endl;

  cout << "g_result (0): " << g_res.at(0) << endl;
  cout << "g_result (1): " << g_res.at(1) << endl;
  // Usage from C++
//   usage_cplusplus();

  return 0;
}
