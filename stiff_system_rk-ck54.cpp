/*
Copyright <2017> <Benjamin Santos>

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

#include <iostream>
#include <fstream>
#include <utility>
#include <vector>

#include <boost/numeric/odeint.hpp>

#include <boost/phoenix/core.hpp>

#include <boost/phoenix/core.hpp>
#include <boost/phoenix/operator.hpp>

using namespace std;
using namespace boost::numeric::odeint;
namespace phoenix = boost::phoenix;

// vector type
typedef double value_type;
typedef vector< value_type > state_type;


// system definition
struct stiff_system{
  
  // constructor, sets a and b
  stiff_system(double a_=-101.0, double b_=-100.0): a(a_), b(b_){  }
  
  // rhs of ode system
  inline
  void operator()( const state_type &x, state_type &dxdt, const value_type t) const {
    dxdt[0] = a*x[0] + b*x[1];
    dxdt[1] = x[0];
  }
  double a;
  double b;
};


int main(int argc, char **argv) {

  state_type xini( 2.0 , 1.0 );

  // system
  //   stiff_system ssys();
  stiff_system* ssys = new stiff_system();

  auto stepper = make_controlled(1.0e-6, 1.0e-6,
                                runge_kutta_cash_karp54< state_type >());


  size_t num_of_steps = integrate_adaptive(stepper, *ssys, xini,
                                           0.01, 50.0, 0.01,
                        cout << phoenix::arg_names::arg2 << '\t'
                             << phoenix::arg_names::arg1[0] << '\t'
                             << phoenix::arg_names::arg1[1] << '\t'
                             << phoenix::arg_names::arg1[1]
                            *phoenix::arg_names::arg1[1] << "\n" );

  cerr << "\n[ii] Number of steps: " <<  num_of_steps << endl;

  // deallocate memory for ssys
  delete(ssys);
  
  return 0;
}
