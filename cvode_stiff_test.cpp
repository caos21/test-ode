/*
The contents of this file is free and unencumbered software released into the
public domain. For more information, please refer to <http://unlicense.org/>
*/

#include <iostream>
#include <cmath>

#include <cvode/cvode.h>
#include <cvode/cvode_dense.h>
#include <cvode/cvode_diag.h>
#include <nvector/nvector_serial.h>
#include <sundials/sundials_types.h>

// typedef int (*CVRhsFn)(realtype t, N_Vector y, N_Vector ydot, void *user_data);


int function(realtype t, N_Vector x, N_Vector dxdt, void *user_data);

class Solver{
public:
  Solver(N_Vector xini_,
         realtype ti_,
         realtype a_=-101.,
         realtype b_=-100.,
         realtype reltol_=1.e-6,
         realtype abstol_=1.e-6) :
         xini(xini_),
         ti(ti_),
         a(a_),
         b(b_),
         reltol(reltol_),
         abstol(abstol_) {
    std::cerr << "\n Initial condition : "
              << NV_Ith_S(xini,0) << '\t' << NV_Ith_S(xini,1) << '\n'; 
    cvode_mem = CVodeCreate(CV_BDF, CV_NEWTON);    
    flag = CVodeInit(cvode_mem, function, ti, xini);
    std::cerr << "\n[ii] flag : " << flag << "\n";
    flag = CVodeSStolerances(cvode_mem, reltol, abstol);
    std::cerr << "\n[ii] flag : " << flag << "\n";
    flag = CVodeSetUserData(cvode_mem, static_cast< void* >(this));
    std::cerr << "\n[ii] flag : " << flag << "\n";
    flag = CVDense(cvode_mem, 2);
    std::cerr << "\n[ii] flag : " << flag << "\n";
    flag = CVDlsSetDenseJacFn(cvode_mem, NULL);
//     flag = CVDiag(cvode_mem);
    std::cerr << "\n[ii] flag : " << flag << "\n";
    
    // allocate vector for output
    y = NULL;
    y = N_VNew_Serial(2);
  }
  
  ~Solver() {
    N_VDestroy_Serial(y);
    CVodeFree(&cvode_mem);
  }
  
  // declaring friend function grants access to a and b from class
  friend int function(realtype t, N_Vector x, N_Vector dxdt, void *user_data);  
  
  realtype a;
  realtype b;
  realtype reltol;
  realtype abstol;
  realtype ti;
  N_Vector xini;
  N_Vector y;

  void compute(realtype tf, realtype dt) {
    realtype t, tout;
    for(iout=1, tout=ti+dt; tout <= tf; iout++, tout += dt) {
      flag = CVode(cvode_mem, tout, y, &t, CV_NORMAL);
      std::cout << tout << '\t' << NV_Ith_S(y,0) << '\t' << NV_Ith_S(y,1) << '\n'; 
      if (flag != CV_SUCCESS) {
        std::cerr << "\n[ee] Terminate. Flag : " << flag << "\n";
        std::terminate();
      }      
    }
    
  }
  
  void compute_step(realtype ts, realtype dt) {
    flag = CVode(cvode_mem, ts+dt, y, &t, CV_NORMAL);
    std::cerr << t << '\t' << NV_Ith_S(y,0) << '\t' << NV_Ith_S(y,1) << '\n'; 
    if (flag != CV_SUCCESS) {
      std::cerr << "\n[ee] Terminate. Flag : " << flag << "\n";
      std::terminate();
    }    
  }
  
private:
  void *cvode_mem = NULL;
  int flag;
  realtype t;
  int iout;
};


int main(int argc, char **argv) {

  // initial conditions, in boost vector_type x( 2 , 1.0 );
  N_Vector xini = NULL;
  xini = N_VNew_Serial(2);
  NV_Ith_S(xini, 0) = RCONST(2.0);
  NV_Ith_S(xini, 1) = RCONST(1.0);
 
  // Instantiate the solver
  Solver sol(xini, 0.0);
  
  // final time
  realtype tf = RCONST(50.0);
  // delta time
  realtype dt = RCONST(0.01);
  sol.compute(tf, dt);


  // or solve one step
//   sol.compute_step(0.0, dt);
  
  // free memory for xini
  N_VDestroy_Serial(xini);
  
  return 0;  
}

// function to integrate
int function(realtype t, N_Vector x, N_Vector dxdt, void *user_data) {
  
  // static_cast user_data to class Solver 
  Solver* s = static_cast< Solver* >(user_data); 
  
  realtype x0, x1;
  
  x0 = NV_Ith_S(x, 0);
  x1 = NV_Ith_S(x, 1);
  
  NV_Ith_S(dxdt, 0) = s->a*x0 + s->b*x1; // in boost dxdt[ 0 ] = -101.0 * x[ 0 ] - 100.0 * x[ 1 ];
  NV_Ith_S(dxdt, 1) = x0; // in boost dxdt[ 1 ] = x[ 0 ];

  return(0);
} 
