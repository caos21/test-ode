# test-ode
Sundials CVODE and BOOST test

## Requirements

### Ubuntu

- Install cmake, boost and sundials
 ```
 sudo apt install cmake g++ libboost-all-dev libsundials-serial-dev
 ```
 
 - Optional, xmgrace for plotting
 ```
 sudo apt install xmgrace
 ```
 
 ## Compile
 
 - Create build
 ```
 cmake -Bbuild -H.
 ```
 
 - Change to directory build
 ```
 cd build
 ```
 
 - Compile and link
 ```
 make
 ```
 
 - Two executables should be created `stiff_system` from boost documentation and `cvode_stiff_test`.
 
 ## Run
 
 - Both programs produce output to `stdout`.
 
 - The following reproduce this graph,
 ![alt text]( https://github.com/caos21/test-ode/blob/master/comp.png  "Results")
 
 
 1. Run cvode test and redirect output to file, e.g. `cvode.dat` 
 ```
 ./cvode_stiff_test > cvode.dat
 ```
 
 2. Run boost test and redirect output to file, e.g. `boost.dat` 
 ```
 ./stiff_system > boost.dat 
 ```
 
 3. Plot results
 ```
 xmgrace -par comp.par -block cvode.dat -bxy 1:4  -block boost.dat -bxy 1:3
 ```
 
