cmake_minimum_required (VERSION 2.6)

# sundials test
set(CVODE_LIB sundials_cvode)
set(NVECS_LIB sundials_nvecserial)
set(EXTRA_LINK_LIBS m)

set(SUNDIALS_LIBS ${CVODE_LIB} ${NVECS_LIB} ${EXTRA_LINK_LIBS})

add_executable(cvode_stiff_test cvode_stiff_test.cpp)
target_link_libraries(cvode_stiff_test ${SUNDIALS_LIBS})


# boost example
add_executable(stiff_system stiff_system.cpp)

