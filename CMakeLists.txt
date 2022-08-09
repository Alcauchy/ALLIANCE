cmake_minimum_required(VERSION 3.20)
project(Alliance_alpha_1_0 C)

set(CMAKE_C_STANDARD 11)
find_package(MPI REQUIRED)
include_directories(${MPI_INCLUDE_PATH})
#find_package(HDF5 REQUIRED)
#set(HDF5_PREFER_PARALLEL TRUE)
include_directories("/home/alcauchy/hdf5/hdf5/include")
set(CMAKE_C_FLAGS_RELEASE "-O3")
set(CMAKE_C_FLAGS_DEBUG "-O3")
add_executable(Alliance_alpha_1_0
        src/main.c
        src/utils_mpi.h src/utils_mpi.c
        src/utils_fftw.c src/utils_fftw.h
        src/utils_hdf.h src/utils_hdf.c
        src/parameters_io.c src/parameters_io.h
        src/diagnostics.h src/diagnostics.c
        src/array.h src/array.c
        src/solver.h src/solver.c
        src/fields.h src/fields.c
        src/equation.h src/equation.c
        src/init.c src/init.h
        src/space_config.c src/space_config.h
        src/distrib.h src/distrib.c
        src/variables.h src/variables.c
        src/utils_tests.h src/utils_tests.c)
target_link_libraries(Alliance_alpha_1_0 "/home/alcauchy/hdf5/hdf5/lib/libhdf5.so")
target_link_libraries(Alliance_alpha_1_0 "/usr/local/lib/libfftw3_mpi.a")
target_link_libraries(Alliance_alpha_1_0 "/usr/local/lib/libfftw3.a")
target_link_libraries(Alliance_alpha_1_0 fftw3_mpi)
target_link_libraries(Alliance_alpha_1_0 fftw3)
target_link_libraries(Alliance_alpha_1_0 m)