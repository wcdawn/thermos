cmake_minimum_required(VERSION 3.14)
project(thermos VERSION 0.1.0 DESCRIPTION "A one-dimensional heat conduction solver." LANGUAGES Fortran)

add_executable(${PROJECT_NAME}.x)

target_sources(${PROJECT_NAME}.x
  PRIVATE
    src/analysis.f90
    src/conductivity_function.f90
    src/constants.f90
    src/exception_handler.f90
    src/fileio.f90
    src/finite_difference.f90
    src/finite_element.f90
    src/geometry.f90
    src/input.f90
    src/kind.f90
    src/linalg.f90
    src/main.f90
    src/output.f90
    src/quadrature.f90
    src/source_function.f90
)

# Note that the CMake language does not allow setting the Fortran_STANDARD from the build instructions
# and it must be specified via the command-line argument.
target_compile_options(${PROJECT_NAME}.x PRIVATE -Wall -Wextra -std=f2018 -pedantic)
target_compile_options(${PROJECT_NAME}.x PRIVATE $<$<CONFIG:Debug>:-O0;-g;-fbacktrace;-ffpe-trap=invalid,overflow,zero;-fcheck=all>)
target_compile_options(${PROJECT_NAME}.x PRIVATE $<$<CONFIG:Release>:-O3>)
target_compile_options(${PROJECT_NAME}.x PRIVATE $<$<CONFIG:>:-O3>) # include default/empty option
