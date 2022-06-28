cmake_minimum_required(VERSION 3.1)

project(oscillation_omega)
enable_language(Fortran)

function(add_fortran_library LIB)
    add_library(${LIB} ${ARGN})

    # set module path to LIB_DIR/mod
    get_target_property(LIB_DIR ${LIB} BINARY_DIR)
    set_target_properties(${LIB} PROPERTIES Fortran_MODULE_DIRECTORY ${LIB_DIR}/mod)

    # making LIB_DIR/mod available for libraries linking LIB 
    target_include_directories(${LIB} INTERFACE ${LIB_DIR}/mod)
endfunction(add_fortran_library)

if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
    set(dialect "-ffree-form -std=f2008 -fimplicit-none")
    set(bounds "-fbounds-check")
endif()
if(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
    set(dialect "-stand f08 -free -implicitnone")
    set(bounds "-check bounds")
endif()

set(CMAKE_Fortran_FLAGS_DEBUG   "-Wall -O0 -g3 -fbounds-check")
set(CMAKE_Fortran_FLAGS_RELEASE "-O3")

set(TARGET oscillation_omega)
set(SOURCE_FILES src/oscillation_omega.f90)

set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/library)

add_fortran_library(discontinuities SHARED src/discontinuities.f90)
add_fortran_library(legendre_quadrature SHARED src/legendre_quadrature.f90)
add_fortran_library(read_model SHARED src/read_model.f90)
add_fortran_library(ylm SHARED src/ylm.f90)
add_fortran_library(gravitational_potential SHARED src/compute_gravitational_potential.f90)
add_fortran_library(moment_of_inertia SHARED src/compute_moment_of_inertia.f90)

target_link_libraries(gravitational_potential discontinuities ylm read_model)
target_link_libraries(moment_of_inertia gravitational_potential legendre_quadrature)

add_executable(${TARGET} ${SOURCE_FILES})
target_link_libraries(${TARGET} gravitational_potential moment_of_inertia)

set_target_properties(${TARGET} PROPERTIES LINKER_LANGUAGE Fortran)
