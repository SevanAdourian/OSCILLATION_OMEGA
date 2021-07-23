# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.19

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Disable VCS-based implicit rules.
% : %,v


# Disable VCS-based implicit rules.
% : RCS/%


# Disable VCS-based implicit rules.
% : RCS/%,v


# Disable VCS-based implicit rules.
% : SCCS/s.%


# Disable VCS-based implicit rules.
% : s.%


.SUFFIXES: .hpux_make_needs_suffix_list


# Command-line flag to silence nested $(MAKE).
$(VERBOSE)MAKESILENT = -s

#Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:

.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/local/Cellar/cmake/3.19.3/bin/cmake

# The command to remove a file.
RM = /usr/local/Cellar/cmake/3.19.3/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /Users/sevan/Documents/PhD/Research/OMEGA_PHD/scripts/FORTRAN/OMEGA_OSCILLATION

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/sevan/Documents/PhD/Research/OMEGA_PHD/scripts/FORTRAN/OMEGA_OSCILLATION/build

# Include any dependencies generated for this target.
include CMakeFiles/legendre_quadrature.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/legendre_quadrature.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/legendre_quadrature.dir/flags.make

CMakeFiles/legendre_quadrature.dir/src/legendre_quadrature.f90.o: CMakeFiles/legendre_quadrature.dir/flags.make
CMakeFiles/legendre_quadrature.dir/src/legendre_quadrature.f90.o: ../src/legendre_quadrature.f90
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/sevan/Documents/PhD/Research/OMEGA_PHD/scripts/FORTRAN/OMEGA_OSCILLATION/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building Fortran object CMakeFiles/legendre_quadrature.dir/src/legendre_quadrature.f90.o"
	/usr/local/bin/gfortran $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -c /Users/sevan/Documents/PhD/Research/OMEGA_PHD/scripts/FORTRAN/OMEGA_OSCILLATION/src/legendre_quadrature.f90 -o CMakeFiles/legendre_quadrature.dir/src/legendre_quadrature.f90.o

CMakeFiles/legendre_quadrature.dir/src/legendre_quadrature.f90.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing Fortran source to CMakeFiles/legendre_quadrature.dir/src/legendre_quadrature.f90.i"
	/usr/local/bin/gfortran $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -E /Users/sevan/Documents/PhD/Research/OMEGA_PHD/scripts/FORTRAN/OMEGA_OSCILLATION/src/legendre_quadrature.f90 > CMakeFiles/legendre_quadrature.dir/src/legendre_quadrature.f90.i

CMakeFiles/legendre_quadrature.dir/src/legendre_quadrature.f90.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling Fortran source to assembly CMakeFiles/legendre_quadrature.dir/src/legendre_quadrature.f90.s"
	/usr/local/bin/gfortran $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -S /Users/sevan/Documents/PhD/Research/OMEGA_PHD/scripts/FORTRAN/OMEGA_OSCILLATION/src/legendre_quadrature.f90 -o CMakeFiles/legendre_quadrature.dir/src/legendre_quadrature.f90.s

# Object files for target legendre_quadrature
legendre_quadrature_OBJECTS = \
"CMakeFiles/legendre_quadrature.dir/src/legendre_quadrature.f90.o"

# External object files for target legendre_quadrature
legendre_quadrature_EXTERNAL_OBJECTS =

library/liblegendre_quadrature.dylib: CMakeFiles/legendre_quadrature.dir/src/legendre_quadrature.f90.o
library/liblegendre_quadrature.dylib: CMakeFiles/legendre_quadrature.dir/build.make
library/liblegendre_quadrature.dylib: CMakeFiles/legendre_quadrature.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/Users/sevan/Documents/PhD/Research/OMEGA_PHD/scripts/FORTRAN/OMEGA_OSCILLATION/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking Fortran shared library library/liblegendre_quadrature.dylib"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/legendre_quadrature.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/legendre_quadrature.dir/build: library/liblegendre_quadrature.dylib

.PHONY : CMakeFiles/legendre_quadrature.dir/build

CMakeFiles/legendre_quadrature.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/legendre_quadrature.dir/cmake_clean.cmake
.PHONY : CMakeFiles/legendre_quadrature.dir/clean

CMakeFiles/legendre_quadrature.dir/depend:
	cd /Users/sevan/Documents/PhD/Research/OMEGA_PHD/scripts/FORTRAN/OMEGA_OSCILLATION/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/sevan/Documents/PhD/Research/OMEGA_PHD/scripts/FORTRAN/OMEGA_OSCILLATION /Users/sevan/Documents/PhD/Research/OMEGA_PHD/scripts/FORTRAN/OMEGA_OSCILLATION /Users/sevan/Documents/PhD/Research/OMEGA_PHD/scripts/FORTRAN/OMEGA_OSCILLATION/build /Users/sevan/Documents/PhD/Research/OMEGA_PHD/scripts/FORTRAN/OMEGA_OSCILLATION/build /Users/sevan/Documents/PhD/Research/OMEGA_PHD/scripts/FORTRAN/OMEGA_OSCILLATION/build/CMakeFiles/legendre_quadrature.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/legendre_quadrature.dir/depend
