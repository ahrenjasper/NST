# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.14

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list


# Suppress display of executed commands.
$(VERBOSE).SILENT:


# A target that is always out of date.
cmake_force:

.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /blues/gpfs/home/kmoore/miniconda/envs/pacc-env-mini/bin/cmake

# The command to remove a file.
RM = /blues/gpfs/home/kmoore/miniconda/envs/pacc-env-mini/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /lcrc/project/CMRP/pacc/NST

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /lcrc/project/CMRP/pacc/NST/build

# Include any dependencies generated for this target.
include CMakeFiles/nstlibs.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/nstlibs.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/nstlibs.dir/flags.make

CMakeFiles/nstlibs.dir/src/msxfreq.f.o: CMakeFiles/nstlibs.dir/flags.make
CMakeFiles/nstlibs.dir/src/msxfreq.f.o: ../src/msxfreq.f
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/lcrc/project/CMRP/pacc/NST/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building Fortran object CMakeFiles/nstlibs.dir/src/msxfreq.f.o"
	/home/kmoore/miniconda/envs/pacc-env-mini/bin/x86_64-conda_cos6-linux-gnu-gfortran $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -c /lcrc/project/CMRP/pacc/NST/src/msxfreq.f -o CMakeFiles/nstlibs.dir/src/msxfreq.f.o

CMakeFiles/nstlibs.dir/src/msxfreq.f.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing Fortran source to CMakeFiles/nstlibs.dir/src/msxfreq.f.i"
	/home/kmoore/miniconda/envs/pacc-env-mini/bin/x86_64-conda_cos6-linux-gnu-gfortran $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -E /lcrc/project/CMRP/pacc/NST/src/msxfreq.f > CMakeFiles/nstlibs.dir/src/msxfreq.f.i

CMakeFiles/nstlibs.dir/src/msxfreq.f.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling Fortran source to assembly CMakeFiles/nstlibs.dir/src/msxfreq.f.s"
	/home/kmoore/miniconda/envs/pacc-env-mini/bin/x86_64-conda_cos6-linux-gnu-gfortran $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -S /lcrc/project/CMRP/pacc/NST/src/msxfreq.f -o CMakeFiles/nstlibs.dir/src/msxfreq.f.s

CMakeFiles/nstlibs.dir/src/opt.f.o: CMakeFiles/nstlibs.dir/flags.make
CMakeFiles/nstlibs.dir/src/opt.f.o: ../src/opt.f
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/lcrc/project/CMRP/pacc/NST/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building Fortran object CMakeFiles/nstlibs.dir/src/opt.f.o"
	/home/kmoore/miniconda/envs/pacc-env-mini/bin/x86_64-conda_cos6-linux-gnu-gfortran $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -c /lcrc/project/CMRP/pacc/NST/src/opt.f -o CMakeFiles/nstlibs.dir/src/opt.f.o

CMakeFiles/nstlibs.dir/src/opt.f.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing Fortran source to CMakeFiles/nstlibs.dir/src/opt.f.i"
	/home/kmoore/miniconda/envs/pacc-env-mini/bin/x86_64-conda_cos6-linux-gnu-gfortran $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -E /lcrc/project/CMRP/pacc/NST/src/opt.f > CMakeFiles/nstlibs.dir/src/opt.f.i

CMakeFiles/nstlibs.dir/src/opt.f.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling Fortran source to assembly CMakeFiles/nstlibs.dir/src/opt.f.s"
	/home/kmoore/miniconda/envs/pacc-env-mini/bin/x86_64-conda_cos6-linux-gnu-gfortran $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -S /lcrc/project/CMRP/pacc/NST/src/opt.f -o CMakeFiles/nstlibs.dir/src/opt.f.s

CMakeFiles/nstlibs.dir/src/proj.f.o: CMakeFiles/nstlibs.dir/flags.make
CMakeFiles/nstlibs.dir/src/proj.f.o: ../src/proj.f
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/lcrc/project/CMRP/pacc/NST/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building Fortran object CMakeFiles/nstlibs.dir/src/proj.f.o"
	/home/kmoore/miniconda/envs/pacc-env-mini/bin/x86_64-conda_cos6-linux-gnu-gfortran $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -c /lcrc/project/CMRP/pacc/NST/src/proj.f -o CMakeFiles/nstlibs.dir/src/proj.f.o

CMakeFiles/nstlibs.dir/src/proj.f.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing Fortran source to CMakeFiles/nstlibs.dir/src/proj.f.i"
	/home/kmoore/miniconda/envs/pacc-env-mini/bin/x86_64-conda_cos6-linux-gnu-gfortran $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -E /lcrc/project/CMRP/pacc/NST/src/proj.f > CMakeFiles/nstlibs.dir/src/proj.f.i

CMakeFiles/nstlibs.dir/src/proj.f.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling Fortran source to assembly CMakeFiles/nstlibs.dir/src/proj.f.s"
	/home/kmoore/miniconda/envs/pacc-env-mini/bin/x86_64-conda_cos6-linux-gnu-gfortran $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -S /lcrc/project/CMRP/pacc/NST/src/proj.f -o CMakeFiles/nstlibs.dir/src/proj.f.s

CMakeFiles/nstlibs.dir/src/dd.f.o: CMakeFiles/nstlibs.dir/flags.make
CMakeFiles/nstlibs.dir/src/dd.f.o: ../src/dd.f
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/lcrc/project/CMRP/pacc/NST/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Building Fortran object CMakeFiles/nstlibs.dir/src/dd.f.o"
	/home/kmoore/miniconda/envs/pacc-env-mini/bin/x86_64-conda_cos6-linux-gnu-gfortran $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -c /lcrc/project/CMRP/pacc/NST/src/dd.f -o CMakeFiles/nstlibs.dir/src/dd.f.o

CMakeFiles/nstlibs.dir/src/dd.f.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing Fortran source to CMakeFiles/nstlibs.dir/src/dd.f.i"
	/home/kmoore/miniconda/envs/pacc-env-mini/bin/x86_64-conda_cos6-linux-gnu-gfortran $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -E /lcrc/project/CMRP/pacc/NST/src/dd.f > CMakeFiles/nstlibs.dir/src/dd.f.i

CMakeFiles/nstlibs.dir/src/dd.f.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling Fortran source to assembly CMakeFiles/nstlibs.dir/src/dd.f.s"
	/home/kmoore/miniconda/envs/pacc-env-mini/bin/x86_64-conda_cos6-linux-gnu-gfortran $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -S /lcrc/project/CMRP/pacc/NST/src/dd.f -o CMakeFiles/nstlibs.dir/src/dd.f.s

CMakeFiles/nstlibs.dir/src/airy.f.o: CMakeFiles/nstlibs.dir/flags.make
CMakeFiles/nstlibs.dir/src/airy.f.o: ../src/airy.f
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/lcrc/project/CMRP/pacc/NST/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_5) "Building Fortran object CMakeFiles/nstlibs.dir/src/airy.f.o"
	/home/kmoore/miniconda/envs/pacc-env-mini/bin/x86_64-conda_cos6-linux-gnu-gfortran $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -c /lcrc/project/CMRP/pacc/NST/src/airy.f -o CMakeFiles/nstlibs.dir/src/airy.f.o

CMakeFiles/nstlibs.dir/src/airy.f.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing Fortran source to CMakeFiles/nstlibs.dir/src/airy.f.i"
	/home/kmoore/miniconda/envs/pacc-env-mini/bin/x86_64-conda_cos6-linux-gnu-gfortran $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -E /lcrc/project/CMRP/pacc/NST/src/airy.f > CMakeFiles/nstlibs.dir/src/airy.f.i

CMakeFiles/nstlibs.dir/src/airy.f.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling Fortran source to assembly CMakeFiles/nstlibs.dir/src/airy.f.s"
	/home/kmoore/miniconda/envs/pacc-env-mini/bin/x86_64-conda_cos6-linux-gnu-gfortran $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -S /lcrc/project/CMRP/pacc/NST/src/airy.f -o CMakeFiles/nstlibs.dir/src/airy.f.s

# Object files for target nstlibs
nstlibs_OBJECTS = \
"CMakeFiles/nstlibs.dir/src/msxfreq.f.o" \
"CMakeFiles/nstlibs.dir/src/opt.f.o" \
"CMakeFiles/nstlibs.dir/src/proj.f.o" \
"CMakeFiles/nstlibs.dir/src/dd.f.o" \
"CMakeFiles/nstlibs.dir/src/airy.f.o"

# External object files for target nstlibs
nstlibs_EXTERNAL_OBJECTS =

libnstlibs.a: CMakeFiles/nstlibs.dir/src/msxfreq.f.o
libnstlibs.a: CMakeFiles/nstlibs.dir/src/opt.f.o
libnstlibs.a: CMakeFiles/nstlibs.dir/src/proj.f.o
libnstlibs.a: CMakeFiles/nstlibs.dir/src/dd.f.o
libnstlibs.a: CMakeFiles/nstlibs.dir/src/airy.f.o
libnstlibs.a: CMakeFiles/nstlibs.dir/build.make
libnstlibs.a: CMakeFiles/nstlibs.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/lcrc/project/CMRP/pacc/NST/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_6) "Linking Fortran static library libnstlibs.a"
	$(CMAKE_COMMAND) -P CMakeFiles/nstlibs.dir/cmake_clean_target.cmake
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/nstlibs.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/nstlibs.dir/build: libnstlibs.a

.PHONY : CMakeFiles/nstlibs.dir/build

CMakeFiles/nstlibs.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/nstlibs.dir/cmake_clean.cmake
.PHONY : CMakeFiles/nstlibs.dir/clean

CMakeFiles/nstlibs.dir/depend:
	cd /lcrc/project/CMRP/pacc/NST/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /lcrc/project/CMRP/pacc/NST /lcrc/project/CMRP/pacc/NST /lcrc/project/CMRP/pacc/NST/build /lcrc/project/CMRP/pacc/NST/build /lcrc/project/CMRP/pacc/NST/build/CMakeFiles/nstlibs.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/nstlibs.dir/depend
