# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.5

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
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/marsland/Desktop/TIdentityCMake_Dev

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/marsland/Desktop/TIdentityCMake_Dev/build

# Include any dependencies generated for this target.
include test/CMakeFiles/testIdentity.dir/depend.make

# Include the progress variables for this target.
include test/CMakeFiles/testIdentity.dir/progress.make

# Include the compile flags for this target's objects.
include test/CMakeFiles/testIdentity.dir/flags.make

test/CMakeFiles/testIdentity.dir/testIdentity.o: test/CMakeFiles/testIdentity.dir/flags.make
test/CMakeFiles/testIdentity.dir/testIdentity.o: ../test/testIdentity.C
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/marsland/Desktop/TIdentityCMake_Dev/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object test/CMakeFiles/testIdentity.dir/testIdentity.o"
	cd /home/marsland/Desktop/TIdentityCMake_Dev/build/test && /usr/bin/c++   $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/testIdentity.dir/testIdentity.o -c /home/marsland/Desktop/TIdentityCMake_Dev/test/testIdentity.C

test/CMakeFiles/testIdentity.dir/testIdentity.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/testIdentity.dir/testIdentity.i"
	cd /home/marsland/Desktop/TIdentityCMake_Dev/build/test && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/marsland/Desktop/TIdentityCMake_Dev/test/testIdentity.C > CMakeFiles/testIdentity.dir/testIdentity.i

test/CMakeFiles/testIdentity.dir/testIdentity.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/testIdentity.dir/testIdentity.s"
	cd /home/marsland/Desktop/TIdentityCMake_Dev/build/test && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/marsland/Desktop/TIdentityCMake_Dev/test/testIdentity.C -o CMakeFiles/testIdentity.dir/testIdentity.s

test/CMakeFiles/testIdentity.dir/testIdentity.o.requires:

.PHONY : test/CMakeFiles/testIdentity.dir/testIdentity.o.requires

test/CMakeFiles/testIdentity.dir/testIdentity.o.provides: test/CMakeFiles/testIdentity.dir/testIdentity.o.requires
	$(MAKE) -f test/CMakeFiles/testIdentity.dir/build.make test/CMakeFiles/testIdentity.dir/testIdentity.o.provides.build
.PHONY : test/CMakeFiles/testIdentity.dir/testIdentity.o.provides

test/CMakeFiles/testIdentity.dir/testIdentity.o.provides.build: test/CMakeFiles/testIdentity.dir/testIdentity.o


# Object files for target testIdentity
testIdentity_OBJECTS = \
"CMakeFiles/testIdentity.dir/testIdentity.o"

# External object files for target testIdentity
testIdentity_EXTERNAL_OBJECTS =

../bin/testIdentity: test/CMakeFiles/testIdentity.dir/testIdentity.o
../bin/testIdentity: test/CMakeFiles/testIdentity.dir/build.make
../bin/testIdentity: ../lib/libTIdentity.so
../bin/testIdentity: /home/marsland/alice/sw/ubuntu1604_x86-64/ROOT/v5-34-30-alice10-1/lib/libCore.so
../bin/testIdentity: /home/marsland/alice/sw/ubuntu1604_x86-64/ROOT/v5-34-30-alice10-1/lib/libCint.so
../bin/testIdentity: /home/marsland/alice/sw/ubuntu1604_x86-64/ROOT/v5-34-30-alice10-1/lib/libMathMore.so
../bin/testIdentity: /home/marsland/alice/sw/ubuntu1604_x86-64/ROOT/v5-34-30-alice10-1/lib/libMathCore.so
../bin/testIdentity: /home/marsland/alice/sw/ubuntu1604_x86-64/ROOT/v5-34-30-alice10-1/lib/libRIO.so
../bin/testIdentity: /home/marsland/alice/sw/ubuntu1604_x86-64/ROOT/v5-34-30-alice10-1/lib/libHist.so
../bin/testIdentity: /home/marsland/alice/sw/ubuntu1604_x86-64/ROOT/v5-34-30-alice10-1/lib/libTree.so
../bin/testIdentity: /home/marsland/alice/sw/ubuntu1604_x86-64/ROOT/v5-34-30-alice10-1/lib/libNet.so
../bin/testIdentity: /home/marsland/alice/sw/ubuntu1604_x86-64/ROOT/v5-34-30-alice10-1/lib/libMatrix.so
../bin/testIdentity: test/CMakeFiles/testIdentity.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/marsland/Desktop/TIdentityCMake_Dev/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable ../../bin/testIdentity"
	cd /home/marsland/Desktop/TIdentityCMake_Dev/build/test && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/testIdentity.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
test/CMakeFiles/testIdentity.dir/build: ../bin/testIdentity

.PHONY : test/CMakeFiles/testIdentity.dir/build

test/CMakeFiles/testIdentity.dir/requires: test/CMakeFiles/testIdentity.dir/testIdentity.o.requires

.PHONY : test/CMakeFiles/testIdentity.dir/requires

test/CMakeFiles/testIdentity.dir/clean:
	cd /home/marsland/Desktop/TIdentityCMake_Dev/build/test && $(CMAKE_COMMAND) -P CMakeFiles/testIdentity.dir/cmake_clean.cmake
.PHONY : test/CMakeFiles/testIdentity.dir/clean

test/CMakeFiles/testIdentity.dir/depend:
	cd /home/marsland/Desktop/TIdentityCMake_Dev/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/marsland/Desktop/TIdentityCMake_Dev /home/marsland/Desktop/TIdentityCMake_Dev/test /home/marsland/Desktop/TIdentityCMake_Dev/build /home/marsland/Desktop/TIdentityCMake_Dev/build/test /home/marsland/Desktop/TIdentityCMake_Dev/build/test/CMakeFiles/testIdentity.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : test/CMakeFiles/testIdentity.dir/depend
