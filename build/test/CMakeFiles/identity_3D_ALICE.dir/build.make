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
include test/CMakeFiles/identity_3D_ALICE.dir/depend.make

# Include the progress variables for this target.
include test/CMakeFiles/identity_3D_ALICE.dir/progress.make

# Include the compile flags for this target's objects.
include test/CMakeFiles/identity_3D_ALICE.dir/flags.make

test/CMakeFiles/identity_3D_ALICE.dir/identity_3D_ALICE.o: test/CMakeFiles/identity_3D_ALICE.dir/flags.make
test/CMakeFiles/identity_3D_ALICE.dir/identity_3D_ALICE.o: ../test/identity_3D_ALICE.C
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/marsland/Desktop/TIdentityCMake_Dev/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object test/CMakeFiles/identity_3D_ALICE.dir/identity_3D_ALICE.o"
	cd /home/marsland/Desktop/TIdentityCMake_Dev/build/test && /usr/bin/c++   $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/identity_3D_ALICE.dir/identity_3D_ALICE.o -c /home/marsland/Desktop/TIdentityCMake_Dev/test/identity_3D_ALICE.C

test/CMakeFiles/identity_3D_ALICE.dir/identity_3D_ALICE.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/identity_3D_ALICE.dir/identity_3D_ALICE.i"
	cd /home/marsland/Desktop/TIdentityCMake_Dev/build/test && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/marsland/Desktop/TIdentityCMake_Dev/test/identity_3D_ALICE.C > CMakeFiles/identity_3D_ALICE.dir/identity_3D_ALICE.i

test/CMakeFiles/identity_3D_ALICE.dir/identity_3D_ALICE.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/identity_3D_ALICE.dir/identity_3D_ALICE.s"
	cd /home/marsland/Desktop/TIdentityCMake_Dev/build/test && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/marsland/Desktop/TIdentityCMake_Dev/test/identity_3D_ALICE.C -o CMakeFiles/identity_3D_ALICE.dir/identity_3D_ALICE.s

test/CMakeFiles/identity_3D_ALICE.dir/identity_3D_ALICE.o.requires:

.PHONY : test/CMakeFiles/identity_3D_ALICE.dir/identity_3D_ALICE.o.requires

test/CMakeFiles/identity_3D_ALICE.dir/identity_3D_ALICE.o.provides: test/CMakeFiles/identity_3D_ALICE.dir/identity_3D_ALICE.o.requires
	$(MAKE) -f test/CMakeFiles/identity_3D_ALICE.dir/build.make test/CMakeFiles/identity_3D_ALICE.dir/identity_3D_ALICE.o.provides.build
.PHONY : test/CMakeFiles/identity_3D_ALICE.dir/identity_3D_ALICE.o.provides

test/CMakeFiles/identity_3D_ALICE.dir/identity_3D_ALICE.o.provides.build: test/CMakeFiles/identity_3D_ALICE.dir/identity_3D_ALICE.o


# Object files for target identity_3D_ALICE
identity_3D_ALICE_OBJECTS = \
"CMakeFiles/identity_3D_ALICE.dir/identity_3D_ALICE.o"

# External object files for target identity_3D_ALICE
identity_3D_ALICE_EXTERNAL_OBJECTS =

../bin/identity_3D_ALICE: test/CMakeFiles/identity_3D_ALICE.dir/identity_3D_ALICE.o
../bin/identity_3D_ALICE: test/CMakeFiles/identity_3D_ALICE.dir/build.make
../bin/identity_3D_ALICE: ../lib/libTIdentity.so
../bin/identity_3D_ALICE: /home/marsland/alice/sw/ubuntu1604_x86-64/ROOT/v5-34-30-alice10-1/lib/libCore.so
../bin/identity_3D_ALICE: /home/marsland/alice/sw/ubuntu1604_x86-64/ROOT/v5-34-30-alice10-1/lib/libCint.so
../bin/identity_3D_ALICE: /home/marsland/alice/sw/ubuntu1604_x86-64/ROOT/v5-34-30-alice10-1/lib/libMathMore.so
../bin/identity_3D_ALICE: /home/marsland/alice/sw/ubuntu1604_x86-64/ROOT/v5-34-30-alice10-1/lib/libMathCore.so
../bin/identity_3D_ALICE: /home/marsland/alice/sw/ubuntu1604_x86-64/ROOT/v5-34-30-alice10-1/lib/libRIO.so
../bin/identity_3D_ALICE: /home/marsland/alice/sw/ubuntu1604_x86-64/ROOT/v5-34-30-alice10-1/lib/libHist.so
../bin/identity_3D_ALICE: /home/marsland/alice/sw/ubuntu1604_x86-64/ROOT/v5-34-30-alice10-1/lib/libTree.so
../bin/identity_3D_ALICE: /home/marsland/alice/sw/ubuntu1604_x86-64/ROOT/v5-34-30-alice10-1/lib/libNet.so
../bin/identity_3D_ALICE: /home/marsland/alice/sw/ubuntu1604_x86-64/ROOT/v5-34-30-alice10-1/lib/libMatrix.so
../bin/identity_3D_ALICE: test/CMakeFiles/identity_3D_ALICE.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/marsland/Desktop/TIdentityCMake_Dev/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable ../../bin/identity_3D_ALICE"
	cd /home/marsland/Desktop/TIdentityCMake_Dev/build/test && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/identity_3D_ALICE.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
test/CMakeFiles/identity_3D_ALICE.dir/build: ../bin/identity_3D_ALICE

.PHONY : test/CMakeFiles/identity_3D_ALICE.dir/build

test/CMakeFiles/identity_3D_ALICE.dir/requires: test/CMakeFiles/identity_3D_ALICE.dir/identity_3D_ALICE.o.requires

.PHONY : test/CMakeFiles/identity_3D_ALICE.dir/requires

test/CMakeFiles/identity_3D_ALICE.dir/clean:
	cd /home/marsland/Desktop/TIdentityCMake_Dev/build/test && $(CMAKE_COMMAND) -P CMakeFiles/identity_3D_ALICE.dir/cmake_clean.cmake
.PHONY : test/CMakeFiles/identity_3D_ALICE.dir/clean

test/CMakeFiles/identity_3D_ALICE.dir/depend:
	cd /home/marsland/Desktop/TIdentityCMake_Dev/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/marsland/Desktop/TIdentityCMake_Dev /home/marsland/Desktop/TIdentityCMake_Dev/test /home/marsland/Desktop/TIdentityCMake_Dev/build /home/marsland/Desktop/TIdentityCMake_Dev/build/test /home/marsland/Desktop/TIdentityCMake_Dev/build/test/CMakeFiles/identity_3D_ALICE.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : test/CMakeFiles/identity_3D_ALICE.dir/depend

