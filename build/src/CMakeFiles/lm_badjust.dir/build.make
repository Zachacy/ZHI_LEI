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
CMAKE_SOURCE_DIR = /home/wang/final_end

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/wang/final_end/build

# Include any dependencies generated for this target.
include src/CMakeFiles/lm_badjust.dir/depend.make

# Include the progress variables for this target.
include src/CMakeFiles/lm_badjust.dir/progress.make

# Include the compile flags for this target's objects.
include src/CMakeFiles/lm_badjust.dir/flags.make

src/CMakeFiles/lm_badjust.dir/LM_badjust.cpp.o: src/CMakeFiles/lm_badjust.dir/flags.make
src/CMakeFiles/lm_badjust.dir/LM_badjust.cpp.o: ../src/LM_badjust.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/wang/final_end/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object src/CMakeFiles/lm_badjust.dir/LM_badjust.cpp.o"
	cd /home/wang/final_end/build/src && g++   $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/lm_badjust.dir/LM_badjust.cpp.o -c /home/wang/final_end/src/LM_badjust.cpp

src/CMakeFiles/lm_badjust.dir/LM_badjust.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/lm_badjust.dir/LM_badjust.cpp.i"
	cd /home/wang/final_end/build/src && g++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/wang/final_end/src/LM_badjust.cpp > CMakeFiles/lm_badjust.dir/LM_badjust.cpp.i

src/CMakeFiles/lm_badjust.dir/LM_badjust.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/lm_badjust.dir/LM_badjust.cpp.s"
	cd /home/wang/final_end/build/src && g++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/wang/final_end/src/LM_badjust.cpp -o CMakeFiles/lm_badjust.dir/LM_badjust.cpp.s

src/CMakeFiles/lm_badjust.dir/LM_badjust.cpp.o.requires:

.PHONY : src/CMakeFiles/lm_badjust.dir/LM_badjust.cpp.o.requires

src/CMakeFiles/lm_badjust.dir/LM_badjust.cpp.o.provides: src/CMakeFiles/lm_badjust.dir/LM_badjust.cpp.o.requires
	$(MAKE) -f src/CMakeFiles/lm_badjust.dir/build.make src/CMakeFiles/lm_badjust.dir/LM_badjust.cpp.o.provides.build
.PHONY : src/CMakeFiles/lm_badjust.dir/LM_badjust.cpp.o.provides

src/CMakeFiles/lm_badjust.dir/LM_badjust.cpp.o.provides.build: src/CMakeFiles/lm_badjust.dir/LM_badjust.cpp.o


# Object files for target lm_badjust
lm_badjust_OBJECTS = \
"CMakeFiles/lm_badjust.dir/LM_badjust.cpp.o"

# External object files for target lm_badjust
lm_badjust_EXTERNAL_OBJECTS =

../lib/liblm_badjust.a: src/CMakeFiles/lm_badjust.dir/LM_badjust.cpp.o
../lib/liblm_badjust.a: src/CMakeFiles/lm_badjust.dir/build.make
../lib/liblm_badjust.a: src/CMakeFiles/lm_badjust.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/wang/final_end/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX static library ../../lib/liblm_badjust.a"
	cd /home/wang/final_end/build/src && $(CMAKE_COMMAND) -P CMakeFiles/lm_badjust.dir/cmake_clean_target.cmake
	cd /home/wang/final_end/build/src && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/lm_badjust.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
src/CMakeFiles/lm_badjust.dir/build: ../lib/liblm_badjust.a

.PHONY : src/CMakeFiles/lm_badjust.dir/build

src/CMakeFiles/lm_badjust.dir/requires: src/CMakeFiles/lm_badjust.dir/LM_badjust.cpp.o.requires

.PHONY : src/CMakeFiles/lm_badjust.dir/requires

src/CMakeFiles/lm_badjust.dir/clean:
	cd /home/wang/final_end/build/src && $(CMAKE_COMMAND) -P CMakeFiles/lm_badjust.dir/cmake_clean.cmake
.PHONY : src/CMakeFiles/lm_badjust.dir/clean

src/CMakeFiles/lm_badjust.dir/depend:
	cd /home/wang/final_end/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/wang/final_end /home/wang/final_end/src /home/wang/final_end/build /home/wang/final_end/build/src /home/wang/final_end/build/src/CMakeFiles/lm_badjust.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : src/CMakeFiles/lm_badjust.dir/depend
