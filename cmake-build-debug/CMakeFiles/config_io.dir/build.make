# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.20

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
CMAKE_COMMAND = /Applications/CLion.app/Contents/bin/cmake/mac/bin/cmake

# The command to remove a file.
RM = /Applications/CLion.app/Contents/bin/cmake/mac/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = "/Users/zlw/CLionProjects/MCCPS luban"

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = "/Users/zlw/CLionProjects/MCCPS luban/cmake-build-debug"

# Include any dependencies generated for this target.
include CMakeFiles/config_io.dir/depend.make
# Include the progress variables for this target.
include CMakeFiles/config_io.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/config_io.dir/flags.make

CMakeFiles/config_io.dir/config_io.cpp.o: CMakeFiles/config_io.dir/flags.make
CMakeFiles/config_io.dir/config_io.cpp.o: ../config_io.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir="/Users/zlw/CLionProjects/MCCPS luban/cmake-build-debug/CMakeFiles" --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/config_io.dir/config_io.cpp.o"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/config_io.dir/config_io.cpp.o -c "/Users/zlw/CLionProjects/MCCPS luban/config_io.cpp"

CMakeFiles/config_io.dir/config_io.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/config_io.dir/config_io.cpp.i"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E "/Users/zlw/CLionProjects/MCCPS luban/config_io.cpp" > CMakeFiles/config_io.dir/config_io.cpp.i

CMakeFiles/config_io.dir/config_io.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/config_io.dir/config_io.cpp.s"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S "/Users/zlw/CLionProjects/MCCPS luban/config_io.cpp" -o CMakeFiles/config_io.dir/config_io.cpp.s

# Object files for target config_io
config_io_OBJECTS = \
"CMakeFiles/config_io.dir/config_io.cpp.o"

# External object files for target config_io
config_io_EXTERNAL_OBJECTS =

config_io: CMakeFiles/config_io.dir/config_io.cpp.o
config_io: CMakeFiles/config_io.dir/build.make
config_io: CMakeFiles/config_io.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir="/Users/zlw/CLionProjects/MCCPS luban/cmake-build-debug/CMakeFiles" --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable config_io"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/config_io.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/config_io.dir/build: config_io
.PHONY : CMakeFiles/config_io.dir/build

CMakeFiles/config_io.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/config_io.dir/cmake_clean.cmake
.PHONY : CMakeFiles/config_io.dir/clean

CMakeFiles/config_io.dir/depend:
	cd "/Users/zlw/CLionProjects/MCCPS luban/cmake-build-debug" && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" "/Users/zlw/CLionProjects/MCCPS luban" "/Users/zlw/CLionProjects/MCCPS luban" "/Users/zlw/CLionProjects/MCCPS luban/cmake-build-debug" "/Users/zlw/CLionProjects/MCCPS luban/cmake-build-debug" "/Users/zlw/CLionProjects/MCCPS luban/cmake-build-debug/CMakeFiles/config_io.dir/DependInfo.cmake" --color=$(COLOR)
.PHONY : CMakeFiles/config_io.dir/depend
