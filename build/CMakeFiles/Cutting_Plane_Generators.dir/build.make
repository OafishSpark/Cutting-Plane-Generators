# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.22

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
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /mnt/d/Projects/Cutting-Plane-Generators

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /mnt/d/Projects/Cutting-Plane-Generators/build

# Include any dependencies generated for this target.
include CMakeFiles/Cutting_Plane_Generators.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include CMakeFiles/Cutting_Plane_Generators.dir/compiler_depend.make

# Include the progress variables for this target.
include CMakeFiles/Cutting_Plane_Generators.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/Cutting_Plane_Generators.dir/flags.make

CMakeFiles/Cutting_Plane_Generators.dir/main.cpp.o: CMakeFiles/Cutting_Plane_Generators.dir/flags.make
CMakeFiles/Cutting_Plane_Generators.dir/main.cpp.o: ../main.cpp
CMakeFiles/Cutting_Plane_Generators.dir/main.cpp.o: CMakeFiles/Cutting_Plane_Generators.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/mnt/d/Projects/Cutting-Plane-Generators/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/Cutting_Plane_Generators.dir/main.cpp.o"
	/usr/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/Cutting_Plane_Generators.dir/main.cpp.o -MF CMakeFiles/Cutting_Plane_Generators.dir/main.cpp.o.d -o CMakeFiles/Cutting_Plane_Generators.dir/main.cpp.o -c /mnt/d/Projects/Cutting-Plane-Generators/main.cpp

CMakeFiles/Cutting_Plane_Generators.dir/main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/Cutting_Plane_Generators.dir/main.cpp.i"
	/usr/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /mnt/d/Projects/Cutting-Plane-Generators/main.cpp > CMakeFiles/Cutting_Plane_Generators.dir/main.cpp.i

CMakeFiles/Cutting_Plane_Generators.dir/main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/Cutting_Plane_Generators.dir/main.cpp.s"
	/usr/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /mnt/d/Projects/Cutting-Plane-Generators/main.cpp -o CMakeFiles/Cutting_Plane_Generators.dir/main.cpp.s

CMakeFiles/Cutting_Plane_Generators.dir/src/cuts/cutter.cpp.o: CMakeFiles/Cutting_Plane_Generators.dir/flags.make
CMakeFiles/Cutting_Plane_Generators.dir/src/cuts/cutter.cpp.o: ../src/cuts/cutter.cpp
CMakeFiles/Cutting_Plane_Generators.dir/src/cuts/cutter.cpp.o: CMakeFiles/Cutting_Plane_Generators.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/mnt/d/Projects/Cutting-Plane-Generators/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object CMakeFiles/Cutting_Plane_Generators.dir/src/cuts/cutter.cpp.o"
	/usr/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/Cutting_Plane_Generators.dir/src/cuts/cutter.cpp.o -MF CMakeFiles/Cutting_Plane_Generators.dir/src/cuts/cutter.cpp.o.d -o CMakeFiles/Cutting_Plane_Generators.dir/src/cuts/cutter.cpp.o -c /mnt/d/Projects/Cutting-Plane-Generators/src/cuts/cutter.cpp

CMakeFiles/Cutting_Plane_Generators.dir/src/cuts/cutter.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/Cutting_Plane_Generators.dir/src/cuts/cutter.cpp.i"
	/usr/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /mnt/d/Projects/Cutting-Plane-Generators/src/cuts/cutter.cpp > CMakeFiles/Cutting_Plane_Generators.dir/src/cuts/cutter.cpp.i

CMakeFiles/Cutting_Plane_Generators.dir/src/cuts/cutter.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/Cutting_Plane_Generators.dir/src/cuts/cutter.cpp.s"
	/usr/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /mnt/d/Projects/Cutting-Plane-Generators/src/cuts/cutter.cpp -o CMakeFiles/Cutting_Plane_Generators.dir/src/cuts/cutter.cpp.s

CMakeFiles/Cutting_Plane_Generators.dir/src/cuts/gmi.cpp.o: CMakeFiles/Cutting_Plane_Generators.dir/flags.make
CMakeFiles/Cutting_Plane_Generators.dir/src/cuts/gmi.cpp.o: ../src/cuts/gmi.cpp
CMakeFiles/Cutting_Plane_Generators.dir/src/cuts/gmi.cpp.o: CMakeFiles/Cutting_Plane_Generators.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/mnt/d/Projects/Cutting-Plane-Generators/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building CXX object CMakeFiles/Cutting_Plane_Generators.dir/src/cuts/gmi.cpp.o"
	/usr/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/Cutting_Plane_Generators.dir/src/cuts/gmi.cpp.o -MF CMakeFiles/Cutting_Plane_Generators.dir/src/cuts/gmi.cpp.o.d -o CMakeFiles/Cutting_Plane_Generators.dir/src/cuts/gmi.cpp.o -c /mnt/d/Projects/Cutting-Plane-Generators/src/cuts/gmi.cpp

CMakeFiles/Cutting_Plane_Generators.dir/src/cuts/gmi.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/Cutting_Plane_Generators.dir/src/cuts/gmi.cpp.i"
	/usr/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /mnt/d/Projects/Cutting-Plane-Generators/src/cuts/gmi.cpp > CMakeFiles/Cutting_Plane_Generators.dir/src/cuts/gmi.cpp.i

CMakeFiles/Cutting_Plane_Generators.dir/src/cuts/gmi.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/Cutting_Plane_Generators.dir/src/cuts/gmi.cpp.s"
	/usr/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /mnt/d/Projects/Cutting-Plane-Generators/src/cuts/gmi.cpp -o CMakeFiles/Cutting_Plane_Generators.dir/src/cuts/gmi.cpp.s

CMakeFiles/Cutting_Plane_Generators.dir/src/parser/parser.cpp.o: CMakeFiles/Cutting_Plane_Generators.dir/flags.make
CMakeFiles/Cutting_Plane_Generators.dir/src/parser/parser.cpp.o: ../src/parser/parser.cpp
CMakeFiles/Cutting_Plane_Generators.dir/src/parser/parser.cpp.o: CMakeFiles/Cutting_Plane_Generators.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/mnt/d/Projects/Cutting-Plane-Generators/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Building CXX object CMakeFiles/Cutting_Plane_Generators.dir/src/parser/parser.cpp.o"
	/usr/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/Cutting_Plane_Generators.dir/src/parser/parser.cpp.o -MF CMakeFiles/Cutting_Plane_Generators.dir/src/parser/parser.cpp.o.d -o CMakeFiles/Cutting_Plane_Generators.dir/src/parser/parser.cpp.o -c /mnt/d/Projects/Cutting-Plane-Generators/src/parser/parser.cpp

CMakeFiles/Cutting_Plane_Generators.dir/src/parser/parser.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/Cutting_Plane_Generators.dir/src/parser/parser.cpp.i"
	/usr/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /mnt/d/Projects/Cutting-Plane-Generators/src/parser/parser.cpp > CMakeFiles/Cutting_Plane_Generators.dir/src/parser/parser.cpp.i

CMakeFiles/Cutting_Plane_Generators.dir/src/parser/parser.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/Cutting_Plane_Generators.dir/src/parser/parser.cpp.s"
	/usr/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /mnt/d/Projects/Cutting-Plane-Generators/src/parser/parser.cpp -o CMakeFiles/Cutting_Plane_Generators.dir/src/parser/parser.cpp.s

CMakeFiles/Cutting_Plane_Generators.dir/src/linalg/linalg.cpp.o: CMakeFiles/Cutting_Plane_Generators.dir/flags.make
CMakeFiles/Cutting_Plane_Generators.dir/src/linalg/linalg.cpp.o: ../src/linalg/linalg.cpp
CMakeFiles/Cutting_Plane_Generators.dir/src/linalg/linalg.cpp.o: CMakeFiles/Cutting_Plane_Generators.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/mnt/d/Projects/Cutting-Plane-Generators/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_5) "Building CXX object CMakeFiles/Cutting_Plane_Generators.dir/src/linalg/linalg.cpp.o"
	/usr/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/Cutting_Plane_Generators.dir/src/linalg/linalg.cpp.o -MF CMakeFiles/Cutting_Plane_Generators.dir/src/linalg/linalg.cpp.o.d -o CMakeFiles/Cutting_Plane_Generators.dir/src/linalg/linalg.cpp.o -c /mnt/d/Projects/Cutting-Plane-Generators/src/linalg/linalg.cpp

CMakeFiles/Cutting_Plane_Generators.dir/src/linalg/linalg.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/Cutting_Plane_Generators.dir/src/linalg/linalg.cpp.i"
	/usr/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /mnt/d/Projects/Cutting-Plane-Generators/src/linalg/linalg.cpp > CMakeFiles/Cutting_Plane_Generators.dir/src/linalg/linalg.cpp.i

CMakeFiles/Cutting_Plane_Generators.dir/src/linalg/linalg.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/Cutting_Plane_Generators.dir/src/linalg/linalg.cpp.s"
	/usr/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /mnt/d/Projects/Cutting-Plane-Generators/src/linalg/linalg.cpp -o CMakeFiles/Cutting_Plane_Generators.dir/src/linalg/linalg.cpp.s

# Object files for target Cutting_Plane_Generators
Cutting_Plane_Generators_OBJECTS = \
"CMakeFiles/Cutting_Plane_Generators.dir/main.cpp.o" \
"CMakeFiles/Cutting_Plane_Generators.dir/src/cuts/cutter.cpp.o" \
"CMakeFiles/Cutting_Plane_Generators.dir/src/cuts/gmi.cpp.o" \
"CMakeFiles/Cutting_Plane_Generators.dir/src/parser/parser.cpp.o" \
"CMakeFiles/Cutting_Plane_Generators.dir/src/linalg/linalg.cpp.o"

# External object files for target Cutting_Plane_Generators
Cutting_Plane_Generators_EXTERNAL_OBJECTS =

Cutting_Plane_Generators: CMakeFiles/Cutting_Plane_Generators.dir/main.cpp.o
Cutting_Plane_Generators: CMakeFiles/Cutting_Plane_Generators.dir/src/cuts/cutter.cpp.o
Cutting_Plane_Generators: CMakeFiles/Cutting_Plane_Generators.dir/src/cuts/gmi.cpp.o
Cutting_Plane_Generators: CMakeFiles/Cutting_Plane_Generators.dir/src/parser/parser.cpp.o
Cutting_Plane_Generators: CMakeFiles/Cutting_Plane_Generators.dir/src/linalg/linalg.cpp.o
Cutting_Plane_Generators: CMakeFiles/Cutting_Plane_Generators.dir/build.make
Cutting_Plane_Generators: CMakeFiles/Cutting_Plane_Generators.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/mnt/d/Projects/Cutting-Plane-Generators/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_6) "Linking CXX executable Cutting_Plane_Generators"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/Cutting_Plane_Generators.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/Cutting_Plane_Generators.dir/build: Cutting_Plane_Generators
.PHONY : CMakeFiles/Cutting_Plane_Generators.dir/build

CMakeFiles/Cutting_Plane_Generators.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/Cutting_Plane_Generators.dir/cmake_clean.cmake
.PHONY : CMakeFiles/Cutting_Plane_Generators.dir/clean

CMakeFiles/Cutting_Plane_Generators.dir/depend:
	cd /mnt/d/Projects/Cutting-Plane-Generators/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /mnt/d/Projects/Cutting-Plane-Generators /mnt/d/Projects/Cutting-Plane-Generators /mnt/d/Projects/Cutting-Plane-Generators/build /mnt/d/Projects/Cutting-Plane-Generators/build /mnt/d/Projects/Cutting-Plane-Generators/build/CMakeFiles/Cutting_Plane_Generators.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/Cutting_Plane_Generators.dir/depend
