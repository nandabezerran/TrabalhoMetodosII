# CMAKE generated file: DO NOT EDIT!
# Generated by "MinGW Makefiles" Generator, CMake Version 3.9

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

SHELL = cmd.exe

# The CMake executable.
CMAKE_COMMAND = "C:\Users\fbeze\AppData\Roaming\JetBrains\CLion 2017.3.4\bin\cmake\bin\cmake.exe"

# The command to remove a file.
RM = "C:\Users\fbeze\AppData\Roaming\JetBrains\CLion 2017.3.4\bin\cmake\bin\cmake.exe" -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = C:\Users\fbeze\Desktop\Faculdade\Metodos_Numericos_II\TrabalhoMetodosII\Trabalho_04

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = C:\Users\fbeze\Desktop\Faculdade\Metodos_Numericos_II\TrabalhoMetodosII\Trabalho_04\cmake-build-debug

# Include any dependencies generated for this target.
include CMakeFiles/Trabalho_04.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/Trabalho_04.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/Trabalho_04.dir/flags.make

CMakeFiles/Trabalho_04.dir/main.c.obj: CMakeFiles/Trabalho_04.dir/flags.make
CMakeFiles/Trabalho_04.dir/main.c.obj: ../main.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=C:\Users\fbeze\Desktop\Faculdade\Metodos_Numericos_II\TrabalhoMetodosII\Trabalho_04\cmake-build-debug\CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building C object CMakeFiles/Trabalho_04.dir/main.c.obj"
	C:\PROGRA~2\MINGW-~1\I686-6~1.0-P\mingw32\bin\gcc.exe $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles\Trabalho_04.dir\main.c.obj   -c C:\Users\fbeze\Desktop\Faculdade\Metodos_Numericos_II\TrabalhoMetodosII\Trabalho_04\main.c

CMakeFiles/Trabalho_04.dir/main.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/Trabalho_04.dir/main.c.i"
	C:\PROGRA~2\MINGW-~1\I686-6~1.0-P\mingw32\bin\gcc.exe $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E C:\Users\fbeze\Desktop\Faculdade\Metodos_Numericos_II\TrabalhoMetodosII\Trabalho_04\main.c > CMakeFiles\Trabalho_04.dir\main.c.i

CMakeFiles/Trabalho_04.dir/main.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/Trabalho_04.dir/main.c.s"
	C:\PROGRA~2\MINGW-~1\I686-6~1.0-P\mingw32\bin\gcc.exe $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S C:\Users\fbeze\Desktop\Faculdade\Metodos_Numericos_II\TrabalhoMetodosII\Trabalho_04\main.c -o CMakeFiles\Trabalho_04.dir\main.c.s

CMakeFiles/Trabalho_04.dir/main.c.obj.requires:

.PHONY : CMakeFiles/Trabalho_04.dir/main.c.obj.requires

CMakeFiles/Trabalho_04.dir/main.c.obj.provides: CMakeFiles/Trabalho_04.dir/main.c.obj.requires
	$(MAKE) -f CMakeFiles\Trabalho_04.dir\build.make CMakeFiles/Trabalho_04.dir/main.c.obj.provides.build
.PHONY : CMakeFiles/Trabalho_04.dir/main.c.obj.provides

CMakeFiles/Trabalho_04.dir/main.c.obj.provides.build: CMakeFiles/Trabalho_04.dir/main.c.obj


# Object files for target Trabalho_04
Trabalho_04_OBJECTS = \
"CMakeFiles/Trabalho_04.dir/main.c.obj"

# External object files for target Trabalho_04
Trabalho_04_EXTERNAL_OBJECTS =

Trabalho_04.exe: CMakeFiles/Trabalho_04.dir/main.c.obj
Trabalho_04.exe: CMakeFiles/Trabalho_04.dir/build.make
Trabalho_04.exe: CMakeFiles/Trabalho_04.dir/linklibs.rsp
Trabalho_04.exe: CMakeFiles/Trabalho_04.dir/objects1.rsp
Trabalho_04.exe: CMakeFiles/Trabalho_04.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=C:\Users\fbeze\Desktop\Faculdade\Metodos_Numericos_II\TrabalhoMetodosII\Trabalho_04\cmake-build-debug\CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking C executable Trabalho_04.exe"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles\Trabalho_04.dir\link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/Trabalho_04.dir/build: Trabalho_04.exe

.PHONY : CMakeFiles/Trabalho_04.dir/build

CMakeFiles/Trabalho_04.dir/requires: CMakeFiles/Trabalho_04.dir/main.c.obj.requires

.PHONY : CMakeFiles/Trabalho_04.dir/requires

CMakeFiles/Trabalho_04.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles\Trabalho_04.dir\cmake_clean.cmake
.PHONY : CMakeFiles/Trabalho_04.dir/clean

CMakeFiles/Trabalho_04.dir/depend:
	$(CMAKE_COMMAND) -E cmake_depends "MinGW Makefiles" C:\Users\fbeze\Desktop\Faculdade\Metodos_Numericos_II\TrabalhoMetodosII\Trabalho_04 C:\Users\fbeze\Desktop\Faculdade\Metodos_Numericos_II\TrabalhoMetodosII\Trabalho_04 C:\Users\fbeze\Desktop\Faculdade\Metodos_Numericos_II\TrabalhoMetodosII\Trabalho_04\cmake-build-debug C:\Users\fbeze\Desktop\Faculdade\Metodos_Numericos_II\TrabalhoMetodosII\Trabalho_04\cmake-build-debug C:\Users\fbeze\Desktop\Faculdade\Metodos_Numericos_II\TrabalhoMetodosII\Trabalho_04\cmake-build-debug\CMakeFiles\Trabalho_04.dir\DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/Trabalho_04.dir/depend

