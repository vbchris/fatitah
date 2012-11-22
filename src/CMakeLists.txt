# This is the root Vaango CMakeLists file.

project(Vaango)
cmake_minimum_required(VERSION 2.8 FATAL_ERROR)

# for the moment, force the user to select the build type here.  Maybe 
# we can enable all types again if dynamic compilation goes away.
# This needs to be before PROJECT, or changing the value will append to CMAKE_CONFIGURATION_TYPES, instead
# of replacing it. - BW
set(CMAKE_BUILD_TYPE Release CACHE string "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel.")
set(CMAKE_CONFIGURATION_TYPES ${BUILD_TYPE} CACHE string "Choose the type of build, options are: Debug Release RelWithDebInfo MinSizeRel." FORCE)

set(PACKAGE_SRC_PATH ${Vaango_SOURCE_DIR} CACHE IINTERNAL "Path to Vaango source tree.")
set(PACKAGE_LIB_PATH ${Vaango_BINARY_DIR}/lib CACHE INTERNAL "Path to Vaango Libraries.")

#-----------------------------------------------------------------------------
# Test for some required system information.
include (CMakeBackwardCompatibilityC)

#----------------------------------------------------------------------------
# shared/static
option(BUILD_SHARED_LIBS "Build shared libraries" ON)
mark_as_advanced(BUILD_SHARED_LIBS)

#----------------------------------------------------------------------------
# Output directories.
set (LIBRARY_OUTPUT_PATH ${Vaango_BINARY_DIR}/lib CACHE INTERNAL "Single output directory for building all libraries.")
set (EXECUTABLE_OUTPUT_PATH ${Vaango_BINARY_DIR} CACHE INTERNAL "Single output directory for building all executables.")

# compiler flags
if(CXX_FLAGS)
  set(CMAKE_CXX_FLAGS ${CXX_FLAGS} CACHE string "Default compiler flags" FORCE)
endif(CXX_FLAGS)

mark_as_advanced(LIBRARY_OUTPUT_PATH EXECUTABLE_OUTPUT_PATH)
#mark_as_advanced(LIBRARY_OUTPUT_PATH EXECUTABLE_OUTPUT_PATH CXX_TEST_PATH)

set(BUILD_FROM_TOP_LEVEL 0)

if(WIN32)
  if(NOT BUILD_SHARED_LIBS)
    add_definitions(-DBUILD_CORE_STATIC)
    add_definitions(-DGLEW_STATIC)
  endif(NOT BUILD_SHARED_LIBS)
  #  This defines M_PI for windows
  add_definitions(-D_USE_MATH_DEFINES)
  # defines Windows 2000 features like mouse leaving window
  add_definitions(-D_WIN32_WINNT=0x0500)
  # turns off VS 2005 warnings about strcpy et.al. being unsafe
  add_definitions(-D_CRT_SECURE_NO_DEPRECATE)
  add_definitions(-DWIN32_LEAN_AND_MEAN) # makes sure things like rad1/grp1 don't get defined
  add_definitions(-DNOMINMAX) # makes sure min and max don't get defined
  
  # win32 l3 warnings are a little ridiculous
  set(CMAKE_CXX_WARNING_LEVEL 1)
  if(CMAKE_CXX_FLAGS MATCHES "/W[0-4]")
    string(REGEX REPLACE "/W[0-4]" "/W1"
            CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
  else(CMAKE_CXX_FLAGS MATCHES "/W[0-4]")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /W1")
  endif(CMAKE_CXX_FLAGS MATCHES "/W[0-4]")
endif(WIN32)

#----------------------------------------------------------------------------------
# Configure options
#----------------------------------------------------------------------------------
option(MPM          "Build MPM, MPMArches and MPMICE" ON)
option(ICE          "Build ICE and MPMICE" ON)
option(ARCHES       "Build Arches" OFF)
option(WASATCH      "Build Wasatch" OFF)
option(RADIATION    "Build Arches/Radiation" OFF)
option(FORTRAN      "Build with fortran" OFF)

# Turn off Arches, Wasatch, ICE, MPM, or Radiation (for sci_defs/uintah_defs.h)
set(DEF_MPM "")
set(DEF_ICE "")
set(DEF_ARCHES "")
set(DEF_MODELS_RADIATION "")
set(DEF_WASATCH "")
set(DEF_NO_FORTRAN "")

if (MPM)
  set(BUILD_MPM "yes")
  set(USE_MPM   "yes")
else ()
  set(DEF_MPM "#define NO_MPM")
  set(BUILD_MPM "no")
  set(USE_MPM   "no")
endif()

if (ICE)
  set(BUILD_ICE "yes")
  set(USE_ICE   "yes")
else()
  set(DEF_ICE "#define NO_ICE")
  set(BUILD_ICE "no")
  set(USE_ICE   "no")
endif()

if (ARCHES)
  set(BUILD_ARCHES "yes")
  set(USE_ARCHES   "yes")
else()
  set(DEF_ARCHES "#define NO_ARCHES")
  set(BUILD_ARCHES "no")
  set(USE_ARCHES "no")
endif()

if (RADIATION)
  set(BUILD_MODELS_RADIATION "yes")
  set(USE_MODELS_RADIATION "yes")
else()
  set(DEF_MODELS_RADIATION "#define NO_MODELS_RADIATION")
  set(BUILD_MODELS_RADIATION "no")
  set(USE_MODELS_RADIATION "no")
endif()

if (WASATCH)
  set(BUILD_WASATCH "yes")
  set(USE_WASATCH   "yes")
else ()
  set(DEF_WASATCH "#define NO_WASATCH")
  set(BUILD_WASATCH "no")
  set(USE_WASATCH "no")
endif ()

# NO_FORTRAN is used in configVars.mk[.in] and DEF_NO_FORTRAN is used in uintah_testdefs.h.in
if (FORTRAN)
  set(NO_FORTRAN "no")
else ()
  set(NO_FORTRAN "yes")
  set(DEF_ARCHES "#define NO_ARCHES")
  set(DEF_MODELS_RADIATION "#define NO_MODELS_RADIATION")
  set(DEF_NO_FORTRAN "#define NO_FORTRAN")
endif ()

#--------------------------------------------------------------------------------
# Name mangling etc. for fortran files
# Most fortran compilers need an "_" for C function names.
# Some have strange underscore needs.  This flag is used in the
# tools/fspec.pl.in file to generate the correct
# cname of the fortran function.
#--------------------------------------------------------------------------------
enable_language(Fortran)

string(REGEX MATCH gfortran HAVE_GFORTRAN ${CMAKE_Fortran_COMPILER})
string(REGEX MATCH xlf HAVE_XLF ${CMAKE_Fortran_COMPILER})
string(REGEX MATCH pg77 HAVE_PG77 ${CMAKE_Fortran_COMPILER})
string(REGEX MATCH g77 HAVE_G77 ${CMAKE_Fortran_COMPILER})
string(REGEX MATCH ifort HAVE_ifORT ${CMAKE_Fortran_COMPILER})
string(REGEX MATCH f77 HAVE_F77 ${CMAKE_Fortran_COMPILER})

# DON'T USE find_library - many machine's paths are compiler specific.  The *right* thing would be
#   to compile a program looking for the fortran lib, but this is probably good enough
if(HAVE_GFORTRAN)
  set(F_LIBRARY gfortran CACHE string "fortran library")
  find_library(F_LIBRARY NAMES gfortran)
  set(FORTRAN_UNDERSCORE end CACHE string "What type of fortran underscore style - linux,end,none")
  set(DEF_FORTRAN_UNDERSCORE "#define FORTRAN_UNDERSCORE_END")
elseif(HAVE_XLF)
  set(F_LIBRARY xlf90 CACHE string "fortran library")
  find_library(F_LIBRARY NAMES xlf90)
  set(FORTRAN_UNDERSCORE none CACHE string "What type of fortran underscore style - linux,end,none")
  set(DEF_FORTRAN_UNDERSCORE "#define FORTRAN_UNDERSCORE_NONE")
elseif(HAVE_PGF77)
  set(F_LIBRARY pgftnrtl CACHE string "fortran library")
  find_library(F_LIBRARY NAMES pgftnrtl)
  set(FORTRAN_UNDERSCORE end CACHE string "What type of fortran underscore style - linux,end,none")
  set(DEF_FORTRAN_UNDERSCORE "#define FORTRAN_UNDERSCORE_END")
elseif(HAVE_G77)
  set(F_LIBRARY g2c CACHE string "fortran library")
  find_library(F_LIBRARY NAMES g2c)
  set(FORTRAN_UNDERSCORE linux CACHE string "What type of fortran underscore style - linux,end,none")
  set(DEF_FORTRAN_UNDERSCORE "#define FORTRAN_UNDERSCORE_LINUX")
elseif(HAVE_ifort)
  set(F_LIBRARY ifcore CACHE string "fortran library")
  find_library(F_LIBRARY NAMES ifcore)
  set(FORTRAN_UNDERSCORE end CACHE string "What type of fortran underscore style - linux,end,none")
  set(DEF_FORTRAN_UNDERSCORE "#define FORTRAN_UNDERSCORE_END")
elseif(HAVE_F77)
  # sgi needs ftn fortran
  #if(CMAKE_SYSTEM MATCHES "IRIX.*")
  #set(FORTRAN_UNDERSCORE end CACHE string "What type of fortran underscore style - linux,end,none")
  # dec needs for Ufor
  #if(CMAKE_SYSTEM MATCHES "OSF.*")
  # find_library(F_LIBRARY NAMES ........)
  #set(FORTRAN_UNDERSCORE end CACHE string "What type of fortran underscore style - linux,end,none")
  set(FORTRAN_UNDERSCORE end CACHE string "What type of fortran underscore style - linux,end,none")
  set(FORTRAN_LIBRARY "" CACHE string "fortran library")
  set(DEF_FORTRAN_UNDERSCORE "#define FORTRAN_UNDERSCORE_END")
endif()

# f77 on redstorm currently an exception - doesn't need it 
if(NOT F_LIBRARY AND NOT HAVE_F77)
  message(FATAL_ERROR "Cannot find fortran library")
endif(NOT F_LIBRARY AND NOT HAVE_F77)

find_program(PERL perl)
set(FSPEC ${CMAKE_CURRENT_BINARY_DIR}/genfspec.pl CACHE FILE "Path to fspec")
configure_file(${Vaango_SOURCE_DIR}/tools/fspec.pl.in ${FSPEC})

# Adds fortran file to list of srcs and creates the .h file for it
macro(FORTRAN_DEPENDENCY fortran_file lib)

  # should pass in a fortran_file and package it belongs to.  I.e.,
  # FORTRAN_DEPENDENCY(fortran/bcscalar_fort.h CCA_Components_Arches)
  string(REPLACE .F _fort.h header ${CMAKE_CURRENT_BINARY_DIR}/${fortran_file})
  string(REPLACE .F .fspec fspec_file ${fortran_file})
  file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/fortran)
  set(${lib}_FORTRAN_SRCS ${${lib}_FORTRAN_SRCS} ${fortran_file})
  # if everybody has cmake 2.4.5, we can do IS_NEWER_THAN, and not create these files on every cmake
  if(${CMAKE_CURRENT_SOURCE_DIR}/${fspec_file} IS_NEWER_THAN ${header})
    if(NOT EXISTS ${header})
      execute_process(COMMAND perl ${FSPEC} ${CMAKE_CURRENT_SOURCE_DIR}/${fspec_file} ${header})
    endif(NOT EXISTS ${header})
  endif(${CMAKE_CURRENT_SOURCE_DIR}/${fspec_file} IS_NEWER_THAN ${header})

endmacro(FORTRAN_DEPENDENCY)

# Find the MPI package to use
find_package(MPI REQUIRED)
if (MPI_FOUND)
  include_directories(${MPI_INCLUDE_PATH})
  set(MPI_DIR "" CACHE PATH "Path to MPI")
  set(MPI_LIBRARY ${MPI_LIBRARIES})
else()
  message(STATUS "No mpi found")
  set(MPI_DIR "")
  set(MPI_LIBRARY "")
endif()

# Find the LibXML2 package to use
find_package(LibXml2 REQUIRED)
if (LIBXML2_FOUND)
  include_directories(${LIBXML2_INCLUDE_DIR})
  set(LIBXML2_DIR "" CACHE PATH "Path to MPI")
  set(XML2_INCLUDE ${LIBXML2_INCLUDE_DIR})
  set(XML2_LIBRARY ${LIBXML2_LIBRARIES})
else()
  message(STATUS "No libxml2 found")
  set(XML2_INCLUDE "")
  set(XML2_LIBRARY "")
endif()

#-----------------------------------------------------------------------------
# Configure time cross-platform header file generation
#-----------------------------------------------------------------------------

# 32-64 bit ?
if (CMAKE_SIZEOF_VOID_P EQUAL 8)
  set(NBITS 64)
else ()
  set(NBITS 32)
endif()

# SCI_MALLOC on?
if (SCI_MALLOC_ON)
  set(DEF_MALLOC "")
  set(DEF_MALLOC_LINENUMBERS "")
else()
  set(DEF_MALLOC "#define DISABLE_SCI_MALLOC 1")
  set(DEF_MALLOC_LINENUMBERS "")
endif()

# SCI Assertions
if (CMAKE_BUILD_TYPE STREQUAL "Release")
  set(THROW_LEVEL  0)
  set(ASSERTION_LEVEL  0)
else()
  set(THROW_LEVEL  0)
  set(ASSERTION_LEVEL  3)
endif()

# SCI Template definitions
set(TEMPLATE_TAG "")
set(TEMPLATE_BOX "<>")

# Threads library  (for populating sci_defs/threads_def.h)
set(CMAKE_THREAD_PREFER_PTHREAD true)
find_package(Threads REQUIRED)
if (${Threads_FOUND})
  set(DEF_THREAD "#define SCI_PTHREAD 1")
  set(THREAD_LIB_FLAG "${CMAKE_THREAD_LIBS_INIT}")
  set(THREAD_LIBRARY "${CMAKE_THREAD_LIBS_INIT}")
  if (CMAKE_USE_PTHREADS_INIT)
    set(THREADS_TYPE  "pthreads")
  elseif (CMAKE_USE_WIN32_THREADS_INIT)
    set(THREADS_TYPE  "win32")
  else()
    set(THREADS_TYPE  "unix")
  endif()
  message(STATUS "Threads library ${THREADS_TYPE} found with flag ${CMAKE_THREAD_LIBS_INIT}")
else()
  set(DEF_THREAD "#define SCI_NOTHREAD 1")
  set(THREAD_LIB_FLAG "")
  set(THREAD_LIBRARY "")
  message(STATUS "No threads library found.  Software may work incorrectly.")
endif()
set(THREAD_IMPL "Thread_${THREADS_TYPE}.cc")

# Check availability of semaphore.h
check_include_file(semaphore.h HAVE_SEMAPHORE_H)

# Check atomic builtins for threads support in gcc
include(CheckCXXSourceCompiles)

CHECK_CXX_SOURCE_COMPILES("
#ifdef _MSC_VER
#include <windows.h>
#endif
int main() {
#ifdef _MSC_VER
        volatile LONG val = 1;
        MemoryBarrier();
        InterlockedCompareExchange(&val, 0, 1);
        InterlockedIncrement(&val);
        InterlockedDecrement(&val);
#else
        volatile unsigned long val = 1;
        __sync_synchronize();
        __sync_val_compare_and_swap(&val, 1, 0);
        __sync_add_and_fetch(&val, 1);
        __sync_sub_and_fetch(&val, 1);
#endif
        return 0;
      }
" GCC_HAS_ATOMICS)

if( NOT GCC_HAS_ATOMICS )
  message(STATUS "Warning: GCC will be built thread-unsafe because atomic builtins are missing")
  set(REFCOUNT_IMPL "RefCounted_gcc.cc")
  set(ATOMIC_IMPL "AtomicCounter_gcc.cc")
else()
  set(REFCOUNT_IMPL "RefCounted.cc")
  set(ATOMIC_IMPL "AtomicCounter_default.cc")
endif()

# time implementation 
if (WIN32)
  set(TIME_IMPL "Time_win32.cc")
else()
  set(TIME_IMPL "Time_unix.cc")
endif()

# Find dl library (for dlopen())
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")
find_package(DL)
if (${DL_FOUND})
  set(DL_LIBRARY "${DL_LIBRARIES}")
  message(STATUS "DL library found")
else()
  set(DL_LIBRARY "")
  message(STATUS "DL library not found")
endif()

# Check included header files and create configuration files (for populating sci_defs/config_defs.h)
include(CheckIncludeFiles)
#check_include_files("ext/algorithm" HAVE_EXT_ALGORITHM)
CHECK_CXX_SOURCE_COMPILES("
#include <ext/algorithm>
#include <vector>
int main() {
        std::vector<int> vec;
        return 0;
      }
" HAVE_EXT_ALGORITHM)
if (HAVE_EXT_ALGORITHM)
  message(STATUS "ext/algorithm found.")
endif()
check_include_files(inttypes.h HAVE_INTTYPES_H)
check_include_files(limits.h HAVE_LIMITS_H)
check_include_files(memory.h HAVE_MEMORY_H)
check_include_files(stdint.h HAVE_STDINT_H)
check_include_files(stdlib.h HAVE_STDLIB_H)
check_include_files(strings.h HAVE_STRINGS_H)
check_include_files(string.h HAVE_STRING_H)
check_include_files("sys/select.h" HAVE_SYS_SELECT_H)
check_include_files("sys/stat.h" HAVE_SYS_STAT_H)
check_include_files("sys/time.h" HAVE_SYS_TIME_H)
check_include_files("sys/types.h" HAVE_SYS_TYPES_H)
check_include_files(unistd.h HAVE_UNISTD_H)
if (${HAVE_STDLIB_H})
  set(STDC_HEADERS 1)
endif()
message(STATUS "STDC_HEADERS ${STDC_HEADERS}")

# Check type sizes for config_defs.h
check_type_size("void*" SIZEOF_VOIDP)
check_type_size("long long" SIZEOF_LONG_LONG)
message(STATUS "SIZEOF_VOIDP ${SIZEOF_VOIDP}")
message(STATUS "SIZEOF_LONG_LONG ${SIZEOF_LONG_LONG}")

# Package info for config_defs.h
set(PACKAGE_BUGREPORT "")
set(PACKAGE_NAME "Vaango")
set(PACKAGE_STRING "Vaango 1.0")
set(PACKAGE_TARNAME "vaango")
set(PACKAGE_URL "")
set(PACKAGE_VERSION "1.0")

# Check for location of hash_map for hashmap_defs.h
message(STATUS "Looking for hash_map")
set(HAVE_HASHMAP "None")
if (NOT OSX)
  CHECK_CXX_SOURCE_COMPILES("
  #include <tr1/unordered_map>
  int main() { std::tr1::unordered_map<int, int> xx; return 0; }
  " DEF_HASHMAP)
  if(DEF_HASHMAP)
    message(STATUS "tr1/unordered_map found")
    set(HAVE_HASHMAP "HAVE_TR1_HASHMAP")
  endif()
endif()
## WARNING - the hashmap checks below don't work yet. BB
if (HAVE_HASHMAP STREQUAL "None")
  CHECK_CXX_SOURCE_COMPILES("
  #include <ext/hash_map>
  int main() { __gnu_cxx::hash_map<int, int> xx; return 0; }
  " DEF_HASHMAP)
  if(DEF_HASHMAP)
    message(STATUS "__gnu_cxx ext/hash_map found")
    set(HAVE_HASHMAP "HAVE_GNU_HASHMAP")
  endif()
endif()
if (HAVE_HASHMAP STREQUAL "None")
  message(STATUS "__gnu_cxx hash_map not found. Trying hash_map.")
CHECK_CXX_SOURCE_COMPILES("
#include <hash_map>
int main() { std::hash_map<int, int> xx; return 0; }
" DEF_HASHMAP)
  if(DEF_HASHMAP)
    message(STATUS "Standard hash_map found")
    set(HAVE_HASHMAP "HAVE_STD_HASHMAP")
  endif()
endif()
if (HAVE_HASHMAP STREQUAL "None")
  message(STATUS "Standard hash_map not found. Trying ext/hash_map.")
  CHECK_CXX_SOURCE_COMPILES("
  #include <ext/hash_map>
  int main() { std::hash_map<int, int> xx; return 0; }
  " DEF_HASHMAP)
  if(DEF_HASHMAP)
    message(STATUS "ext/hash_map found")
    set(HAVE_HASHMAP "HAVE_EXT_HASHMAP")
  else()
    message(STATUS "ext/hash_map not found. Trying stdext/hash_map.")
    CHECK_CXX_SOURCE_COMPILES("
    #include <hash_map>
    int main() { stdext::hash_map<int, int> xx; return 0; }
    " DEF_HASHMAP)
    if(DEF_HASHMAP)
      message(STATUS "stdext/hash_map found")
      set(HAVE_HASHMAP "HAVE_STDEXT_HASHMAP")
    else()
      message(STATUS "stdext/hash_map not found. Trying tr1/unordered_map.")
    endif()
  endif()
endif()
if (DEF_HASHMAP)
  set(DEF_HASHMAP "#define ${HAVE_HASHMAP} 1")
else()
  message(STATUS "hash_map not found, will use map instead.")
  set(DEF_HASHMAP "")
endif()

# Check for required standard headers
check_include_files(unistd.h HAVE_UNISTD_H)
check_include_files(dirent.h HAVE_DIRENT_H)
check_include_files(sys/time.h HAVE_SYS_TIME_H)
check_include_files(sys/param.h HAVE_SYS_PARAM_H)
check_include_files(values.h HAVE_VALUES_H)
check_include_files(malloc.h HAVE_MALLOC_H)
check_include_files(netdb.h HAVE_NETDB_H)
check_include_files(sys/select.h HAVE_SYS_SELECT_H)
check_include_files(sys/socket.h HAVE_SYS_SOCKET_H)
check_include_files(sys/mman.h HAVE_SYS_MMAN_H)
check_include_files(sys/ioctl.h HAVE_SYS_IOCTL_H)
check_include_files(sys/resource.h HAVE_SYS_RESOURCE_H)
check_include_files(sys/wait.h HAVE_SYS_WAIT_H)
check_include_files(sys/utsname.h HAVE_SYS_UTSNAME_H)
check_include_files(rpc/types.h HAVE_RPC_TYPES_H)
check_include_files(netinet/in.h HAVE_NETINET_IN_H)

# Check for standard template library headers
CHECK_CXX_SOURCE_COMPILES("
#include <algorithm>
#include <vector>
#include <map>
#include <iostream>
#include <string>
int main() {
        std::vector<int> vec;
        return 0;
      }
" HAVE_STL_HEADERS)
if (HAVE_STL_HEADERS)
  message(STATUS "Standard template library headers found.")
endif()

# Create the cross-platform header files
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)
include_directories(${CMAKE_CURRENT_BINARY_DIR}/include)
include_directories(${CMAKE_CURRENT_BINARY_DIR}/include/sci_defs)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/include/sci_defs/config_testdefs.h.in ${CMAKE_CURRENT_BINARY_DIR}/include/sci_defs/config_defs.h)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/include/sci_defs/audio_testdefs.h.in ${CMAKE_CURRENT_BINARY_DIR}/include/sci_defs/audio_defs.h)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/include/sci_defs/babel_testdefs.h.in ${CMAKE_CURRENT_BINARY_DIR}/include/sci_defs/babel_defs.h)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/include/sci_defs/bits_testdefs.h.in ${CMAKE_CURRENT_BINARY_DIR}/include/sci_defs/bits_defs.h)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/include/sci_defs/blas_testdefs.h.in ${CMAKE_CURRENT_BINARY_DIR}/include/sci_defs/blas_defs.h)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/include/sci_defs/boost_testdefs.h.in ${CMAKE_CURRENT_BINARY_DIR}/include/sci_defs/boost_defs.h)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/include/sci_defs/chromium_testdefs.h.in ${CMAKE_CURRENT_BINARY_DIR}/include/sci_defs/chromium_defs.h)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/include/sci_defs/collab_vis_testdefs.h.in ${CMAKE_CURRENT_BINARY_DIR}/include/sci_defs/collab_vis_defs.h)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/include/sci_defs/compile_testdefs.h.in ${CMAKE_CURRENT_BINARY_DIR}/include/sci_defs/compile_defs.h)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/include/sci_defs/crypto_testdefs.h.in ${CMAKE_CURRENT_BINARY_DIR}/include/sci_defs/crypto_defs.h)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/include/sci_defs/cuda_testdefs.h.in ${CMAKE_CURRENT_BINARY_DIR}/include/sci_defs/cuda_defs.h)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/include/sci_defs/dataflow_testdefs.h.in ${CMAKE_CURRENT_BINARY_DIR}/include/sci_defs/dataflow_defs.h)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/include/sci_defs/dynamic_cast_testdefs.h.in ${CMAKE_CURRENT_BINARY_DIR}/include/sci_defs/dynamic_cast_defs.h)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/include/sci_defs/environment_testdefs.h.in ${CMAKE_CURRENT_BINARY_DIR}/include/sci_defs/environment_defs.h)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/include/sci_defs/error_testdefs.h.in ${CMAKE_CURRENT_BINARY_DIR}/include/sci_defs/error_defs.h)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/include/sci_defs/exe_testdefs.h.in ${CMAKE_CURRENT_BINARY_DIR}/include/sci_defs/exe_defs.h)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/include/sci_defs/framework_testdefs.h.in ${CMAKE_CURRENT_BINARY_DIR}/include/sci_defs/framework_defs.h)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/include/sci_defs/globus_testdefs.h.in ${CMAKE_CURRENT_BINARY_DIR}/include/sci_defs/globus_defs.h)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/include/sci_defs/gperftools_testdefs.h.in ${CMAKE_CURRENT_BINARY_DIR}/include/sci_defs/gperftools_defs.h)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/include/sci_defs/hashmap_testdefs.h.in ${CMAKE_CURRENT_BINARY_DIR}/include/sci_defs/hashmap_defs.h)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/include/sci_defs/hdf5_testdefs.h.in ${CMAKE_CURRENT_BINARY_DIR}/include/sci_defs/hdf5_defs.h)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/include/sci_defs/hypre_testdefs.h.in ${CMAKE_CURRENT_BINARY_DIR}/include/sci_defs/hypre_defs.h)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/include/sci_defs/ieeefp_testdefs.h.in ${CMAKE_CURRENT_BINARY_DIR}/include/sci_defs/ieeefp_defs.h)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/include/sci_defs/image_testdefs.h.in ${CMAKE_CURRENT_BINARY_DIR}/include/sci_defs/image_defs.h)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/include/sci_defs/kepler_testdefs.h.in ${CMAKE_CURRENT_BINARY_DIR}/include/sci_defs/kepler_defs.h)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/include/sci_defs/lapack_testdefs.h.in ${CMAKE_CURRENT_BINARY_DIR}/include/sci_defs/lapack_defs.h)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/include/sci_defs/loki_testdefs.h.in ${CMAKE_CURRENT_BINARY_DIR}/include/sci_defs/loki_defs.h)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/include/sci_defs/malloc_testdefs.h.in ${CMAKE_CURRENT_BINARY_DIR}/include/sci_defs/malloc_defs.h)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/include/sci_defs/mdsplus_testdefs.h.in ${CMAKE_CURRENT_BINARY_DIR}/include/sci_defs/mdsplus_defs.h)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/include/sci_defs/mpi_testdefs.h.in ${CMAKE_CURRENT_BINARY_DIR}/include/sci_defs/mpi_defs.h)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/include/sci_defs/osx_testdefs.h.in ${CMAKE_CURRENT_BINARY_DIR}/include/sci_defs/osx_defs.h)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/include/sci_defs/papi_testdefs.h.in ${CMAKE_CURRENT_BINARY_DIR}/include/sci_defs/papi_defs.h)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/include/sci_defs/petsc_testdefs.h.in ${CMAKE_CURRENT_BINARY_DIR}/include/sci_defs/petsc_defs.h)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/include/sci_defs/ptolemy_testdefs.h.in ${CMAKE_CURRENT_BINARY_DIR}/include/sci_defs/ptolemy_defs.h)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/include/sci_defs/ruby_testdefs.h.in ${CMAKE_CURRENT_BINARY_DIR}/include/sci_defs/ruby_defs.h)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/include/sci_defs/scisock_testdefs.h.in ${CMAKE_CURRENT_BINARY_DIR}/include/sci_defs/scisock_defs.h)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/include/sci_defs/ssl_testdefs.h.in ${CMAKE_CURRENT_BINARY_DIR}/include/sci_defs/ssl_defs.h)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/include/sci_defs/stat64_testdefs.h.in ${CMAKE_CURRENT_BINARY_DIR}/include/sci_defs/stat64_defs.h)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/include/sci_defs/tao_testdefs.h.in ${CMAKE_CURRENT_BINARY_DIR}/include/sci_defs/tao_defs.h)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/include/sci_defs/teem_testdefs.h.in ${CMAKE_CURRENT_BINARY_DIR}/include/sci_defs/teem_defs.h)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/include/sci_defs/template_testdefs.h.in ${CMAKE_CURRENT_BINARY_DIR}/include/sci_defs/template_defs.h)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/include/sci_defs/tena_testdefs.h.in ${CMAKE_CURRENT_BINARY_DIR}/include/sci_defs/tena_defs.h)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/include/sci_defs/tena_testdefs.in ${CMAKE_CURRENT_BINARY_DIR}/include/sci_defs/tena_defs.h)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/include/sci_defs/thread_testdefs.h.in ${CMAKE_CURRENT_BINARY_DIR}/include/sci_defs/thread_defs.h)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/include/sci_defs/uintah_testdefs.h.in ${CMAKE_CURRENT_BINARY_DIR}/include/sci_defs/uintah_defs.h)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/include/sci_defs/vdt_testdefs.h.in ${CMAKE_CURRENT_BINARY_DIR}/include/sci_defs/vdt_defs.h)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/include/sci_defs/vtk_testdefs.h.in ${CMAKE_CURRENT_BINARY_DIR}/include/sci_defs/vtk_defs.h)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/include/sci_defs/wx_testdefs.h.in ${CMAKE_CURRENT_BINARY_DIR}/include/sci_defs/wx_defs.h)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/include/sci_defs/z_testdefs.h.in ${CMAKE_CURRENT_BINARY_DIR}/include/sci_defs/z_defs.h)

#-----------------------------------------------------------------------------
# Packages to enable (must enable packages when building core for now)
set(LOAD_PACKAGE "Vaango" CACHE string "Building Vaango package by default.")
set(HANDLED_PACKAGE_Vaango "YES")

include_directories(
  ${PETSC_INCLUDES}
  ${HYPRE_INCLUDE}
  ${XML2_INCLUDE}
)

include_directories(
  ${Vaango_SOURCE_DIR}
  ${CMAKE_CURRENT_BINARY_DIR}
)

add_subdirectory(Core)
add_subdirectory(CCA)
add_subdirectory(StandAlone)

