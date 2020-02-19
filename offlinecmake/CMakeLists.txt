# From https://root.cern.ch/faq/can-i-integrate-root-my-cmake-build

# CMakeLists.txt for event package. It creates a library with dictionary and a main program
cmake_minimum_required(VERSION 3.0 FATAL_ERROR)
project(myanalysis)

# You need to tell CMake where to find the ROOT installation. This can be done in a number of ways:
#   - ROOT built with classic configure/make use the provided $ROOTSYS/etc/cmake/FindROOT.cmake
#   - ROOT built with CMake. Add in CMAKE_PREFIX_PATH the installation prefix for ROOT
list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})

#---Locate the ROOT package and defines a number of variables (e.g. ROOT_INCLUDE_DIRS)
find_package(ROOT REQUIRED COMPONENTS Core RIO RooFit RooFitCore RooStats Hist Tree Matrix Physics MathCore Net)

#---Define useful ROOT functions and macros (e.g. ROOT_GENERATE_DICTIONARY)
include(${ROOT_USE_FILE})

include_directories(${CMAKE_SOURCE_DIR} ${ROOT_INCLUDE_DIRS} $ENV{ONLINE_MAIN}/include)
add_definitions(${ROOT_CXX_FLAGS})

#include source and include dir
include_directories(${PROJECT_SOURCE_DIR}/include)
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.h)


link_directories($ENV{ONLINE_MAIN}/lib)


#---Create a shared library with geneated dictionary
ROOT_GENERATE_DICTIONARY(G__libDataStruct ${PROJECT_SOURCE_DIR}/include/libDataStruct.h LINKDEF ${PROJECT_SOURCE_DIR}/include/libDataStructLinkDef.h)
add_library(libDataStruct SHARED ${PROJECT_SOURCE_DIR}/include/libDataStruct.h  G__libDataStruct.cxx)
target_link_libraries(libDataStruct ${ROOT_LIBRARIES})

#---Create a shared library with geneated dictionary
ROOT_GENERATE_DICTIONARY(G__libmyanalysis ${PROJECT_SOURCE_DIR}/include/myanalysis.h LINKDEF ${PROJECT_SOURCE_DIR}/include/myanalysisLinkDef.h)
add_library(myanalysis SHARED ${sources} ${headers}  G__libmyanalysis.cxx)
target_link_libraries(myanalysis pmonitor Event NoRootEvent message)
target_link_libraries(myanalysis libDataStruct)

#----------------------------------------------------------------------------
# Copy all resources including scripts, examples, text files, etc... to the build directory, i.e. the directory in which we
# build project. This is so that we can us the resources because it
# relies on these scripts being in the current working directory.
#
set(PROJ_DATA
  resource1.txt
  resource2.txt
  )


foreach(_script ${PROJ_DATA})
  configure_file(
    ${PROJECT_SOURCE_DIR}/${_script}
    ${PROJECT_BINARY_DIR}/${_script}
    COPYONLY
    )
endforeach()
