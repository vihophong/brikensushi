# From https://root.cern.ch/faq/can-i-integrate-root-my-cmake-build

# CMakeLists.txt for event package. It creates a library with dictionary and a main program
cmake_minimum_required(VERSION 3.0 FATAL_ERROR)
project(main)

# You need to tell CMake where to find the ROOT installation. This can be done in a number of ways:
#   - ROOT built with classic configure/make use the provided $ROOTSYS/etc/cmake/FindROOT.cmake
#   - ROOT built with CMake. Add in CMAKE_PREFIX_PATH the installation prefix for ROOT
list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})

#---Locate the ROOT package and defines a number of variables (e.g. ROOT_INCLUDE_DIRS)
find_package(ROOT REQUIRED COMPONENTS Core RIO Hist Tree Matrix Physics MathCore Net)
#roofit enabled
#find_package(ROOT REQUIRED COMPONENTS Core RIO RooFit RooFitCore RooStats Hist Tree Matrix Physics MathCore Net)


#---Define useful ROOT functions and macros (e.g. ROOT_GENERATE_DICTIONARY)
include(${ROOT_USE_FILE})

include_directories(${CMAKE_SOURCE_DIR} ${ROOT_INCLUDE_DIRS})
add_definitions(${ROOT_CXX_FLAGS})


#Generate dictionaries
ROOT_GENERATE_DICTIONARY(G__DataStruct ${PROJECT_SOURCE_DIR}/include/DataStruct.hh LINKDEF ${PROJECT_SOURCE_DIR}/include/DataStructLinkDef.hh)
ROOT_GENERATE_DICTIONARY(G__BELEN ${PROJECT_SOURCE_DIR}/include/BELEN.hh LINKDEF ${PROJECT_SOURCE_DIR}/include/BELENLinkDef.hh)
ROOT_GENERATE_DICTIONARY(G__Clover ${PROJECT_SOURCE_DIR}/include/Clover.hh LINKDEF ${PROJECT_SOURCE_DIR}/include/CloverLinkDef.hh)
ROOT_GENERATE_DICTIONARY(G__Beam ${PROJECT_SOURCE_DIR}/include/Beam.hh LINKDEF ${PROJECT_SOURCE_DIR}/include/BeamLinkDef.hh)
ROOT_GENERATE_DICTIONARY(G__BRIKENTreeData ${PROJECT_SOURCE_DIR}/include/BRIKENTreeData.hh LINKDEF ${PROJECT_SOURCE_DIR}/include/BRIKENTreeDataLinkDef.hh)

add_library(DataStruct SHARED ${PROJECT_SOURCE_DIR}/src/DataStruct.cc G__DataStruct.cxx) # Link2Dictionary!
add_library(BELEN SHARED ${PROJECT_SOURCE_DIR}/src/BELEN.cc G__BELEN.cxx) # Link2Dictionary!
add_library(Clover SHARED  ${PROJECT_SOURCE_DIR}/src/Clover.cc G__Clover.cxx) # Link2Dictionary!
add_library(Beam SHARED  ${PROJECT_SOURCE_DIR}/src/Beam.cc G__Beam.cxx) # Link2Dictionary!
add_library(BRIKENTreeData SHARED  ${PROJECT_SOURCE_DIR}/src/BRIKENTreeData.cc G__BRIKENTreeData.cxx) # Link2Dictionary!
target_link_libraries(DataStruct ${ROOT_LIBRARIES})   # Link2Dictionary!
target_link_libraries(BELEN ${ROOT_LIBRARIES})   # Link2Dictionary!
target_link_libraries(Clover ${ROOT_LIBRARIES})   # Link2Dictionary!
target_link_libraries(Beam ${ROOT_LIBRARIES})   # Link2Dictionary!
target_link_libraries(BRIKENTreeData ${ROOT_LIBRARIES})   # Link2Dictionary!


#include source and include dir
include_directories(${PROJECT_SOURCE_DIR}/include)
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

#---Create  a main program using the library
add_executable(merger
    ${sources} ${headers}
    main.cc)
target_link_libraries(merger ${ROOT_LIBRARIES})
target_link_libraries(merger DataStruct)
target_link_libraries(merger BELEN)
target_link_libraries(merger Clover)
target_link_libraries(merger Beam)
target_link_libraries(merger BRIKENTreeData)


add_executable(belenreader
    ${sources} ${headers}
    belenreader.cc)
target_link_libraries(belenreader ${ROOT_LIBRARIES})
target_link_libraries(belenreader DataStruct)
target_link_libraries(belenreader BELEN)
target_link_libraries(belenreader Clover)
target_link_libraries(belenreader Beam)
target_link_libraries(belenreader BRIKENTreeData)

#----------------------------------------------------------------------------
# Copy all resources including scripts, examples, text files, etc... to the build directory, i.e. the directory in which we
# build project. This is so that we can us the resources because it
# relies on these scripts being in the current working directory.
#
set(PROJ_DATA
  resource1.txt
  He3_mapping.txt
  )


foreach(_script ${PROJ_DATA})
  configure_file(
    ${PROJECT_SOURCE_DIR}/${_script}
    ${PROJECT_BINARY_DIR}/${_script}
    COPYONLY
    )
endforeach()

