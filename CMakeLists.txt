# CMakeLists.txt for the "event" package. It creates a library with a dictionary and a main program.
# If ROOT is not installed in a default system location you need to tell CMake where to find it.
# Sourcing `thisroot.sh` already sets the required environment variables.
# Otherwise, you must tell the build system where to look for ROOT,
# for example by passing `-DROOT_DIR="/path/to/root/installation` at CMake configuration time.

  cmake_minimum_required(VERSION 3.19.6 FATAL_ERROR)
  project(thinslice)

# Locate the ROOT package and define a number of useful targets and variables.
  find_package(ROOT REQUIRED COMPONENTS RIO MathMore)
  if(NOT ROOT_FOUND)
         message(STATUS "ROOT package not found.")
         return()
  endif()

# Locate the RooUnfold package
  find_package(RooUnfold CONFIG REQUIRED)  

  find_package(PkgConfig REQUIRED)
  list(APPEND CMAKE_PREFIX_PATH $ENV{JSONCPP_FQ_DIR})
#  message(STATUS "CMAKE_PREFIX_PATH: ${CMAKE_PREFIX_PATH}")

  pkg_check_modules(JSONCPP jsoncpp)

# require C++11
  set (CMAKE_CXX_STANDARD 17)

  include_directories(${ROOT_INCLUDE_DIRS})
  include_directories(${RooUnfold_INCLUDE_DIRS})
  include_directories(${JSONCPP_INCLUDE_DIRS})
#   include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)
#   ROOT_GENERATE_DICTIONARY(G__Event Event.h LINKDEF EventLinkDef.h)

  link_directories(${RooUnfold_LIBRARY_DIRS})
  link_directories(${JSONCPP_LIBRARY_DIRS})
#  message(STATUS "JSONCPP_LIBRARY_DIRS: ${JSONCPP_LIBRARY_DIRS}")

  SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
  set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
  SET(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin)
  SET(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/lib)


# Create a shared library with a generated dictionary.
# Passing cmake targets such as `ROOT::RIO` as dependencies (rather than plain
# library names for example via ${ROOT_LIBRARIES}) ensures that properties such as required
# include directories and C++ standard are propagated to our libraries or executables.
# Note: To ensure compatibility with Cling, targets *must* be compiled using the
# same C++ standard as ROOT was compiled with.
  add_library(Event SHARED anavar.C BeamNtuple.C BeamVirtualDetector.C HadAna.cxx ThinSlice.cxx EffEval.cxx ProtonInel.cxx ProtonEnergy.cxx PionXsec.cxx BeamSampleAna.cxx util.cxx FillHist_Helper.cxx RecoDaughter.cxx Unfold.cxx TemplateFitter.cxx BetheBloch.cxx)
  target_link_libraries(Event PUBLIC ${ROOT_LIBRARIES} ${RooUnfold_LIBRARIES} ${JSONCPP_LIBRARIES})
#   target_compile_options(Event INTERFACE ${CMAKE_REQUIRED_FLAGS})
# Create the main program using the library.
  add_executable(RunCrossSection RunCrossSection.cxx)
  target_link_libraries(RunCrossSection Event)
  
  add_executable(BackgroundFit BackgroundFit.cxx)
  target_link_libraries(BackgroundFit Event)

  add_executable(RunProtonInel RunProtonInel.cxx)
  target_link_libraries(RunProtonInel Event)

  add_executable(RunProtonEnergy RunProtonEnergy.cxx)
  target_link_libraries(RunProtonEnergy Event)

  add_executable(RunPionXsec RunPionXsec.cxx)
  target_link_libraries(RunPionXsec Event)

  add_executable(RunBeamSampleAna RunBeamSampleAna.cxx)
  target_link_libraries(RunBeamSampleAna Event)

  add_executable(RunEffEval RunEffEval.cxx)
  target_link_libraries(RunEffEval Event)

  add_executable(RunCalcXS CalcXS.cxx)
  target_link_libraries(RunCalcXS Event)

  add_executable(RunBBStudies BBStudies.cxx)
  target_link_libraries(RunBBStudies Event)

  install(TARGETS RunCrossSection BackgroundFit RunProtonInel RunProtonEnergy RunPionXsec RunBeamSampleAna RunEffEval RunCalcXS RunBBStudies  DESTINATION bin)
  install(TARGETS Event DESTINATION lib)
