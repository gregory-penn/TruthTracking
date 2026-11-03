cmake_minimum_required(VERSION 3.22)

project(TruthTracking LANGUAGES CXX)

# Same C++ standard as k4ActsTracking and k4FWCore
set(CMAKE_CXX_STANDARD 20 CACHE STRING "")
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE RelWithDebInfo CACHE STRING "Build type" FORCE)
endif()

include(GNUInstallDirs)

# Built against the k4ActsTracking installed in the software image. KFRunner and the hit/seed helpers are header-only
# templates, the remaining functions they call are exported by libk4ActsTracking.so. The k4ActsTracking config also
# finds k4FWCore, Acts (Core, PluginDD4hep) and DD4hep.
find_package(Gaudi REQUIRED)
find_package(EDM4HEP REQUIRED)
find_package(k4FWCore REQUIRED)
find_package(Acts REQUIRED COMPONENTS Core)
find_package(TBB REQUIRED)
find_package(k4ActsTracking REQUIRED)

# Installed libraries find their dependencies through the paths they were linked with
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}")

# gaudi_add_module also generates the component list (TruthTracking.components) and the Python configurables
# (from Configurables import TruthTrackFinder); both are installed with the library.
gaudi_add_module(TruthTracking
  SOURCES
    src/TruthTrackFinder.cpp
  LINK
    Gaudi::GaudiKernel
    k4FWCore::k4FWCore
    EDM4HEP::edm4hep
    k4ActsTracking::k4ActsTracking
    Acts::Core
    TBB::tbb
)

install(TARGETS TruthTracking
  LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}"
)
