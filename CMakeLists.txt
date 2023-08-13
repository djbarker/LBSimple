cmake_minimum_required(VERSION 3.0)

project(lbsimple)

list(APPEND CMAKE_MODULE_PATH "/home/dan/workspace/LBSimple")

set(CMAKE_BUILD_TYPE RELEASE)
# set(CMAKE_BUILD_TYPE DEBUG)
# set(CMAKE_BUILD_TYPE RELWITHDEBINFO)

file(GLOB_RECURSE LBM_SOURCES  "src/*.cpp")
file(GLOB_RECURSE LBM_INCLUDES "include/*.h*")


set(CXX_OPT_FLAGS "-O3")
set(CXX_OPT_FLAGS "${CXX_OPT_FLAGS} -march=native") 

set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} ${CXX_OPT_FLAGS} ${CXX_OPENMP_FLAGS}")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_RELEASE} ${CXX_OPT_FLAGS} ${CXX_OPENMP_FLAGS}")
# set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -pg")
# set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -fopt-info")

add_compile_options(--std=c++1y)
add_definitions(-D_GLIBCXX_USE_CXX11_ABI=0)

include_directories(
	${CMAKE_SOURCE_DIR}
	${CMAKE_SOURCE_DIR}/include
	)



find_package(VTK 6.1 REQUIRED NO_MODULE)
if(VTK_FOUND)
	message("found VTK. Version:" ${VTK_VERSION}. VTK_DIR: ${VTK_DIR})
endif()
include(${VTK_USE_FILE})

# find_package( Boost 1.54 COMPONENTS program_options REQUIRED )
# include_directories(${Boost_INCLUDE_DIR})

# add_executable(lbsimple ${LBM_SOURCES})
# target_link_libraries(lbsimple ${VTK_LIBRARIES} ${pybind11_LIBRARIES}) # ${Boost_LIBRARIES})
# target_link_libraries(lbsimple fmt::fmt-header-only)

add_subdirectory("deps/fmt")
include_directories(${fmt_INCLUDE_DIRS})
include_directories("deps/fmt/include") 

find_package(pybind11 CONFIG)
include_directories(${pybind11_INCLUDE_DIRS})
pybind11_add_module(lbsimple ${LBM_SOURCES})

target_link_libraries(lbsimple PRIVATE ${VTK_LIBRARIES})

# get_cmake_property(_variableNames VARIABLES)
# foreach (_variableName ${_variableNames})
# 	message(STATUS "${_variableName}=${${_variableName}}")
# endforeach()

