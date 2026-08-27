# CMakeList.txt: 顶层 CMake 项目文件，在此处执行全局配置
# 并包含子项目。
#
cmake_minimum_required( VERSION 3.21 )
#set( CMAKE_EXPORT_COMPILE_COMMANDS ON )

# set(CMAKE_TOOLCHAIN_FILE C:/vcpkg/scripts/buildsystems/vcpkg.cmake)
# 项目名称
project( Zaran LANGUAGES CXX )
set( PROJECT_VERSION_MAJOR 3 )
set( PROJECT_VERSION_MINOR 10 )
set( PROJECT_VERSION_PATCH 1 )
set( PROJECT_VERSION ${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR}.${PROJECT_VERSION_PATCH} )
set( EXECUTABLE_NAME ${PROJECT_NAME}${PROJECT_VERSION} )
set( CMAKE_CXX_STANDARD 17 )
set( inc_dir
    ${PROJECT_SOURCE_DIR}/inc/Zone/Data
    ${PROJECT_SOURCE_DIR}/inc/Zone/DataManager
    ${PROJECT_SOURCE_DIR}/inc/Zone/Field
    ${PROJECT_SOURCE_DIR}/inc/Zone/Solver
    ${PROJECT_SOURCE_DIR}/inc/Zone/SolverParam
    ${PROJECT_SOURCE_DIR}/inc/Generator/Field
    ${PROJECT_SOURCE_DIR}/inc/Generator/Grid
    ${PROJECT_SOURCE_DIR}/inc/Basic/Common
    ${PROJECT_SOURCE_DIR}/inc/Basic/Data
    ${PROJECT_SOURCE_DIR}/inc/Basic/Math
    ${PROJECT_SOURCE_DIR}/inc/Basic/Grid
    ${PROJECT_SOURCE_DIR}/inc/Basic/Model
    ${PROJECT_SOURCE_DIR}/inc/Basic/PostProcess
    ${PROJECT_SOURCE_DIR}/inc/Basic/PreProcess
    ${PROJECT_SOURCE_DIR}/inc/Basic/RiemannSolver
    ${PROJECT_SOURCE_DIR}/inc/Basic/Physics
    ${PROJECT_SOURCE_DIR}/inc/Basic/Utilities
    ${PROJECT_SOURCE_DIR}/inc/Main
)
add_executable( ${EXECUTABLE_NAME} )
file( GLOB_RECURSE
    src_file
    src/*.cpp
)
message( STATUS "source files:" ${src_file} )
message( STATUS "include directories:" ${inc_dir} )
target_sources( ${EXECUTABLE_NAME} PUBLIC ${src_file})
target_include_directories( ${EXECUTABLE_NAME} PRIVATE ${inc_dir} )
message( STATUS "CMAKE_BUILD_TYPE:" ${CMAKE_BUILD_TYPE} )
# 程序DEBUG模式
option (DEBUG_MODE "ON to enable debug info" OFF)
if (DEBUG_MODE OR CMAKE_BUILD_TYPE STREQUAL "Debug")
    target_compile_definitions( ${EXECUTABLE_NAME} PRIVATE DEBUG_MODE )
endif()
 
# 程序输出路径
set( EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/bin/${CMAKE_BUILD_TYPE} )

# OMP并行
option( USE_OMP "ON to use OMP, OFF to serial" ON )
message( STATUS "USE_OMP:" ${USE_OMP} )

if( USE_OMP )
    target_compile_definitions( ${EXECUTABLE_NAME} PRIVATE USE_OMP )
    find_package( OpenMP REQUIRED )
    target_link_libraries( ${EXECUTABLE_NAME} PRIVATE OpenMP::OpenMP_CXX )
    if(CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
        target_compile_options( ${EXECUTABLE_NAME} PRIVATE /openmp:llvm )
     elseif(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
        target_compile_options( ${EXECUTABLE_NAME} PRIVATE -fopenmp )
     elseif(CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
        target_compile_options( ${EXECUTABLE_NAME} PRIVATE -fopenmp )
     endif()
endif()

# Eigen库
set( Eigen_DIR  ${PROJECT_SOURCE_DIR}/3rdparty/Eigen/eigen3 )
target_include_directories( ${EXECUTABLE_NAME} PRIVATE ${Eigen_DIR} )

# VTK库
set( VTK_DIR  ${PROJECT_SOURCE_DIR}/3rdparty/vtk )
set( VTK_INC_DIR ${VTK_DIR}/include/vtk-9.3 )
set(VTK_BIN_DIR ${VTK_DIR}/bin)
set( VTK_LIB_DIR ${VTK_DIR}/lib)
file( GLOB VTK_LIBS "${VTK_LIB_DIR}/*.lib" )
file( GLOB VTK_DLLS "${VTK_BIN_DIR}/*.dll" )
target_include_directories( ${EXECUTABLE_NAME} PRIVATE ${VTK_INC_DIR} )
target_link_libraries( ${EXECUTABLE_NAME} PRIVATE ${VTK_LIBS} )
file(COPY ${VTK_DLLS} DESTINATION ${EXECUTABLE_OUTPUT_PATH})

# spdlog库
set( spdlog_DIR  ${PROJECT_SOURCE_DIR}/3rdparty/spdlog/include )
message(STATUS "spdlog_DIR: ${spdlog_DIR}")
target_include_directories( ${EXECUTABLE_NAME} PRIVATE ${spdlog_DIR} )

#cgns库
set( CGNS_DIR  ${PROJECT_SOURCE_DIR}/3rdparty/cgns )
set( CGNS_INC_DIR ${CGNS_DIR}/include )
set( CGNS_LIB_DIR ${CGNS_DIR}/lib )
file( GLOB CGNS_LIBS "${CGNS_LIB_DIR}/*.lib" )
file( GLOB CGNS_DLLS "${CGNS_LIB_DIR}/*.dll" )
target_include_directories( ${EXECUTABLE_NAME} PRIVATE ${CGNS_INC_DIR} )
target_link_libraries( ${EXECUTABLE_NAME} PRIVATE ${CGNS_LIBS} )
file(COPY ${CGNS_DLLS} DESTINATION ${EXECUTABLE_OUTPUT_PATH})

#hdf5库
set( HDF5_DIR  ${PROJECT_SOURCE_DIR}/3rdparty/hdf5 )
set( HDF5_INC_DIR ${HDF5_DIR}/include )
set( HDF5_LIB_DIR ${HDF5_DIR}/lib )
file( GLOB HDF5_LIBS "${HDF5_LIB_DIR}/*.lib" )
file( GLOB HDF5_DLLS "${HDF5_LIB_DIR}/*.dll" )
target_include_directories( ${EXECUTABLE_NAME} PRIVATE ${HDF5_INC_DIR} )
target_link_libraries( ${EXECUTABLE_NAME} PRIVATE ${HDF5_LIBS} )
file(COPY ${HDF5_DLLS} DESTINATION ${EXECUTABLE_OUTPUT_PATH})

#tecio库
set( Tecio_DIR  ${PROJECT_SOURCE_DIR}/3rdparty/tecio)
set( Tecio_INC_DIR ${Tecio_DIR}/include )
set( Tecio_LIB_DIR ${Tecio_DIR}/lib )
file( GLOB Tecio_LIBS "${Tecio_LIB_DIR}/*.lib" )
file( GLOB Tecio_DLLS "${Tecio_LIB_DIR}/*.dll" )
target_include_directories( ${EXECUTABLE_NAME} PRIVATE ${Tecio_INC_DIR} )
target_link_libraries( ${EXECUTABLE_NAME} PRIVATE ${Tecio_LIBS} )
file(COPY ${Tecio_DLLS} DESTINATION ${EXECUTABLE_OUTPUT_PATH})

#span-lite库
set( SpanLite_DIR  ${PROJECT_SOURCE_DIR}/3rdparty/span-lite )
set( SpanLite_INC_DIR ${SpanLite_DIR}/include/nonstd )
message(STATUS "SpanLite_INC_DIR: ${SpanLite_INC_DIR}")
target_include_directories( ${EXECUTABLE_NAME} PRIVATE ${SpanLite_INC_DIR})

#toml++库（单头文件，置于 3rdparty/tomlplusplus/toml.hpp）
set( TomlPlusPlus_DIR  ${PROJECT_SOURCE_DIR}/3rdparty/tomlplusplus )
set( TomlPlusPlus_INC_DIR ${TomlPlusPlus_DIR} )
message(STATUS "TomlPlusPlus_INC_DIR: ${TomlPlusPlus_INC_DIR}")
target_include_directories( ${EXECUTABLE_NAME} PRIVATE ${TomlPlusPlus_INC_DIR} )

set( CMAKE_NINJA_FORCE_RESPONSE_FILE "ON" CACHE BOOL "Force Ninja to use response files." )
message( STATUS "CXX_COMPILER_ID:" ${CMAKE_CXX_COMPILER_ID} )
if(CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
    target_compile_options( ${EXECUTABLE_NAME} PRIVATE /bigobj /MP /utf-8 )
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    target_compile_options( ${EXECUTABLE_NAME} PRIVATE -mbig-obj )
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
    target_compile_options( ${EXECUTABLE_NAME} PRIVATE -fopenmp )
endif()
message( STATUS "CMAKE_CXX_FLAGS:" ${CMAKE_CXX_FLAGS} )

#输出target_compile_definitions
get_target_property(compile_definitions ${EXECUTABLE_NAME} COMPILE_DEFINITIONS)
message(STATUS "compile_definitions: ${compile_definitions}")
