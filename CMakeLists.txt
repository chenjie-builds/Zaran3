# CMakeList.txt: 顶层 CMake 项目文件，在此处执行全局配置
# 并包含子项目。
#
cmake_minimum_required( VERSION 3.21 )
set( CMAKE_EXPORT_COMPILE_COMMANDS ON )

# set(CMAKE_TOOLCHAIN_FILE C:/vcpkg/scripts/buildsystems/vcpkg.cmake)
# 项目名称
project( Zaran LANGUAGES CXX )
set( PROJECT_VERSION_MAJOR 3 )
set( PROJECT_VERSION_MINOR 8 )
set( PROJECT_VERSION_PATCH 0 )
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
set( LIBRARY_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/lib/${CMAKE_BUILD_TYPE} )

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
find_package( Eigen3 CONFIG REQUIRED )

# VTK库
find_package( VTK REQUIRED )
# spdlog库
find_package( spdlog CONFIG REQUIRED )
target_include_directories( ${EXECUTABLE_NAME} PRIVATE ${VTK_INCLUDE_DIRS} )
target_include_directories( ${EXECUTABLE_NAME} PRIVATE ${PROJECT_SOURCE_DIR}/3rdparty/tecio/inc )
target_link_libraries( ${EXECUTABLE_NAME} PRIVATE ${PROJECT_SOURCE_DIR}/3rdparty/tecio/lib/tecio.lib )
target_link_libraries( ${EXECUTABLE_NAME} PRIVATE Eigen3::Eigen ${VTK_LIBRARIES} spdlog::spdlog spdlog::spdlog_header_only )

set( CMAKE_NINJA_FORCE_RESPONSE_FILE "ON" CACHE BOOL "Force Ninja to use response files." )
message( STATUS "CXX_COMPILER_ID:" ${CMAKE_CXX_COMPILER_ID} )
if(CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
    target_compile_options( ${EXECUTABLE_NAME} PRIVATE /bigobj /MP /W0 )
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    target_compile_options( ${EXECUTABLE_NAME} PRIVATE -mbig-obj )
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
    target_compile_options( ${EXECUTABLE_NAME} PRIVATE -fopenmp )
endif()
message( STATUS "CMAKE_CXX_FLAGS:" ${CMAKE_CXX_FLAGS} )

#输出target_compile_definitions
get_target_property(compile_definitions ${EXECUTABLE_NAME} COMPILE_DEFINITIONS)
message(STATUS "compile_definitions: ${compile_definitions}")