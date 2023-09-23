# CMakeList.txt: 顶层 CMake 项目文件，在此处执行全局配置
# 并包含子项目。
#
cmake_minimum_required (VERSION 3.11)
#set(CMAKE_TOOLCHAIN_FILE C:/vcpkg/scripts/buildsystems/vcpkg.cmake)
#项目名称
project (Zaran VERSION 20221215.1 LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 17)

FILE(GLOB_RECURSE src_file  src/*.cpp src/*.h)
message( STATUS "source files:" ${src_file} )

add_executable(Zaran ${src_file} )
target_include_directories(Zaran PUBLIC inc/Basic )
target_include_directories(Zaran PUBLIC inc/Controller )
target_include_directories(Zaran PUBLIC inc/FieldData )
target_include_directories(Zaran PUBLIC inc/Fluid )
target_include_directories(Zaran PUBLIC inc/Geometry/Grid )
target_include_directories(Zaran PUBLIC inc/Geometry/Model )
target_include_directories(Zaran PUBLIC inc/Numeric/RiemannSolver )
target_include_directories(Zaran PUBLIC inc/Numeric/Gradient )
target_include_directories(Zaran PUBLIC inc/Math )
target_include_directories(Zaran PUBLIC inc/Physics/Fluid )
target_include_directories(Zaran PUBLIC inc/ReferenceValue )
target_include_directories(Zaran PUBLIC inc/Simulation )
target_include_directories(Zaran PUBLIC inc/Solver )
target_include_directories(Zaran PUBLIC inc/TimeCompute )
target_include_directories(Zaran PUBLIC inc/Utility )
target_include_directories(Zaran PUBLIC inc/Visual )
target_include_directories(Zaran PUBLIC inc/PreProcess )
#程序包含路径

#程序输出路径
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/bin/${CMAKE_BUILD_TYPE})

#OMP并行
option(USE_OMP "ON to use OMP, OFF to serial" OFF)
if(USE_OMP)
add_definitions(-DUSE_OMP)
find_package(OpenMP REQUIRED)
target_link_libraries(Zaran PRIVATE OpenMP::OpenMP_CXX)
endif()
#VTK库
find_package(VTK CONFIG REQUIRED)
target_link_libraries(Zaran PRIVATE ${VTK_LIBRARIES})
set(CMAKE_NINJA_FORCE_RESPONSE_FILE "ON" CACHE BOOL "Force Ninja to use response files.")

#Eigen库
find_package(Eigen3 CONFIG REQUIRED)
target_link_libraries(Zaran PRIVATE Eigen3::Eigen)
#Eigen库
find_package(spdlog CONFIG REQUIRED)
target_link_libraries(Zaran PRIVATE spdlog::spdlog spdlog::spdlog_header_only)
target_compile_options(Zaran PRIVATE $<$<CXX_COMPILER_ID:MSVC>:/bigobj> $<$<CXX_COMPILER_ID:GNU>:-Wa,-mbig-obj>)

#可执行文件后缀
set_target_properties(Zaran PROPERTIES DEBUG_POSTFIX ${PROJECT_VERSION})
set_target_properties(Zaran PROPERTIES RELEASE_POSTFIX ${PROJECT_VERSION}) 