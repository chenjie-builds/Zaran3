# CMakeList.txt: 顶层 CMake 项目文件，在此处执行全局配置
# 并包含子项目。
#
cmake_minimum_required (VERSION 3.11)
#set(CMAKE_TOOLCHAIN_FILE C:/vcpkg/scripts/buildsystems/vcpkg.cmake)
#项目名称
project (Zaran VERSION 20221215.1 LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 17)

add_executable(Zaran )
aux_source_directory(src/Basic src_file)
aux_source_directory(src/Controller src_file)
aux_source_directory(src/FieldData src_file)
aux_source_directory(src/Geometry/Grid src_file)
aux_source_directory(src/Geometry/Model src_file)
aux_source_directory(src/Math src_file)
aux_source_directory(src/Numeric/RiemannSolver src_file)
aux_source_directory(src/PreProcess src_file)
aux_source_directory(src/ReferenceValue src_file)
aux_source_directory(src/Simulation src_file)
aux_source_directory(src/Solver src_file)
aux_source_directory(src/TimeCompute src_file)
aux_source_directory(src/Physics/Fluid src_file)
aux_source_directory(src/Utility src_file)
aux_source_directory(src/Visual src_file)
aux_source_directory(src/Zaran src_file)

aux_source_directory(inc/Basic inc_file)
aux_source_directory(inc/Controller inc_file)
aux_source_directory(inc/FieldData inc_file)
aux_source_directory(inc/Geometry/Grid inc_file)
aux_source_directory(inc/Geometry/Model inc_file)
aux_source_directory(inc/Math inc_file)
aux_source_directory(inc/Numeric/RiemannSolver inc_file)
aux_source_directory(inc/PreProcess inc_file)
aux_source_directory(inc/ReferenceValue inc_file)
aux_source_directory(inc/Simulation inc_file)
aux_source_directory(inc/Solver inc_file)
aux_source_directory(inc/TimeCompute inc_file)
aux_source_directory(inc/Physics/Fluid inc_file)
aux_source_directory(inc/Utility inc_file)
aux_source_directory(inc/Visual inc_file)
aux_source_directory(inc/Zaran inc_file)
message( STATUS "source files:" ${src_file} )

target_sources(Zaran PRIVATE  ${src_file}  ${inc_file})
target_include_directories(Zaran PRIVATE inc/Basic )
target_include_directories(Zaran PRIVATE inc/Controller )
target_include_directories(Zaran PRIVATE inc/FieldData )
target_include_directories(Zaran PRIVATE inc/Geometry/Grid )
target_include_directories(Zaran PRIVATE inc/Geometry/Model )
target_include_directories(Zaran PRIVATE inc/Math )
target_include_directories(Zaran PRIVATE inc/Numeric/RiemannSolver )
target_include_directories(Zaran PRIVATE inc/Physics/Fluid )
target_include_directories(Zaran PRIVATE inc/PreProcess )
target_include_directories(Zaran PRIVATE inc/ReferenceValue )
target_include_directories(Zaran PRIVATE inc/Simulation )
target_include_directories(Zaran PRIVATE inc/Solver )
target_include_directories(Zaran PRIVATE inc/TimeCompute )
target_include_directories(Zaran PRIVATE inc/Utility )
target_include_directories(Zaran PRIVATE inc/Visual )
message( STATUS "include files:" ${src_file} )
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
#Eigen库
find_package(Eigen3 CONFIG REQUIRED)
target_link_libraries(Zaran PRIVATE Eigen3::Eigen)
#VTK库
find_package(VTK CONFIG REQUIRED)
target_link_libraries(Zaran PRIVATE ${VTK_LIBRARIES})
set(CMAKE_NINJA_FORCE_RESPONSE_FILE "ON" CACHE BOOL "Force Ninja to use response files.")

#Eigen库
find_package(spdlog CONFIG REQUIRED)
target_link_libraries(Zaran PRIVATE spdlog::spdlog spdlog::spdlog_header_only)
target_compile_options(Zaran PRIVATE $<$<CXX_COMPILER_ID:MSVC>:/bigobj> $<$<CXX_COMPILER_ID:GNU>:-Wa,-mbig-obj>)

#可执行文件后缀
set_target_properties(Zaran PROPERTIES DEBUG_POSTFIX ${PROJECT_VERSION})
set_target_properties(Zaran PROPERTIES RELEASE_POSTFIX ${PROJECT_VERSION}) 