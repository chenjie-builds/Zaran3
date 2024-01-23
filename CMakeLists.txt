# CMakeList.txt: 顶层 CMake 项目文件，在此处执行全局配置
# 并包含子项目。
#
cmake_minimum_required(VERSION 3.11)

# set(CMAKE_TOOLCHAIN_FILE C:/vcpkg/scripts/buildsystems/vcpkg.cmake)
# 项目名称
project(Zaran LANGUAGES CXX)
set(PROJECT_VERSION_MAJOR 3)
set(PROJECT_VERSION_MINOR 0)
set(PROJECT_VERSION_PATCH 3)
set(PROJECT_VERSION ${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR}.${PROJECT_VERSION_PATCH})
set(EXECUTABLE_NAME ${PROJECT_NAME}${PROJECT_VERSION})
set(CMAKE_CXX_STANDARD 17)

# 递归包含头文件的函数
function(include_sub_directories_recursively root_dir)
    if(IS_DIRECTORY ${root_dir}) # 当前路径是一个目录吗，是的话就加入到包含目录
        if(${root_dir} MATCHES "inc")
            # include_directories(${root_dir})
            set(inc_dir1 ${root_dir} PARENT_SCOPE)
        endif()
    endif()

    file(GLOB ALL_SUB RELATIVE ${root_dir} ${root_dir}/*) # 获得当前目录下的所有文件，让如ALL_SUB列表中

    foreach(sub ${ALL_SUB})
        if(IS_DIRECTORY ${root_dir}/${sub})
            include_sub_directories_recursively(${root_dir}/${sub}) # 对子目录递归调用，包含
            set(inc_dir ${inc_dir} ${inc_dir1} PARENT_SCOPE)
        endif()
    endforeach()
endfunction()

include_sub_directories_recursively(${PROJECT_SOURCE_DIR})
add_executable(${EXECUTABLE_NAME})
file(GLOB_RECURSE
    src_file
    Factory/*.cpp
    Basic/*.cpp
    Field/*.cpp
    Main/*.cpp
    Test/*.cpp
    )
message(STATUS "source files:" ${src_file})
target_sources(${EXECUTABLE_NAME} PRIVATE ${src_file})
target_include_directories(${EXECUTABLE_NAME} PRIVATE ${inc_dir})

# 程序输出路径
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/bin/${CMAKE_BUILD_TYPE})

# OMP并行
option(USE_OMP "ON to use OMP, OFF to serial" ON)

if(USE_OMP)
add_definitions(-DUSE_OMP)
find_package(OpenMP REQUIRED)
target_link_libraries(${EXECUTABLE_NAME} PRIVATE OpenMP::OpenMP_CXX)
endif()

# Eigen库
find_package(Eigen3 CONFIG REQUIRED)
# VTK库
find_package(VTK CONFIG REQUIRED)
find_package(spdlog CONFIG REQUIRED)

target_link_libraries(${EXECUTABLE_NAME} PRIVATE ${PROJECT_SOURCE_DIR}/3rdparty/tecio/lib/tecio.lib)
target_link_libraries(${EXECUTABLE_NAME} PRIVATE Eigen3::Eigen ${VTK_LIBRARIES} spdlog::spdlog spdlog::spdlog_header_only)

set(CMAKE_NINJA_FORCE_RESPONSE_FILE "ON" CACHE BOOL "Force Ninja to use response files.")
target_compile_options(${EXECUTABLE_NAME} PRIVATE $<$<CXX_COMPILER_ID:MSVC>:/bigobj /MP /openmp:llvm> $<$<CXX_COMPILER_ID:GNU>:-Wa,-mbig-obj>)
