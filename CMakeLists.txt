# CMakeList.txt: 顶层 CMake 项目文件，在此处执行全局配置
# 并包含子项目。
#
cmake_minimum_required(VERSION 3.11)

# set(CMAKE_TOOLCHAIN_FILE C:/vcpkg/scripts/buildsystems/vcpkg.cmake)
# 项目名称
project(Zaran VERSION 20221215.1 LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 17)

# 递归包含头文件的函数
function(include_sub_directories_recursively root_dir)
    if(IS_DIRECTORY ${root_dir}) # 当前路径是一个目录吗，是的话就加入到包含目录
        if(${root_dir} MATCHES "inc")
        #   include_directories(${root_dir})
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
add_executable(Zaran)
file(GLOB_RECURSE
    src_file
    Factory/*.cpp
    Basic/*.cpp
    Field/*.cpp
    Main/*.cpp)
message(STATUS "source files:" ${src_file})
target_sources(Zaran PRIVATE ${src_file})
target_include_directories(Zaran PRIVATE ${inc_dir})
# 程序包含路径

# 程序输出路径
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/bin/${CMAKE_BUILD_TYPE})

# OMP并行
option(USE_OMP "ON to use OMP, OFF to serial" ON)

if(USE_OMP)
    add_definitions(-DUSE_OMP)
    find_package(OpenMP REQUIRED)
    target_link_libraries(Zaran PRIVATE OpenMP::OpenMP_CXX)
endif()

# Eigen库
find_package(Eigen3 CONFIG REQUIRED)
target_link_libraries(Zaran PRIVATE Eigen3::Eigen)

# VTK库
find_package(VTK CONFIG REQUIRED)
target_link_libraries(Zaran PRIVATE ${VTK_LIBRARIES})
set(CMAKE_NINJA_FORCE_RESPONSE_FILE "ON" CACHE BOOL "Force Ninja to use response files.")

# Eigen库
find_package(spdlog CONFIG REQUIRED)
target_link_libraries(Zaran PRIVATE spdlog::spdlog spdlog::spdlog_header_only)
target_compile_options(Zaran PRIVATE $<$<CXX_COMPILER_ID:MSVC>:/bigobj> $<$<CXX_COMPILER_ID:GNU>:-Wa,-mbig-obj>)

# 可执行文件后缀
set_target_properties(Zaran PROPERTIES DEBUG_POSTFIX ${PROJECT_VERSION})
set_target_properties(Zaran PROPERTIES RELEASE_POSTFIX ${PROJECT_VERSION})