file( GLOB_RECURSE
    src_file
    Data/*.cpp
    DataManager/*.cpp
    Field/*.cpp
    Solver/*.cpp
    SolverParam/*.cpp
)
add_library( Zone ${src_file} )
target_include_directories(Zone PRIVATE ${inc_dir})

# Eigen库
find_package( Eigen3 CONFIG REQUIRED )

# VTK库
find_package( VTK REQUIRED )
find_package( spdlog CONFIG REQUIRED )

target_link_libraries(Zone PRIVATE Eigen3::Eigen)
target_link_libraries(Zone PRIVATE spdlog::spdlog)
target_link_libraries(Zone PRIVATE ${VTK_LIBRARIES})
target_link_libraries(Zone PRIVATE ${PROJECT_SOURCE_DIR}/3rdparty/tecio/lib/tecio.lib)
