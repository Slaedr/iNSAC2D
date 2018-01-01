
add_library(insac2ds_base gradientschemesins.cpp inviscidflux.cpp structmesh2d.cpp aarray2d.cpp)

add_executable(insac2ds insacs.cpp)
target_link_libraries(insac2ds insac2ds_base)
