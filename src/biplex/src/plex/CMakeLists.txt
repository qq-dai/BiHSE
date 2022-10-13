include_directories(../biGraph ../tools)

# add_library(biplex biplex.h biplex.cpp)
# target_link_libraries(biplex biGraph linearSet)

add_library(biplexv2 biplexv2.h biplexv2.cpp)
target_link_libraries(biplexv2 biGraph linearSet)
