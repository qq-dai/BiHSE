include_directories(../biGraph ../tools)

add_library(BCE BCE.h BCE.cpp)
target_link_libraries(BCE biGraph linearSetThree)