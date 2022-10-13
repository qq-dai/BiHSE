cmake_minimum_required(VERSION 3.6)

project(run)

SET(CMAKE_CXX_COMPILER "g++")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -march=native -mavx -O3")

set(exName1 "MBC")

add_subdirectory(./src/BCE)
add_subdirectory(./src/biGraph)
add_subdirectory(./src/tools)


add_executable(${exName1} ./src/run/run.cpp)

target_link_libraries(${exName1} getArgs)
target_link_libraries(${exName1} biGraph)

target_link_libraries(${exName1} BCE)

set(binPath "bin")
ADD_CUSTOM_COMMAND(TARGET ${exName1}
          POST_BUILD
          COMMAND ${CMAKE_COMMAND} -E make_directory ${binPath}
          COMMAND ${CMAKE_COMMAND} -E copy ${exName1} ${binPath}/.
          COMMAND ${CMAKE_COMMAND} -E remove ${exName1}
)
