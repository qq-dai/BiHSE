cmake_minimum_required(VERSION 3.6)

project(multiPivot)

SET(CMAKE_CXX_COMPILER "g++")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -march=native -mavx -O3")

set(exName1 "run")
# set(exName2 "def")

add_subdirectory(src/biGraph)
add_subdirectory(src/plex)
# add_subdirectory(src/defbiclique)
add_subdirectory(src/tools)

add_executable(${exName1} src/run/run.cpp)
# add_executable(${exName2} src/run/defbicliqueRunner.cpp)

target_link_libraries(${exName1} getArgs)
target_link_libraries(${exName1} biGraph)
target_link_libraries(${exName1} biplexv2)

# target_link_libraries(${exName2} getArgs)
# target_link_libraries(${exName2} biGraph)
# target_link_libraries(${exName2} defBiclique)

set(binPath "bin")
ADD_CUSTOM_COMMAND(TARGET ${exName1}  
# ${exName2} 
    POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy ${exName1} ${binPath}/.
    COMMAND ${CMAKE_COMMAND} -E remove ${exName1}
   #  COMMAND ${CMAKE_COMMAND} -E copy ${exName2} ${binPath}/.
   #  COMMAND ${CMAKE_COMMAND} -E remove ${exName2}
)
