cmake_minimum_required(VERSION 2.8)

INCLUDE_DIRECTORIES(/usr/include/)

    SET(SRC_LIST ./src/main.cpp
        ./src/constant.h
        ./src/constant.cpp
        ./src/map.h
        ./src/map.cpp
        ./src/mole.h
        ./src/mole.cpp
        ./src/overlap_builder.h
        ./src/overlap_builder.cpp
        ./src/overlap.cpp
        ./src/runner.h
       )

ADD_EXECUTABLE(nanoARCS ${SRC_LIST})
TARGET_LINK_LIBRARIES(nanoARCS boost_filesystem boost_system boost_thread log4cxx)
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -O2 -g -Wno-unused-variable -pthread -lboost_system -lboost_filesystem -llog4cxx")
INSTALL(TARGETS nanoARCS
            RUNTIME DESTINATION bin
           )
