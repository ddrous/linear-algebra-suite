cmake_minimum_required(VERSION 3.12)
project(matrices C)

set(CMAKE_C_STANDARD 99)

add_executable(matrices main.c)

target_link_libraries(matrices m)