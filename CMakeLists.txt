cmake_minimum_required(VERSION 3.16)
project(evocirc C)

set(CMAKE_C_STANDARD 99)

set(CMAKE_C_FLAGS "-O3")

add_executable(evocirc main.c circuit.h heap.h types.h stb_perlin.h)