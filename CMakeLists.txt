cmake_minimum_required(VERSION 3.30)
project(laboratoria C)

set(CMAKE_C_STANDARD 11)

add_executable(lab1 lab1/lab1.c)
set_target_properties(lab1 PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/lab1 )

add_executable(lab2 lab2/lab2.c)
set_target_properties(lab2 PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/lab2 )

add_executable(lab3 lab3/lab3.c)
set_target_properties(lab3 PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/lab3 )

add_executable(lab4 lab4/lab4.c)
set_target_properties(lab4 PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/lab4 )

add_executable(lab5 lab5/lab5.c)
set_target_properties(lab5 PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/lab5 )