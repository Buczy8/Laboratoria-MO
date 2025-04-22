cmake_minimum_required(VERSION 3.30)
project(laboratoria C CXX)

set(CMAKE_C_STANDARD 11)
set(CMAKE_CXX_STANDARD 23)

#add_executable(lab1 lab1/lab1.c)
#set_target_properties(lab1 PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/lab1 )
#
#add_executable(lab2 lab2/lab2.c)
#set_target_properties(lab2 PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/lab2 )
#
#add_executable(lab3 lab3/lab3.c)
#set_target_properties(lab3 PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/lab3 )
#
#add_executable(lab4 lab4/lab4.c)
#set_target_properties(lab4 PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/lab4 )
#
#add_executable(lab5 lab5/lab5.c)
#set_target_properties(lab5 PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/lab5 )
#
#add_executable(lab6 lab6/lab6.c)
#set_target_properties(lab6 PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/lab6 )
#
add_executable(lab7 lab7/lab7.c)
set_target_properties(lab7 PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/lab7 )

add_executable(lab8 lab8/lab8.cpp)
set_target_properties(lab8 PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/lab8 )