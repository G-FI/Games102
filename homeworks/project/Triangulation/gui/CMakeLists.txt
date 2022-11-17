cmake_minimum_required(VERSION 3.10)


include_directories("D:/include")
aux_source_directory(. gui_dir)

add_library(GUI STATIC)

target_sources(GUI PRIVATE ${gui_dir})