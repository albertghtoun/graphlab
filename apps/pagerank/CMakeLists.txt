cmake_minimum_required(VERSION 3.5)
project(GraphLab)

#see http://stackoverflow.com/questions/5096881/does-set-target-properties-in-cmake-override-cmake-cxx-flags
add_compile_options(-Wall) # for all targets in current directory
add_executable(pagerank pagerank.cpp)
#add_executable(tsv_to_graphlab_bin tsv_to_graphlab_bin.cpp)


#add_mex_library(pagerankmex pagerankapp.cpp  pagerankmex.cpp)
