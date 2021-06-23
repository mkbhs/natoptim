# NatOptim
## Nature-Inspired Optimization Library for C++

### Install

Install CMake
```
$ sudo apt-get install cmake
```
Clone this git repository in your $HOME directory
```
$ git clone https://github.com/mkbhs/natoptim.git
```
Build NatOptim:
```
$ cd natoptim
$ mkdir build
$ cd build
$ cmake ..
$ make
```
### How to include NatOptim in your project:
In your project directory create a CMakeLists.txt: 
```
cmake_minimum_required(VERSION 3.5)

project(my_project)

find_library(NatOptim_LIB libNatOptim.so PATHS $ENV{HOME}/natoptim/build/)

include_directories($ENV{HOME}/natoptim/include/)

add_executable(my_project my_project.cpp)

target_link_libraries(my_project ${NatOptim_LIB})
``

