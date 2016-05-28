## Compile using CMake

Create the build directory and enter there:

    $ mkdir build && cd build

### configure with ICC

    $ cmake .. -DCMAKE_CXX_COMPILER=icc -DCMAKE_C_COMPILER=icc

### configure with system default compiler
  
    $ cmake .. 

### Enable debug mode in CMake:

    $ cmake -DCMAKE_BUILD_TYPE=Debug ..

now you should be able to compile using `make`.


