### Compilation using makefile
compile with:

    $ make

compile with debug flag

    $ make debug

clean with:

    $ make clean

### Using CMake

Create the build directory and enter there:

    $ mkdir build && cd build

then run cmake

    $ cmake ..

now you should be able to compile using `make`.

To enable the debug mode in CMake:

    $ cmake -DCMAKE_BUILD_TYPE=Debug ..
