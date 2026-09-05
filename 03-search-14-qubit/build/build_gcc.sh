#!/bin/sh
# g++/clang++ release builds (the flags requested in the specification)
g++ -std=c++20 -O3 -march=native -flto -DNDEBUG -pthread stab14.cpp        -o stab14
g++ -std=c++20 -O3 -march=native -flto -DNDEBUG -pthread stab14_ext.cpp    -o stab14_ext
g++ -std=c++20 -O3 -march=native -flto -DNDEBUG -pthread stab14_2stage.cpp -o stab14_2stage
