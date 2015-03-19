#!/bin/bash

rm ./fem_valgrind 2> /dev/null
rm ./valgrind_result.log 2> /dev/null
make -C ~/fem/ -f ~/fem/build/makefile.icpc clean
make -C ~/fem/ -f ~/fem/build/makefile.icpc build dbg=1
cp ~/fem/fem_icpc ./fem_valgrind

/home/jane/valgrind/bin/valgrind --tool=memcheck -v --leak-check=full ./fem_valgrind --gtest_filter=cpu.valgrind_test | tee -a valgrind_result.log
rm ./fem_valgrind
rm valgrind_result.log
