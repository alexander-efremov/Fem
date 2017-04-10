#!/bin/bash

rm ./fem_openmp 2> /dev/null ./openmp_result.log 2> /dev/null
make -C ~/fem/ -f ~/fem/build/makefile.openmp clean
make -C ~/fem/ -f ~/fem/build/makefile.openmp build
cp ~/fem/fem_openmp ./fem_openmp

export OMP_NUM_THREADS=16; ./fem_openmp --gtest_filter=cpu.omp_2_time_test_16_641 | tee -a openmp_result.log
zip -r $(date -d "today" +"%Y%m%d%H%M")"_onmp_log" openmp_result.log fem_openmp 
rm ./fem_openmp ./openmp_result.log