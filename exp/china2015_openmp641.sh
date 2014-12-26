#!/bin/bash  

rm ./fem_openmp 2> /dev/null
rm ./openmp_result.log 2> /dev/null
make -C ~/fem/ -f ~/fem/build/makefile.openmp clean
make -C ~/fem/ -f ~/fem/build/makefile.openmp build
cp ~/fem/fem_openmp ./fem_openmp

export OMP_NUM_THREADS=24; ./fem_openmp --gtest_filter=cpu.china2015_test | tee -a openmp_result.log
zip -r $(date -d "today" +"%Y%m%d%H%M")"_onmp_log" openmp_result.log fem_openmp 
rm ./fem_openmp
rm openmp_result.log
