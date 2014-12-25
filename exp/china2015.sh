#!/bin/bash  

rm ./fem_openmp 2> /dev/null
make -C ~/fem/ -f ~/fem/build/makefile.openmp clean
make -C ~/fem/ -f ~/fem/build/makefile.openmp build
cp ~/fem/fem_openmp ./fem_openmp
echo "Test for 24 threads"
export OMP_NUM_THREADS=24; ./fem_openmp --gtest_filter=cpu.openmp_experiment | tee -a openmp_result.log
echo "Test for 20 threads"
export OMP_NUM_THREADS=20; ./fem_openmp --gtest_filter=cpu.openmp_experiment | tee -a openmp_result.log
echo "Test for 16 threads"
export OMP_NUM_THREADS=16; ./fem_openmp --gtest_filter=cpu.openmp_experiment | tee -a openmp_result.log
echo "Test for 12 threads"
export OMP_NUM_THREADS=12; ./fem_openmp --gtest_filter=cpu.openmp_experiment | tee -a openmp_result.log
echo "Test for 8 threads"
export OMP_NUM_THREADS=8;  ./fem_openmp --gtest_filter=cpu.openmp_experiment | tee -a openmp_result.log
echo "Test for 6 threads"
export OMP_NUM_THREADS=6;  ./fem_openmp --gtest_filter=cpu.openmp_experiment | tee -a openmp_result.log
echo "Test for 4 threads"
export OMP_NUM_THREADS=4;  ./fem_openmp --gtest_filter=cpu.openmp_experiment | tee -a openmp_result.log
echo "Test for 2 threads"
export OMP_NUM_THREADS=2;  ./fem_openmp --gtest_filter=cpu.openmp_experiment | tee -a openmp_result.log
zip -r "log_"$(date -d "today" +"%Y%m%d%H%M") .
rm ./fem_openmp
rm openmp_result.log

