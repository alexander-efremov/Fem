#!/bin/bash  

rm ./fem_openmp 2> /dev/null
rm ./openmp_result.log 2> /dev/null
make -C ~/fem/ -f ~/fem/build/makefile.openmp clean
make -C ~/fem/ -f ~/fem/build/makefile.openmp build
cp ~/fem/fem_openmp ./fem_openmp

rm ./fem_icpc 2> /dev/null
rm ./icpc_result.log 2> /dev/null
make -C ~/fem/ -f ~/fem/build/makefile.icpc clean
make -C ~/fem/ -f ~/fem/build/makefile.icpc build
cp ~/fem/fem_icpc ./fem_icpc

export OMP_NUM_THREADS=20; ./fem_openmp --gtest_filter=cpu.china2015_openmp_test | tee -a openmp_result.log
export OMP_NUM_THREADS=16; ./fem_openmp --gtest_filter=cpu.china2015_openmp_test | tee -a openmp_result.log
export OMP_NUM_THREADS=12; ./fem_openmp --gtest_filter=cpu.china2015_openmp_test | tee -a openmp_result.log
export OMP_NUM_THREADS=8;  ./fem_openmp --gtest_filter=cpu.china2015_openmp_test | tee -a openmp_result.log
export OMP_NUM_THREADS=6;  ./fem_openmp --gtest_filter=cpu.china2015_openmp_test | tee -a openmp_result.log
export OMP_NUM_THREADS=4;  ./fem_openmp --gtest_filter=cpu.china2015_openmp_test | tee -a openmp_result.log
export OMP_NUM_THREADS=2;  ./fem_openmp --gtest_filter=cpu.china2015_openmp_test | tee -a openmp_result.log
zip -r $(date -d "today" +"%Y%m%d%H%M")"_onmp_log" openmp_result.log fem_openmp 
rm ./fem_openmp
rm openmp_result.log

./fem_icpc --gtest_filter=cpu.china2015_seq_test | tee -a icpc_result.log
zip -r $(date -d "today" +"%Y%m%d%H%M")"_icpc_log" icpc_result.log fem_icpc
rm ./fem_icpc
rm icpc_result.log