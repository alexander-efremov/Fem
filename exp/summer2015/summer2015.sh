#!/bin/bash  

rm ./fem_openmp 2> /dev/null ./openmp_result.log 2> /dev/null
make -C ~/fem/ -f ~/fem/build/makefile.openmp clean
make -C ~/fem/ -f ~/fem/build/makefile.openmp build
cp ~/fem/fem_openmp ./fem_openmp

rm ./fem_icpc 2> /dev/null ./icpc_result.log 2> /dev/null
make -C ~/fem/ -f ~/fem/build/makefile.icpc clean
make -C ~/fem/ -f ~/fem/build/makefile.icpc build
cp ~/fem/fem_icpc ./fem_icpc

rm ./fem_cuda 2> /dev/null ./cuda_result.log 2> /dev/null
make -C ~/fem/ -f ~/fem/build/makefile.cuda clean
make -C ~/fem/ -f ~/fem/build/makefile.cuda build
cp ~/fem/fem_cuda ./fem_cuda

#./fem_icpc --gtest_filter=cpu.cpu_2_norm_test | tee -a icpc_result.log
#zip -r $(date -d "today" +"%Y%m%d%H%M")"_icpc_log" icpc_result.log fem_icpc
#rm ./fem_icpc ./icpc_result.log

export OMP_NUM_THREADS=24; ./fem_openmp --gtest_filter=cpu.omp_2_time_test | tee -a openmp_result.log
export OMP_NUM_THREADS=20; ./fem_openmp --gtest_filter=cpu.omp_2_time_test | tee -a openmp_result.log
export OMP_NUM_THREADS=16; ./fem_openmp --gtest_filter=cpu.omp_2_time_test | tee -a openmp_result.log
export OMP_NUM_THREADS=12; ./fem_openmp --gtest_filter=cpu.omp_2_time_test | tee -a openmp_result.log
export OMP_NUM_THREADS=8;  ./fem_openmp --gtest_filter=cpu.omp_2_time_test | tee -a openmp_result.log
export OMP_NUM_THREADS=6;  ./fem_openmp --gtest_filter=cpu.omp_2_time_test | tee -a openmp_result.log
export OMP_NUM_THREADS=4;  ./fem_openmp --gtest_filter=cpu.omp_2_time_test | tee -a openmp_result.log
export OMP_NUM_THREADS=2;  ./fem_openmp --gtest_filter=cpu.omp_2_time_test | tee -a openmp_result.log
zip -r $(date -d "today" +"%Y%m%d%H%M")"_onmp_log" openmp_result.log fem_openmp 
rm ./fem_openmp ./openmp_result.log

#./fem_cuda --gtest_filter=cpu.gpu_2_time_test | tee -a cuda_result.log
#zip -r $(date -d "today" +"%Y%m%d%H%M")"_cuda_log" cuda_result.log fem_cuda
#rm ./fem_cuda ./cuda_result.log