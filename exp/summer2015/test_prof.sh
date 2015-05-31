#!/bin/bash

rm ./fem_cuda 2> /dev/null ./cuda_result.log 2> /dev/null
make -C ~/fem/ -f ~/fem/build/makefile.cuda clean
make -C ~/fem/ -f ~/fem/build/makefile.cuda build
cp ~/fem/fem_cuda ./fem_cuda

./fem_cuda --gtest_filter=cpu.gpu_2_time_test_1 | tee -a cuda_result.log
zip -r $(date -d "today" +"%Y%m%d%H%M")"_cuda_log" cuda_result.log fem_cuda
rm ./fem_cuda ./cuda_result.log
