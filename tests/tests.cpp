#include <iostream>
#include <gtest/gtest.h>
#include "timer.h"
#include "common.h"
#include "utils.h"

class cpu : public testing::Test {
protected:
    double *GetCpuToLevel(int level) {
        double tau = 0.02;
        double time_step_count = 50;
        return solve_cpu_test(C_par_a, C_par_b, C_lbDom, C_rbDom, C_bbDom,
                C_ubDom, tau, time_step_count,  C_numOfOXSt,
                C_numOfOYSt, level);
    }
};

TEST_F(cpu, main_test) {
    const int finishLevel = 1;
    const int startLevel = 0;

    for (int level = startLevel; level < finishLevel; ++level) {
        std::cout << "level = " << level << std::endl;
        ComputeParameters *p = new ComputeParameters(level, true, false);
        std::cout << *p << std::endl;
        double *data = GetCpuToLevel(level);
        print_matrix(p->get_real_x_size(), p->get_real_y_size(), data);

        delete p;
        delete[] data;
    }
}