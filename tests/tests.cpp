#include <iostream>
#include <gtest/gtest.h>
#include "timer.h"
#include "common.h"
#include "utils.h"

class cpu : public testing::Test
{
protected:

    double *GetCpuToLevel(ComputeParameters *p)
    {
        return solve_cpu_test(p->a, p->b, p->lb, p->rb, p->bb,
                              p->ub, p->tau, p->t_count, p->x_size,
                              p->y_size, p->level);
    }
};

TEST_F(cpu, main_test)
{
    const int finishLevel = 1;
    const int startLevel = 0;

    for (int level = startLevel; level < finishLevel; ++level)
    {
        ComputeParameters *p = new ComputeParameters(level);
        std::cout << *p << std::endl;
        double *data = GetCpuToLevel(p);
        print_matrix(p->get_real_x_size(), p->get_real_y_size(), data);

        delete p;
        delete[] data;
    }
}