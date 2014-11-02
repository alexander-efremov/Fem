#include <iostream>
#include <gtest/gtest.h>
#include "common.h"
#include "utils.h"

class cpu : public testing::Test
{
protected:
    double *solve(ComputeParameters *p, int lvl)
    {
        return cpu_solve(p->a, p->b, p->lb, p->rb, p->bb,
                         p->ub, p->tau, p->t_count, p->x_size,
                         p->y_size, lvl);
    }
};

TEST_F(cpu, main_test)
{
    const int first = 0;
    const int last = 1;

    for (int lvl = first; lvl < last; ++lvl)
    {
        std::cout << "level = " << lvl << std::endl;
        ComputeParameters *p = new ComputeParameters();
        std::cout << *p << std::endl;
        double *data = solve(p, lvl);
        print_matrix(data, p->x_length(), p->y_length());
        delete p;
        delete[] data;
    }
}