#include <iostream>
#include <gtest/gtest.h>
#include "common.h"
#include "utils.h"

class cpu : public testing::Test
{
protected:
    double *solve(ComputeParameters *p)
    {
        return cpu_solve(p->a, p->b, p->lb, p->rb, p->bb,
                              p->ub, p->tau, p->t_count, p->x_size,
                              p->y_size, p->level);
    }
};

TEST_F(cpu, main_test)
{
    const int first = 0;
    const int last = 1;
    
    for (int lvl = first; lvl < last; ++lvl)
    {
        ComputeParameters *p = new ComputeParameters(lvl);
        std::cout << *p << std::endl;
        double *data = solve(p);
        print_matrix(p->x_size(), p->y_size(), data);
        delete p;
        delete[] data;
    }
}