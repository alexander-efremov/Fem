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

    void print_parameters(ComputeParameters *tr)
    {
        printf("a = %f\n", tr->a);
        printf("b = %f\n", tr->b);
        printf("lb = %f\n", tr->lb);
        printf("rb = %f\n", tr->rb);
        printf("bb = %f\n", tr->bb);
        printf("ub = %f\n", tr->ub);
        printf("size = %d\n", tr->size);
        printf("tau = %f\n", tr->tau);
        printf("Time levels = %d\n", tr->t_count);
        printf("x size = %d\n", tr->x_size);
        printf("y size = %d\n", tr->y_size);
    }
};

TEST_F(cpu, main_test)
{
    const int first = 0;
    const int last = 1;

    for (int lvl = first; lvl < last; ++lvl)
    {
        printf("level = %d\n", lvl);
        ComputeParameters *p = new ComputeParameters();
        print_parameters(p);
        double *data = solve(p, lvl);
        print_matrix(data, p->x_length(), p->y_length());
        delete p;
        delete[] data;
    }
}