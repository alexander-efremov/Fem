#include <gtest/gtest.h>
#include "common.h"
#include "utils.h"
#include "LowOrdOper.h"

class cpu : public testing::Test
{
protected:

    double *solve(ComputeParameters *p, int lvl)
    {
        return cpu_solve(p->a, p->b, p->lb, p->rb, p->bb,
                         p->ub, p->tau, p->t_count, p->x_size,
                         p->y_size, lvl);
    }
    
    double* get_model_result(ComputeParameters *p, int lvl)
    {
  //      return NULL;
       return solByEqualVolWithVarStepPlusPrint(p->a, p->b, p->lb, p->rb, p->bb,
                         p->ub, p->tau, p->t_count, p->x_size,
                         p->y_size, true, true, lvl+1);
    }
};

TEST_F(cpu, DISABLED_main_test)
{
    const int first = 0;
    const int last = 1;
    ComputeParameters *p = new ComputeParameters();
    
    for (int lvl = first; lvl < last; ++lvl)
    {
        printf("level = %d\n", lvl);
        p->recompute_params(lvl);
        double *data = solve(p, lvl);
        print_matrix(data, p->x_length(), p->y_length());
        delete[] data;
    }
    delete p;
}

TEST_F(cpu, test_to_model)
{
    const int first = 0;
    const int last = 1;
    ComputeParameters *p = new ComputeParameters();
    
    for (int lvl = first; lvl < last; ++lvl)
    {
        printf("level = %d\n", lvl);
        p->recompute_params(lvl);
        double *data = solve(p, lvl);
        double *model = get_model_result(p, lvl);
        for (int i = 0; i < p->get_size(); i++)
        {
            EXPECT_NEAR(model[i], data[i], 1e-12);
        }
        
       // print_matrix(data, p->x_length(), p->y_length());
        delete[] data;
        delete[] model;
    }
    delete p;
}