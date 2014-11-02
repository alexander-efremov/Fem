#ifndef COMMON_H_INCLUDED
#define COMMON_H_INCLUDED
#include <iostream>
#include <iostream>
#include <math.h>

static const double C_pi = 3.14159265358979323846264338327;
static const double C_par_a = 2.; //   -  Solution parameter.
static const double C_par_b = 1.; //-  Item of second parameter from "u_funcion" or "v_funcion".
static const double C_StepsRel = 1. / 5.; //   -  A relation of the time step "C_tau" to the max grid step "maxOfGrSt";
static const double C_timeEnd = 1.; //   -  Finish time.
static const int C_numOfOXSt = 10; //   -  Number of OX steps (segments).
static const int C_numOfOYSt = 10; //   -  Number of OY steps (segments).
static double C_lbDom = 0., C_rbDom = 1.; //   -  Left and right boundaries of rectangular domain.
static double C_bbDom = 0., C_ubDom = 1.; //   -  Bottom and upper boundaries of rectangular domain.

struct ComputeParameters {
private:
     ComputeParameters() {
        t_count = 0;
        level = 0;
        tau = 0.02;
        t_count = 50;
        a = C_par_a;
        b = C_par_b;
        lb = C_lbDom;
        rb = C_rbDom;
        bb = C_bbDom;
        ub = C_ubDom;
    }
public:
    int level;
    double a;
    double b;
    double lb;
    double rb;
    double bb;
    double ub;
    double tau;
    int size;
    int t_count;
    int x_size;
    int y_size;

    ComputeParameters(int lvl) : ComputeParameters() {
        level = lvl;
        double value = pow(2., level);
        x_size = C_numOfOXSt * value;
        y_size = C_numOfOYSt * value;
        size = (x_size + 1) * (y_size + 1);
        tau /= value;
        t_count *= value;
    }

    ~ComputeParameters() {
    }
    
    int x_size() {
        return x_size + 1;
    }

    int y_size() {
        return y_size + 1;
    }

    // получает размер внутренней матрицы

    int get_inner_matrix_size() {
        return (x_size() - 2) * (y_size() - 2);
    }

    int get_real_matrix_size() {
        return x_size() * y_size();
    }

    // получает размер внутренней матрицы

    int get_inner_x_size() {
        return x_size() - 2;
    }

    friend std::ostream &operator<<(std::ostream &output,
            const ComputeParameters &tr) {
        output << "level = " << tr.level << std::endl;
        output << "a = " << tr.a << std::endl;
        output << "b = " << tr.b << std::endl;
        output << "lb = " << tr.lb << std::endl;
        output << "rb = " << tr.rb << std::endl;
        output << "bb = " << tr.bb << std::endl;
        output << "ub = " << tr.ub << std::endl;
        output << "size = " << tr.size << std::endl;
        output << "tau = " << tr.tau << std::endl;
        output << "Time levels = " << tr.t_count << std::endl;
        output << "x size = " << tr.x_size << std::endl;
        output << "y size = " << tr.y_size << std::endl;
        return output;
    }
};


extern double *cpu_solve(
        double a,
        double b,
        double lb,
        double rb,
        double bb,
        double ub,
        double time_step,
        int time_step_count,
        int ox_length,
        int oy_length,
        const int step);

#endif
