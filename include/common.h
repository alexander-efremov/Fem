#ifndef COMMON_H_INCLUDED
#define COMMON_H_INCLUDED
#include <iostream>
#include <math.h>

struct ComputeParameters {
public:

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

    ComputeParameters() {
        const int x_length_default = 10;
        const int y_length_default = 10;
        t_count = 0;
        tau = 0.02;
        t_count = 50;
        a = 2.;
        b = 1.;
        lb = bb = 0.;
        rb = ub = 1.;
        double value = pow(2., 0);
        x_size = x_length_default * value;
        y_size = y_length_default * value;
        size = (x_size + 1) * (y_size + 1);
        tau /= value;
        t_count *= value;
    }

    ~ComputeParameters() {
    }

    int x_length() {
        return x_size + 1;
    }

    int y_length() {
        return y_size + 1;
    }

    // получает размер внутренней матрицы

    int get_inner_matrix_size() {
        return (x_length() - 2) * (y_length() - 2);
    }

    int get_real_matrix_size() {
        return x_length() * y_length();
    }

    // получает размер внутренней матрицы

    int get_inner_x_size() {
        return x_length() - 2;
    }

    friend std::ostream &operator<<(std::ostream &output,
            const ComputeParameters &tr) {

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
