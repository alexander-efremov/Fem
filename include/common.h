#ifndef COMMON_H_INCLUDED
#define COMMON_H_INCLUDED
#include <iostream>
#include <iostream>
#include <math.h>

static const double C_pi = 3.14159265358979323846264338327;
static const double C_par_a = 2.; //   -  Solution parameter.
static const int C_numOfOXSt = 10; //   -  Number of OX steps (segments).
static const double C_par_b = 1.; //-  Item of second parameter from "u_funcion" or "v_funcion".
static const int C_numOfOYSt = 10; //   -  Number of OY steps (segments).
static const double C_StepsRel = 1. / 5.; //   -  A relation of the time step "C_tau" to the max grid step "maxOfGrSt";
static const double C_timeEnd = 1.; //   -  Finish time.
static double C_lbDom = 0., C_rbDom = 1.; //   -  Left and right boundaries of rectangular domain.
static double *masOX; //   -  Massive of OX points. Dimension = C_numOfOXSt +1.
static double C_bbDom = 0., C_ubDom = 1.; //   -  Botton and upper boundaries of rectangular domain.
static double *masOY; //   -  Massive of OY points. Dimension = C_numOfOYSt +1.
static double C_tau; //   -  Time step. Will be computed by relation "C_StepsRel" with small correction.
static int C_numOfTSt; //   -  A number of time steps.

inline void initCompOfGlVar() {
    double maxOfGrSt = 0.;
    int j;

    //   Massive of OX steps. Dimension = C_numOfOXSt +1.
    masOX = new double[ C_numOfOXSt + 1];
    for (j = 0; j < C_numOfOXSt + 1; j++) {
        //   Regular grid for simplification.
        masOX[j] = C_lbDom + (C_rbDom - C_lbDom) * ((double) j / C_numOfOXSt);
    }

    //   Let's find out maximum of OX step.
    for (j = 0; j < C_numOfOXSt; j++) {
        if (maxOfGrSt < (masOX[j + 1] - masOX[j])) {
            maxOfGrSt = (masOX[j + 1] - masOX[j]);
        }
    }

    //   Massive of OY steps. Dimension = C_numOfOYSt +1.
    masOY = new double[ C_numOfOYSt + 1];
    for (j = 0; j < C_numOfOYSt + 1; j++) {
        //   Regular grid for simplification.
        masOY[j] = C_bbDom + (C_ubDom - C_bbDom) * ((double) j / C_numOfOYSt);
    }

    //   Let's find out maximum of OY step.
    for (j = 0; j < C_numOfOYSt; j++) {
        if (maxOfGrSt < (masOY[j + 1] - masOY[j])) {
            maxOfGrSt = (masOY[j + 1] - masOY[j]);
        }
    }

    C_tau = C_StepsRel * maxOfGrSt;

    C_numOfTSt = (int) (C_timeEnd / C_tau);


    if (fabs(C_numOfTSt * C_tau - C_timeEnd) < (C_tau * 1.e-7)) {
        //   O.K. This error is error of computer arithmetic.

        return;
    }

    if (fabs(C_numOfTSt * C_tau - C_timeEnd) >= (C_tau * 1.e-7)) {
        //   We need to change number of time steps and time step.

        //cout<<"\nA number of time steps and time step itself have changed."<< flush;

        C_numOfTSt++;

        C_tau = C_timeEnd / C_numOfTSt;

        //cout<<"\nNow time step = "<< C_tau << " and a number of such time steps = " << C_numOfTSt << flush;

        //cout<<"\nSo finish time is "<< C_tau*C_numOfTSt <<flush;

        return;
    }
}

inline void memClean() {
    delete[] masOX;
    delete[] masOY;
}

struct ComputeParameters {
private:
    int time_i;
    bool _initresult;
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
    int i;
    int j;
    double *result;
    double currentTimeLevel;
    double *diff;
    bool _initdiff;

    ComputeParameters(int level, bool initresult, bool initdiff = false) : currentTimeLevel(1), t_count(0) {
        _initresult = initresult;
        _initdiff = initdiff;

        a = C_par_a;
        b = C_par_b;
        lb = C_lbDom;
        rb = C_rbDom;
        bb = C_bbDom;
        ub = C_ubDom;

        double value = pow(2., level);
        int n = C_numOfOXSt * value;
        tau = C_tau / value;
        t_count = C_numOfTSt * value;

        x_size = n;

        int m = C_numOfOYSt * value;

        y_size = m;
        size = (n + 1) * (m + 1);
        result = NULL;
        if (initresult) {
            result = new double[size];
        }
        if (_initdiff) {
            diff = new double[size];
        }

    }

    ~ComputeParameters() {
        if (_initresult)
            delete[] result;
        if (_initdiff)
            delete[] diff;
    }

    void reset_time_counter() {
        time_i = 1;
    }

    bool can_iterate_over_time_level() {
        return time_i <= t_count;
    }

    void inc_time_level() {
        time_i++;
    }

    int get_real_x_size() {
        return x_size + 1;
    }

    int get_real_y_size() {
        return y_size + 1;
    }

    // получает размер внутренней матрицы

    int get_inner_matrix_size() {
        return (get_real_x_size() - 2) * (get_real_y_size() - 2);
    }

    int get_real_matrix_size() {
        return get_real_x_size() * get_real_y_size();
    }

    // получает размер внутренней матрицы

    int get_inner_x_size() {
        return get_real_x_size() - 2;
    }

    void print_info() {
        std::cout << "current time level " << currentTimeLevel << std::endl;
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


extern double *solve_cpu_test(
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
