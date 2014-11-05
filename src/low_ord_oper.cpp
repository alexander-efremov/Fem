#include "common.h"
#include "point.h"
#include "utils.h"

static int TMP_WALL_CNT = 0;
static const double _PI = 3.14159265358979323846264338327;
static const double _MIN_VALUE = 1.e-12;
static const double _MIN_VALUE_1 = 1.e-14;
static double B;
static double UB;
static double BB;
static double LB;
static double RB;
static double TAU;
static int OX_LEN;
static int OY_LEN;
static int XY_LEN;
static int TIME_STEP_CNT;

double analytical_solution(double t, double x, double y) {
    return 1.1 + sin(t * x * y);
}

double init_side(double x, double y, double t, bound bound) {
    switch (bound) {
        case up:
            return analytical_solution(t, x, y);
        case bottom:
            return analytical_solution(t, x, y);
        case left:
            return analytical_solution(t, x, y);
        case right:
            return analytical_solution(t, x, y);
    }
    return 0;
}

double integrate_first_type(double py,
        double qy,
        double gx,
        double hx,
        double a,
        double b) {
    double integ = (hx - a) * (hx - a) - (gx - a) * (gx - a);
    integ *= (qy - b) * (qy - b) - (py - b) * (py - b);
    return integ / 4;
}

double integrate_second_type(double py,
        double qy,
        double alpha,
        double a,
        double b,
        double beta) {
    double tmp, integ;
    tmp = (qy - alpha) * (a * qy + b - beta) * (a * qy + b - beta) * (a * qy + b - beta);
    tmp -= (py - alpha) * (a * py + b - beta) * (a * py + b - beta) * (a * py + b - beta);
    integ = tmp / (3 * a);
    tmp = (a * qy + b - beta) * (a * qy + b - beta) * (a * qy + b - beta) * (a * qy + b - beta);
    tmp -= (a * py + b - beta) * (a * py + b - beta) * (a * py + b - beta) * (a * py + b - beta);
    return integ - tmp / (12 * a * a);
}

double integrate_rectangle_one_cell(double Py,
        double Qy,
        double Gx,
        double Hx,
        int tl,
        int *indCurSqOx, //   -  Index of current square by Ox axis.
        int *indCurSqOy,
        const double *ox,
        const double *oy,
        double *density) {
    double hx = ox[1] - ox[0];
    double hy = oy[1] - oy[0];
    double result = 0;
    double tmp;
    double rho[2][2];
    double t = TAU * (tl - 1.);
    double x, y;
    if (indCurSqOx[0] >= 0 && indCurSqOy[0] >= 0) {
        rho[0][0] = density[ (OX_LEN + 1) * indCurSqOy[0] + indCurSqOx[0] ];
        rho[0][1] = density[ (OX_LEN + 1) * indCurSqOy[1] + indCurSqOx[0] ];
        rho[1][0] = density[ (OX_LEN + 1) * indCurSqOy[0] + indCurSqOx[1] ];
        rho[1][1] = density[ (OX_LEN + 1) * indCurSqOy[1] + indCurSqOx[1] ];
    } else {
        // TODO: убрать потому что это неверно (надо расчитывать граничные условия)
        x = indCurSqOx[0] * hx;
        y = indCurSqOy[0] * hy;
        rho[0][0] = analytical_solution(t, x, y);
        x = indCurSqOx[0] * hx;
        y = indCurSqOy[1] * hy;
        rho[0][1] = analytical_solution(t, x, y);
        x = indCurSqOx[1] * hx;
        y = indCurSqOy[0] * hy;
        rho[1][0] = analytical_solution(t, x, y);
        x = indCurSqOx[1] * hx;
        y = indCurSqOy[1] * hy;
        rho[1][1] = analytical_solution(t, x, y);

        TMP_WALL_CNT++;

    }

    if (indCurSqOx[1] >= 0 && indCurSqOy[1] >= 0) {
        tmp = integrate_first_type(Py, Qy, Gx, Hx, ox[ indCurSqOx[1] ], oy[ indCurSqOy[1] ]);
    } else {
        tmp = integrate_first_type(Py, Qy, Gx, Hx, hx * indCurSqOx[1], hy * indCurSqOy[1]);
    }
    tmp = tmp / hx / hy;
    result = tmp * rho[0][0]; //   rhoInPrevTL[ indCurSqOx[0] ][ indCurSqOy[0] ];
    if (indCurSqOx[0] >= 0 && indCurSqOy[1] >= 0) {
        tmp = integrate_first_type(Py, Qy, Gx, Hx, ox[ indCurSqOx[0] ], oy[ indCurSqOy[1] ]);
    } else {
        tmp = integrate_first_type(Py, Qy, Gx, Hx, hx * indCurSqOx[0], hy * indCurSqOy[1]);
    }
    tmp = tmp / hx / hy;
    result = result - tmp * rho[1][0]; //   rhoInPrevTL[ indCurSqOx[1] ][ indCurSqOy[0] ];
    if (indCurSqOx[1] >= 0 && indCurSqOy[0] >= 0) {
        tmp = integrate_first_type(Py, Qy, Gx, Hx, ox[ indCurSqOx[1] ], oy[ indCurSqOy[0] ]);
    } else {
        tmp = integrate_first_type(Py, Qy, Gx, Hx, hx * indCurSqOx[1], hy * indCurSqOy[0]);
    }
    tmp = tmp / hx / hy;
    result -= tmp * rho[0][1]; //   rhoInPrevTL[ indCurSqOx[0] ][ indCurSqOy[1] ];
    if (indCurSqOx[0] >= 0 && indCurSqOy[0] >= 0) {
        tmp = integrate_first_type(Py, Qy, Gx, Hx, ox[ indCurSqOx[0] ], oy[ indCurSqOy[0] ]);
    } else {
        tmp = integrate_first_type(Py, Qy, Gx, Hx, hx * indCurSqOx[0], hy * indCurSqOy[0]);
    }
    tmp = tmp / hx / hy;
    return result + tmp * rho[1][1]; //   rhoInPrevTL[ indCurSqOx[1] ][ indCurSqOy[1] ];
}

double integrate_triangle_left_one_cell(
        double Py,
        double Qy,
        double a_SL,
        double b_SL,
        double Hx,
        int tl,
        int *indCurSqOx, //   -  Index of current square by Ox axis.
        int *indCurSqOy, //   -  Index of current square by Oy axis.        
        const double *ox,
        const double *oy,
        double *density) {
    double hx = ox[1] - ox[0];
    double hy = oy[1] - oy[0];
    double result = 0;
    double tmp, bufInteg_D;
    double rho[2][2];
    double t = TAU * (tl - 1.);
    double x, y;
    if (indCurSqOx[0] >= 0 && indCurSqOx[1] <= OX_LEN) {
        if (indCurSqOy[0] >= 0 && indCurSqOy[1] <= OY_LEN) {
            rho[0][0] = density[ (OX_LEN + 1) * indCurSqOy[0] + indCurSqOx[0] ];
            rho[0][1] = density[ (OX_LEN + 1) * indCurSqOy[1] + indCurSqOx[0] ];
            rho[1][0] = density[ (OX_LEN + 1) * indCurSqOy[0] + indCurSqOx[1] ];
            rho[1][1] = density[ (OX_LEN + 1) * indCurSqOy[1] + indCurSqOx[1] ];
        }
    }

    // TODO: убрать потому что это неверно (надо расчитывать граничные условия)
    // норма должна уменьшиться
    if (indCurSqOx[0] < 0 || indCurSqOx[1] > OX_LEN || indCurSqOy[0] < 0 || indCurSqOy[1] > OY_LEN) {
        x = indCurSqOx[0] * hx;
        y = indCurSqOy[0] * hy;
        rho[0][0] = analytical_solution(t, x, y);
        x = indCurSqOx[0] * hx;
        y = indCurSqOy[1] * hy;
        rho[0][1] = analytical_solution(t, x, y);
        x = indCurSqOx[1] * hx;
        y = indCurSqOy[0] * hy;
        rho[1][0] = analytical_solution(t, x, y);
        x = indCurSqOx[1] * hx;
        y = indCurSqOy[1] * hy;
        rho[1][1] = analytical_solution(t, x, y);

        TMP_WALL_CNT++;
    }

    //   1.
    tmp = (Qy - oy[ indCurSqOy[1] ]) * (Qy - oy[ indCurSqOy[1] ]) - (Py - oy[ indCurSqOy[1] ]) * (Py - oy[ indCurSqOy[1] ]);
    if ((indCurSqOx[1] >= 0) && (indCurSqOy[1] >= 0)) {
        tmp = tmp * (Hx - ox[ indCurSqOx[1] ]) * (Hx - ox[ indCurSqOx[1] ]) / 4.;
        bufInteg_D = integrate_second_type(Py, Qy, oy[ indCurSqOy[1] ], a_SL, b_SL, ox[ indCurSqOx[1] ]);
    } else {
        tmp = tmp * (Hx - hx * indCurSqOx[1]) * (Hx - hx * indCurSqOx[1]) / 4.;
        bufInteg_D = integrate_second_type(Py, Qy, hy * indCurSqOy[1], a_SL, b_SL, hx * indCurSqOx[1]);
    }
    tmp -= bufInteg_D / 2.;
    result = tmp * rho[0][0] / hx / hy;

    //   2.
    tmp = (Qy - oy[ indCurSqOy[1] ]) * (Qy - oy[ indCurSqOy[1] ]) - (Py - oy[ indCurSqOy[1] ]) * (Py - oy[ indCurSqOy[1] ]);
    if ((indCurSqOx[0] >= 0) && (indCurSqOy[1] >= 0)) {
        tmp = -1. * tmp * (Hx - ox[ indCurSqOx[0] ]) * (Hx - ox[ indCurSqOx[0] ]) / 4.;
        bufInteg_D = integrate_second_type(Py, Qy, oy[ indCurSqOy[1] ], a_SL, b_SL, ox[ indCurSqOx[0] ]);
    } else {
        tmp = -1. * tmp * (Hx - hx * indCurSqOx[0]) * (Hx - hx * indCurSqOx[0]) / 4.;
        bufInteg_D = integrate_second_type(Py, Qy, hy * indCurSqOy[1], a_SL, b_SL, hx * indCurSqOx[0]);
    }
    tmp = tmp + bufInteg_D / 2.;
    result += tmp * rho[1][0] / hx / hy;

    //   3.
    tmp = (Qy - oy[ indCurSqOy[0] ]) * (Qy - oy[ indCurSqOy[0] ]) - (Py - oy[ indCurSqOy[0] ]) * (Py - oy[ indCurSqOy[0] ]);
    if ((indCurSqOx[1] >= 0) && (indCurSqOy[0] >= 0)) {
        tmp = -1. * tmp * (Hx - ox[ indCurSqOx[1] ]) * (Hx - ox[ indCurSqOx[1] ]) / 4.;
        bufInteg_D = integrate_second_type(Py, Qy, oy[ indCurSqOy[0] ], a_SL, b_SL, ox[ indCurSqOx[1] ]);
    } else {
        tmp = -1. * tmp * (Hx - hx * indCurSqOx[1]) * (Hx - hx * indCurSqOx[1]) / 4.;
        bufInteg_D = integrate_second_type(Py, Qy, hy * indCurSqOy[0], a_SL, b_SL, hx * indCurSqOx[1]);
    }
    tmp = tmp + bufInteg_D / 2.;
    result += tmp * rho[0][1] / hx / hy;

    //   4.
    tmp = (Qy - oy[ indCurSqOy[0] ]) * (Qy - oy[ indCurSqOy[0] ]) - (Py - oy[ indCurSqOy[0] ]) * (Py - oy[ indCurSqOy[0] ]);
    if ((indCurSqOx[0] >= 0) && (indCurSqOy[0] >= 0)) {
        tmp = tmp * (Hx - ox[ indCurSqOx[0] ]) * (Hx - ox[ indCurSqOx[0] ]) / 4.;
        bufInteg_D = integrate_second_type(Py, Qy, oy[ indCurSqOy[0] ], a_SL, b_SL, ox[ indCurSqOx[0] ]);
    } else {
        tmp = tmp * (Hx - hx * indCurSqOx[0]) * (Hx - hx * indCurSqOx[0]) / 4.;
        bufInteg_D = integrate_second_type(Py, Qy, hy * indCurSqOy[0], a_SL, b_SL, hx * indCurSqOx[0]);
    }
    tmp -= bufInteg_D / 2.;
    result += tmp * rho[1][1] / hx / hy;

    return result;
}

double integrate_triangle_right_one_cell(double py,
        double qy,
        double a_SL,
        double b_SL,
        double gx,
        int tl,
        int *indCurSqOx, //   -  Index of current square by Ox axis.
        int *indCurSqOy, //   -  Index of current square by Oy axis.
        const double *ox,
        const double *oy,
        double *density) {
    return -1. * integrate_triangle_left_one_cell(
            py, qy,
            a_SL, b_SL, gx,
            tl,
            indCurSqOx, //   -  Index of current square by Ox axis.
            indCurSqOy, //   -  Index of current square by Oy axis.            
            ox, oy,
            density);
}

double integrate_chanel_slant_right(int tl,
        point_t &bv, int wTrPCI, //   -  Where travel point current (botton vertex) is.
        point_t &uv, int wTrPNI, //   -  Where travel point next (upper vertex) is.
        //
        int *indCurSqOx, //   -  Index by OX axis where bv and uv are.
        //
        double lb, int *indLB, //   -  Left boundary by Ox. Index by OX axis where lb is.
        //
        int *indCurSqOy, //   -  Index of current square by Oy axis.
        //
        const double *ox,
        const double *oy,
        double *density) {
    double mv[2], rv[2]; //   -  Middle and right vertices.
    int wMvI; //   -  Where middle vertex is.
    int indCurSqOxToCh[2]; //   -  Indices of current square by Ox axis to be changed. Under which we want to integrate.
    double h = ox[1] - ox[0];
    double a_SL, b_SL; //   -  Coefficients of slant line: x = a_SL *y  +  b_SL.
    double Gx, Hx; //   -  Left boundary for each integration.
    double result = 0.;
    double tmp;

    //   Let's compute helpful values.

    if (uv[0] <= bv[0]) {
        mv[0] = uv[0];
        mv[1] = uv[1];
        wMvI = wTrPNI;
        rv[0] = bv[0];
        rv[1] = bv[1];
    }

    if (uv[0] > bv[0]) {
        mv[0] = bv[0];
        mv[1] = bv[1];
        wMvI = wTrPCI;
        rv[0] = uv[0];
        rv[1] = uv[1];
    }

    if ((fabs(uv[1] - bv[1])) <= _MIN_VALUE) {
        //   Computation is impossible. Too smale values. Let's return some approximate value.
        //   buf_D  =  (uv[1] - bv[1])  *  ((uv[0] + bv[0]) /2.  -  lb) * rhoInPrevTL[ indCurSqOx[0] ][ indCurSqOy[0] ];
        return fabs(uv[1] - bv[1]); //   fabs(uv[1] - bv[1]);
    }


    //   First step: from "lb" to "masOX[ indCurSqOx[0] ]" by iteration.
    //   integ  += fabs( mv[0] - lb) * fabs(uv[1] - bv[1]);

    indCurSqOxToCh[0] = indLB[0];
    indCurSqOxToCh[1] = indCurSqOxToCh[0] + 1;

    for (int j = indLB[0]; j < indCurSqOx[0]; j++) {
        //   If this is first cell we should integrate under rectangle only.
        if (indCurSqOxToCh[0] >= 0) {
            Gx = ox[ indCurSqOxToCh[0] ];
            Hx = ox[ indCurSqOxToCh[1] ];
        }


        if (indCurSqOxToCh[0] < 0) {
            Gx = h * indCurSqOxToCh[0];
            Hx = h * indCurSqOxToCh[1];
        }

        if (j == indLB[0]) {
            Gx = lb;
        }

        tmp = integrate_rectangle_one_cell(
                bv[1], //   -  double Py,
                uv[1], //   -  double Qy,
                Gx, //   -  double Gx,
                Hx, //   -  double Hx,
                //
                tl,
                indCurSqOxToCh, //   -  Index of current square by Ox axis.
                indCurSqOy, //   -  Index of current square by Oy axis.
                //
                ox,
                oy,
                density);

        result += tmp;

        indCurSqOxToCh[0] += 1;
        indCurSqOxToCh[1] = indCurSqOxToCh[0] + 1;
    }

    //   Integration. Second step: under [ indCurSqOx[0]; indCurSqOx[1] ] square.

    //   A. Under rectangle.
    if (wMvI == 1) {
        if (indCurSqOx[0] == indLB[0]) {
            Gx = lb;
        }

        if (indCurSqOx[0] > indLB[0]) {
            if (indCurSqOx[0] >= 0) {
                Gx = ox[ indCurSqOx[0] ];
            }

            if (indCurSqOx[0] < 0) {
                Gx = h * indCurSqOx[0];
            }
        }

        tmp = integrate_rectangle_one_cell(bv[1], //   -  double Py,
                uv[1], //   -  double Qy,
                //
                Gx, //   -  double Gx,
                mv[0], //   -  double Hx,
                //
                tl,
                indCurSqOx, //   -  Index of current square by Ox axis.
                indCurSqOy, //   -  Index of current square by Oy axis.
                //
                ox,
                oy,
                density);

        result += tmp;

    }

    //   B. Under triangle.

    if (fabs(uv[1] - bv[1]) > _MIN_VALUE) {
        //   integ += fabs(uv[1] - bv[1]) * (rv[0] - mv[0]) /2.;
        //   Coefficients of slant line: x = a_SL *y  +  b_SL.
        a_SL = (uv[0] - bv[0]) / (uv[1] - bv[1]);
        b_SL = bv[0] - a_SL * bv[1];


        //   Integration under one cell triangle.

        if (fabs(a_SL) > _MIN_VALUE) {
            tmp = integrate_triangle_right_one_cell(
                    bv[1], //   -  double Py,
                    uv[1], //   -  double Qy,
                    //
                    a_SL,
                    b_SL,
                    mv[0], //   -  double Gx,
                    //
                    tl,
                    indCurSqOx, //   -  Index of current square by Ox axis.
                    indCurSqOy, //   -  Index of current square by Oy axis.
                    //
                    ox,
                    oy,
                    density);

            result += tmp;
        }
    }

    return result;
}

double integrate_chanel_slant_left(
        int tl,
        point_t &bv, int wTrPCI, //   -  Where travel point current (bottom vertex) is.
        point_t &uv, int wTrPNI, //   -  Where travel point next (upper vertex) is.
        //
        int *indCurSqOx, //   -  Index by OX axis where bv and uv are.
        //
        double rb, int *indRB, //   -  Right boundary by Ox. Index by OX axis where rb is.
        //
        int *indCurSqOy, //   -  Index of current square by Oy axis.
        //
        const double *ox,
        const double *oy,
        double *density) {
    double lv[2], mv[2]; //   -  Left and middle vertices.
    int wMvI; //   -  Where middle vertex is.
    int indCurSqOxToCh[2]; //   -  Indices of current square by Ox axis to be changed. Under which we want to integrate.
    double h = ox[1] - ox[0];
    double a_SL, b_SL; //   -  Coefficients of slant line: x = a_SL *y  +  b_SL.
    double gx, hx; //   -  Left and right boundary for each integration.
    double result = 0.;
    double tmp;
    int j;

    //   Let's compute helpful values.

    if (uv[0] <= bv[0]) {
        lv[0] = uv[0];
        lv[1] = uv[1];
        mv[0] = bv[0];
        mv[1] = bv[1];
        wMvI = wTrPCI;
    }

    if (uv[0] > bv[0]) {
        lv[0] = bv[0];
        lv[1] = bv[1];
        mv[0] = uv[0];
        mv[1] = uv[1];
        wMvI = wTrPNI;
    }

    if ((fabs(uv[1] - bv[1])) <= _MIN_VALUE) {
        //   Computation is impossible. Too smale values. Let's return some approximate value.
        //   buf_D  =  (uv[1] - bv[1])  *  (rb  - (uv[0] + bv[0]) /2.) * rhoInPrevTL[ indCurSqOx[0] ][ indCurSqOy[0] ];
        return fabs(uv[1] - bv[1]); //   fabs(uv[1] - bv[1]);
    }

    //   Integration. First step: under [ indCurSqOx[0]; indCurSqOx[1] ] square.

    //   A. Under triangle.

    if (fabs(uv[1] - bv[1]) > _MIN_VALUE) {
        //   Coefficients of slant line: x = a_SL *y  +  b_SL.
        a_SL = (uv[0] - bv[0]) / (uv[1] - bv[1]);
        b_SL = bv[0] - a_SL * bv[1];

        //   Integration under one cell triangle.
        if (fabs(a_SL) > _MIN_VALUE) {
            tmp = integrate_triangle_left_one_cell(
                    bv[1], //   -  double Py
                    uv[1], //   -  double Qy                    
                    a_SL, b_SL, mv[0], //   -  double Hx,                   
                    tl,
                    indCurSqOx, //   -  Index of current square by Ox axis.
                    indCurSqOy, //   -  Index of current square by Oy axis.                   
                    ox,
                    oy,
                    density);
            result += tmp;
        }
    }

    //   B. Under rectangle. Need to be checking.
    if (wMvI == 1) {
        if (indCurSqOx[0] == indRB[0]) {
            hx = rb;
        }

        if (indCurSqOx[0] < indRB[0]) {
            if (indCurSqOx[1] >= 0) {
                hx = ox[ indCurSqOx[1] ];
            }

            if (indCurSqOx[1] < 0) {
                hx = h * indCurSqOx[1];
            }
        }

        tmp = integrate_rectangle_one_cell(bv[1], //   -  double Py,
                uv[1], //   -  double Qy,                
                mv[0], //   -  double Gx,
                hx, //   -  double Hx,                
                tl,
                indCurSqOx, //   -  Index of current square by Ox axis.
                indCurSqOy, //   -  Index of current square by Oy axis.                
                ox, oy,
                density);

        result += tmp;
    }

    //   Second step: from "masOX[ indCurSqOx[1] ]" to "rb" by iteration.


    indCurSqOxToCh[0] = indCurSqOx[0] + 1;
    indCurSqOxToCh[1] = indCurSqOxToCh[0] + 1;

    for (j = indCurSqOx[0] + 1; j < indRB[0] + 1; j++) {
        //   If this is first cell we should integrate under triangle only.

        if (indCurSqOxToCh[1] > 0) {
            gx = ox[ indCurSqOxToCh[0] ];
            hx = ox[ indCurSqOxToCh[1] ];
        }


        if (indCurSqOxToCh[1] <= 0) {
            gx = h * indCurSqOxToCh[0];
            hx = h * indCurSqOxToCh[1];
        }


        if (j == indRB[0]) {
            hx = rb;
        }


        tmp = integrate_rectangle_one_cell(bv[1], //   -  double Py,
                uv[1], //   -  double Qy,
                //
                gx, //   -  double Gx,
                hx, //   -  double Hx,
                //
                tl, //   -  Index of current time layer.
                //
                indCurSqOxToCh, //   -  Index of current square by Ox axis.
                indCurSqOy, //   -  Index of current square by Oy axis.
                //
                ox, oy,
                density);

        result += tmp;

        indCurSqOxToCh[0] += 1;
        indCurSqOxToCh[1] = indCurSqOxToCh[0] + 1;
    }

    return result;
}

double integrate_right_triangle_bottom_left(
        point_t &bv,
        point_t &uv,
        int tl,
        const double *ox,
        const double *oy,
        double *density) {
    point_t trPC; //   -  Travel point current;
    int wTrPCI = 0; //   -  Where travel point current is?
    point_t trPN; //   -  Travel point next;
    int wTrPNI = 0; //   -  Where travel point next is?
    double ang; //   -  Angle of slant line. Should be greater zero.
    int indCurSqOx[2], indCurSqOy[2]; //   -  Index of current square by Ox and Oy axes.
    int indRB[2]; //   -  Index of right boundary.
    double distOx, distOy; //   -  Distance to near Ox and Oy straight lines.
    bool isTrDone = false; //   -  Is travel done.
    double hx = ox[1] - ox[0];
    double hy = oy[1] - oy[0];
    double result = 0.; //   -  Value which we are computing.
    double tmp;
    //   Initial data.
    trPC = bv;
    if ((fabs(bv[0] - uv[0])) < _MIN_VALUE) {
        //   This triangle has very small width. I guess further computation isn't correct.
        return fabs(bv[0] - uv[0]);
    }
    ang = (uv[1] - bv[1]) / (bv[0] - uv[0]);
    if (fabs(ang) < _MIN_VALUE) {
        //   This triangle has very small height. I guess further computation isn't correct.
        return fabs(ang);
    }
    indCurSqOx[0] = (int) ((trPC[0] - _MIN_VALUE_1) / hx); //   -  If trPC[0] is in grid edge I want it will be between in the left side of indCurSqOx[1].
    if ((trPC[0] - _MIN_VALUE_1) <= 0) {
        indCurSqOx[0] -= 1; //   -  The case when "trPC[0]" ia negative.
    }
    indCurSqOx[1] = indCurSqOx[0] + 1; //   -  It's important only in rare case then trPC is in grid edge.
    indRB[0] = indCurSqOx[0];
    indRB[1] = indRB[0] + 1;
    indCurSqOy[0] = (int) ((trPC[1] + _MIN_VALUE_1) / hy); //   -  If trPC[1] is in grid edge I want it will be between indCurSqOx[0] and indCurSqOx[1].
    if ((trPC[1] + _MIN_VALUE_1) <= 0) {
        indCurSqOy[0] -= 1; //   -  The case when "trPC[0]" ia negative.
    }
    indCurSqOy[1] = indCurSqOy[0] + 1; //   -  It's important only in rare case then trPC is in grid edge.
    if (indCurSqOx[0] >= 0) {
        distOx = trPC[0] - ox[ indCurSqOx[0] ];
    }
    if (indCurSqOx[0] < 0) {
        distOx = fabs(trPC[0] - hx * indCurSqOx[0]);
    }
    if (indCurSqOy[1] >= 0) {
        distOy = oy[ indCurSqOy[1] ] - trPC[1];
    }
    if (indCurSqOy[1] < 0) {
        distOy = fabs(hy * indCurSqOy[1] - trPC[1]);
    }
    do {
        //   a. First case.
        if ((distOy / distOx) <= ang) {
            //   Across with straight line parallel Ox axis.
            wTrPNI = 1;
            if (indCurSqOy[1] >= 0) {
                trPN.y = oy[ indCurSqOy[1] ];
            }
            if (indCurSqOy[1] < 0) {
                trPN.y = hy * indCurSqOy[1];
            }
            trPN.x = bv[0] - (trPN.y - bv[1]) / ang;
        }
        //   b. Second case.
        if ((distOy / distOx) > ang) {
            //   Across with straight line parallel Oy axis.
            wTrPNI = 2;
            if (indCurSqOx[0] >= 0) {
                trPN.x = ox[ indCurSqOx[0] ];
            }
            if (indCurSqOx[0] < 0) {
                trPN.x = hx * indCurSqOx[0];
            }
            trPN.y = bv[1] - ang * (trPN.x - bv[0]);
        }
        //   c. Cheking.
        if (trPN.x < (uv[0] + _MIN_VALUE_1)) {
            trPN.x = uv[0];
            trPN.y = uv[1];
            isTrDone = true;
            wTrPNI = 0;
        }
        //   d. Integration.
        tmp = integrate_chanel_slant_left(
                tl, //   -  Index of current time layer.
                //
                trPC, wTrPCI, //   -  double *bv,
                trPN, wTrPNI, //   -  double *uv,
                //
                indCurSqOx, //   -  Indices where trPC and trPN are.
                //
                bv[0], indRB, //   -  double rb  =  Right boundary by Ox.
                //
                indCurSqOy, //   -  Index of current square by Oy axis.
                //
                ox,
                oy,
                density);
        result += tmp;
        //   e. Updating.
        if (isTrDone == false) {
            //   We will compute more. We need to redefine some values.
            wTrPCI = wTrPNI;
            trPC = trPN;
            if (wTrPNI == 1) {
                indCurSqOy[0] += 1;
                indCurSqOy[1] += 1;
            }
            if (wTrPNI == 2) {
                indCurSqOx[0] -= 1;
                indCurSqOx[1] -= 1;
            }
            if (indCurSqOx[0] >= 0) {
                distOx = trPC.x - ox[ indCurSqOx[0] ];
            }
            if (indCurSqOx[0] < 0) {
                distOx = fabs(trPC.x - hx * indCurSqOx[0]);
            }
            if (indCurSqOy[1] >= 0) {
                distOy = oy[ indCurSqOy[1] ] - trPC.y;
            }
            if (indCurSqOy[1] < 0) {
                distOy = fabs(hy * indCurSqOy[1] - trPC.y);
            }
        }
    } while (!isTrDone);
    return result;
}

double integrate_right_triangle_bottom_right(
        point_t &bv,
        point_t &uv,
        int tl,
        const double *ox,
        const double *oy,
        double *density) {

    if (fabs(bv[0] - uv[0]) < _MIN_VALUE) return fabs(bv[0] - uv[0]);
    double ang = (uv[1] - bv[1]) / (uv[0] - bv[0]); //   -  Angle of slant line. Should be greater zero.
    if (fabs(ang) < _MIN_VALUE) return fabs(ang);

    point_t trPC, trPN; //   -  Travel point current; Travel point next;
    int wTrPCI = 0, wTrPNI = 0; //   -  Where travel point current is? Where travel point next is?   


    int indCurSqOx[2], indCurSqOy[2]; //   -  Index of current square by Ox and Oy axes.
    int indLB[2]; //   -  Index of left boundary.
    double distOx, distOy; //   -  Distance to near Ox and Oy straight lines.
    bool isTrDone = false; //   -  Is travel done.
    double hx = ox[1] - ox[0];
    double hy = oy[1] - oy[0];
    double result = 0.; //   -  Value which we are computing.
    double tmp;

    trPC = bv;

    indCurSqOx[0] = (int) ((trPC.x + _MIN_VALUE_1) / hx); //   -  If trPC.x is in grid edge I want it will be between in the right side.

    if (trPC.x + _MIN_VALUE_1 <= 0)
        indCurSqOx[0] -= 1; //   -  The case when "trPC.x" is negative.

    indCurSqOx[1] = indCurSqOx[0] + 1; //   -  It's important only in rare case then trPC is in grid edge.
    indLB[0] = indCurSqOx[0];
    indLB[1] = indLB[0] + 1;
    indCurSqOy[0] = (int) ((trPC.y + _MIN_VALUE_1) / hy); //   -  If trPC.y is in grid edge I want it will be in the upper side.
    if ((trPC.y + _MIN_VALUE_1) <= 0) {
        indCurSqOy[0] -= 1; //   -  The case when "trPC.x" ia negative.
    }
    indCurSqOy[1] = indCurSqOy[0] + 1; //   -  It's important only in rare case then trPC is in grid edge.

    if (indCurSqOx[1] >= 0) {
        distOx = fabs(ox[ indCurSqOx[1] ] - trPC.x);
    }
    if (indCurSqOx[1] < 0) {
        distOx = fabs(hx * indCurSqOx[1] - trPC.x);
    }
    if (indCurSqOy[1] >= 0) {
        distOy = fabs(oy[ indCurSqOy[1] ] - trPC.y);
    }
    if (indCurSqOy[1] < 0) {
        distOy = fabs(hy * indCurSqOy[1] - trPC.y);
    }
    do {
        //   a. First case.
        if ((distOy / distOx) <= ang) {
            //   Across with straight line parallel Ox axis.
            wTrPNI = 1;
            if (indCurSqOy[1] >= 0) {
                trPN.y = oy[ indCurSqOy[1] ];
            }
            if (indCurSqOy[1] < 0) {
                trPN.y = hy * indCurSqOy[1];
            }
            trPN.x = bv[0] + (trPN.y - bv[1]) / ang;
        }
        //   b. Second case.
        if ((distOy / distOx) > ang) {
            //   Across with straight line parallel Oy axis.
            wTrPNI = 2;
            if (indCurSqOx[1] >= 0) {
                trPN.x = ox[ indCurSqOx[1] ];
            }
            if (indCurSqOx[1] < 0) {
                trPN.x = hx * indCurSqOx[1];
            }
            trPN.y = bv[1] + ang * (trPN.x - bv[0]);
        }
        //   c. Checking.
        if (trPN.x > (uv[0] - _MIN_VALUE_1)) {
            trPN = uv;
            isTrDone = true;
            wTrPNI = 0;
        }
        //   d. Integration.
        tmp = integrate_chanel_slant_right(
                tl,
                trPC, wTrPCI, //   -  double *bv,
                trPN, wTrPNI, //   -  double *uv,
                //
                indCurSqOx, //   -  Indices where trPC and trPN are.
                //
                bv.x, indLB, //   -  double lb  =  Left boundary by Ox.
                //
                indCurSqOy, //   -  Index of current square by Oy axis.
                //
                ox,
                oy,
                density);
        result += tmp;
        //   e. Updating.
        if (isTrDone == false) {
            //   We will compute more. We need to redefine some values.
            wTrPCI = wTrPNI;
            trPC = trPN;
            if (wTrPNI == 1) {
                indCurSqOy[0] += 1;
                indCurSqOy[1] += 1;
            }
            if (wTrPNI == 2) {
                indCurSqOx[0] += 1;
                indCurSqOx[1] += 1;
            }
            if (indCurSqOx[1] >= 0) {
                distOx = fabs(ox[ indCurSqOx[1] ] - trPC.x);
            }
            if (indCurSqOx[1] < 0) {
                distOx = fabs(hx * indCurSqOx[1] - trPC.x);
            }
            if (indCurSqOy[1] >= 0) {
                distOy = fabs(oy[ indCurSqOy[1] ] - trPC.y);
            }
            if (indCurSqOy[1] < 0) {
                distOy = fabs(hy * indCurSqOy[1] - trPC.y);
            }
        }
    } while (!isTrDone);
    return result;
}

double integrate_right_triangle_upper_left(
        point_t &bv,
        point_t &uv,
        int tl,
        const double *ox,
        const double *oy,
        double *density) {
    point_t trPC; //   -  Travel point current;
    int wTrPCI = 0; //   -  Where travel point current is?
    point_t trPN; //   -  Travel point next;
    int wTrPNI = 0; //   -  Where travel point next is?
    double ang; //   -  Angle of slant line. Should be greater zero.
    int indCurSqOx[2], indCurSqOy[2]; //   -  Index of current square by Ox and Oy axes.
    int indRB[2]; //   -  Index of right boundary.
    double distOx, distOy; //   -  Distance to near Ox and Oy straight lines.
    bool isTrDone = false; //   -  Is travel done.
    double hx = ox[1] - ox[0];
    double hy = oy[1] - oy[0];
    double result = 0.; //   -  Value which we are computing.
    double tmp;
    //   Initial data.
    trPC.x = bv[0];
    trPC.y = bv[1];
    if ((fabs(bv[0] - uv[0])) < _MIN_VALUE) return fabs(bv[0] - uv[0]);

    ang = (uv[1] - bv[1]) / (uv[0] - bv[0]);
    if (fabs(ang) < _MIN_VALUE) return fabs(ang);

    //   The follow equations are quite important.
    indCurSqOx[0] = (int) ((trPC.x + _MIN_VALUE_1) / hx); //   -  If trPC.x is in grid edge I want it will be in the right side.
    if ((trPC.x + _MIN_VALUE_1) <= 0) {
        indCurSqOx[0] -= 1; //   -  The case when "trPC.x" ia negative.
    }
    indCurSqOx[1] = indCurSqOx[0] + 1; //   -  It's important only in rare case then trPC is in grid edge.
    indCurSqOy[0] = (int) ((trPC.y + _MIN_VALUE_1) / hy); //   -  If trPC.y is in grid edge I want it will be in the upper square.
    if ((trPC.y + _MIN_VALUE_1) <= 0) {
        indCurSqOy[0] -= 1; //   -  The case when "trPC.x" ia negative.
    }
    indCurSqOy[1] = indCurSqOy[0] + 1;
    indRB[0] = (int) ((uv[0] - _MIN_VALUE_1) / hy); //   -  If uv[0] is in grid edge I want it will be in the left side.
    if ((uv[0] - _MIN_VALUE_1) <= 0) {
        indRB[0] -= 1; //   -  The case when "trPC.x" ia negative.
    }
    indRB[1] = indRB[0] + 1;
    if (indCurSqOx[1] >= 0) {
        distOx = ox[ indCurSqOx[1] ] - trPC.x;
    }
    if (indCurSqOx[1] < 0) {
        distOx = fabs(hx * indCurSqOx[1] - trPC.x);
    }
    if (indCurSqOy[1] >= 0) {
        distOy = oy[ indCurSqOy[1] ] - trPC.y;
    }
    if (indCurSqOy[1] < 0) {
        distOy = fabs(hy * indCurSqOy[1] - trPC.y);
    }
    do {
        //   a. First case.
        if ((distOy / distOx) <= ang) {
            //   Across with straight line parallel Ox axis.
            wTrPNI = 1;
            if (indCurSqOy[1] >= 0) {
                trPN.y = oy[ indCurSqOy[1] ];
            }
            if (indCurSqOy[1] < 0) {
                trPN.y = hy * indCurSqOy[1];
            }
            trPN.x = bv[0] + (trPN.y - bv[1]) / ang;
        }
        //   b. Second case.
        if ((distOy / distOx) > ang) {
            //   Across with straight line parallel Oy axis.
            wTrPNI = 2;
            if (indCurSqOx[1] >= 0) {
                trPN.x = ox[ indCurSqOx[1] ];
            }
            if (indCurSqOx[1] < 0) {
                trPN.x = hx * indCurSqOx[1];
            }
            trPN.y = bv[1] + ang * (trPN.x - bv[0]);
        }
        //   c. Cheking.
        if (trPN.x > (uv[0] - _MIN_VALUE_1)) {
            trPN = uv;
            isTrDone = true;
            wTrPNI = 0;
        }
        //   d. Integration.
        tmp = integrate_chanel_slant_left(
                tl,
                trPC, wTrPCI, //   -  double *bv,
                trPN, wTrPNI, //   -  double *uv,
                //
                indCurSqOx, //   -  Indices where trPC and trPN are.
                //
                uv.x, indRB, //   -  double rb  =  Right boundary by Ox.
                //
                indCurSqOy, //   -  Index of current square by Oy axis.
                //
                ox,
                oy,
                density);
        result += tmp;
        //   e. Updating.
        if (isTrDone == false) {
            //   We will compute more. We need to redefine some values.
            wTrPCI = wTrPNI;
            trPC = trPN;
            if (wTrPNI == 1) {
                indCurSqOy[0] += 1;
                indCurSqOy[1] += 1;
            }
            if (wTrPNI == 2) {
                indCurSqOx[0] += 1;
                indCurSqOx[1] += 1;
            }
            if (indCurSqOx[1] >= 0) {
                distOx = fabs(ox[ indCurSqOx[1] ] - trPC.x);
            }
            if (indCurSqOx[1] < 0) {
                distOx = fabs(hx * indCurSqOx[1] - trPC.x);
            }
            if (indCurSqOy[1] >= 0) {
                distOy = fabs(oy[ indCurSqOy[1] ] - trPC.y);
            }
            if (indCurSqOy[1] < 0) {
                distOy = fabs(hy * indCurSqOy[1] - trPC.y);
            }
        }
    } while (!isTrDone);
    return result;
}

double integrate_right_triangle_upper_right(
        point_t &bv,
        point_t &uv,
        int tl,
        const double *ox,
        const double *oy,
        double *density) {
    point_t trPC; //   -  Travel point current;
    int wTrPCI = 0; //   -  Where travel point current is?
    point_t trPN; //   -  Travel point next;
    int wTrPNI = 0; //   -  Where travel point next is?
    double ang; //   -  Angle of slant line. Should be greater zero.
    int indCurSqOx[2], indCurSqOy[2]; //   -  Index of current square by Ox and Oy axes.
    int indLB[2]; //   -  Index of left boundary.
    double distOx, distOy; //   -  Distance to near Ox and Oy straight lines.
    bool isTrDone = false; //   -  Is travel done.
    double hx = ox[1] - ox[0];
    double hy = oy[1] - oy[0];
    double result = 0.;
    double tmp;
    //   Initial data.
    trPC = bv;
    if ((fabs(bv[0] - uv[0])) < _MIN_VALUE) {
        //   This triangle has very small width. I guess further computation isn't correct.
        return fabs(bv[0] - uv[0]);
    }
    ang = (uv[1] - bv[1]) / (bv[0] - uv[0]);
    if (fabs(ang) < _MIN_VALUE) {
        //   This triangle has very small height. I guess further computation isn't correct.
        return fabs(ang);
    }
    indCurSqOx[0] = (int) ((trPC.x - _MIN_VALUE_1) / hx); //   -  If trPC.x is in grid edge I want it will be between in the left side.
    if ((trPC.x - _MIN_VALUE_1) <= 0) {
        indCurSqOx[0] -= 1; //   -  The case when "trPC.x" ia negative.
    }
    indCurSqOx[1] = indCurSqOx[0] + 1; //   -  It's important only in rare case then trPC is in grid edge.
    indLB[0] = (int) ((uv[0] + _MIN_VALUE_1) / hx);
    if ((uv[0] + _MIN_VALUE_1) <= 0) {
        indLB[0] -= 1; //   -  The case when "trPC.x" ia negative.
    }
    indLB[1] = indLB[0] + 1;
    indCurSqOy[0] = (int) ((trPC.y + _MIN_VALUE_1) / hy); //   -  If trPC.y is in grid edge I want it will be in the upper side.
    if ((trPC.y + _MIN_VALUE_1) <= 0) {
        indCurSqOy[0] -= 1; //   -  The case when "trPC.x" ia negative.
    }
    indCurSqOy[1] = indCurSqOy[0] + 1; //   -  It's important only in rare case then trPC is in grid edge.
    if (indCurSqOx[0] >= 0) {
        distOx = fabs(trPC.x - ox[ indCurSqOx[0] ]);
    }
    if (indCurSqOx[0] < 0) {
        distOx = fabs(trPC.x - hx * indCurSqOx[0]);
    }
    if (indCurSqOy[1] >= 0) {
        distOy = fabs(oy[ indCurSqOy[1] ] - trPC.y);
    }
    if (indCurSqOy[1] < 0) {
        distOy = fabs(hy * indCurSqOy[1] - trPC.y);
    }
    do {
        //   a. First case.
        if ((distOy / distOx) <= ang) {
            //   Across with straight line parallel Ox axis.
            wTrPNI = 1;
            if (indCurSqOy[1] >= 0) {
                trPN.y = oy[ indCurSqOy[1] ];
            }
            if (indCurSqOy[1] < 0) {
                trPN.y = hy * indCurSqOy[1];
            }
            trPN.x = bv[0] - (trPN.y - bv[1]) / ang;
        }
        //   b. Second case.
        if ((distOy / distOx) > ang) {
            //   Across with straight line parallel Oy axis.
            wTrPNI = 2;
            if (indCurSqOx[0] >= 0) {
                trPN.x = ox[ indCurSqOx[0] ];
            }
            if (indCurSqOx[0] < 0) {
                trPN.x = hx * indCurSqOx[0];
            }
            trPN.y = bv[1] - ang * (trPN.x - bv[0]);
        }
        //   c. Checking.
        if (trPN.x < (uv[0] + _MIN_VALUE_1)) {
            trPN = uv;
            isTrDone = true;
            wTrPNI = 0;
        }
        //   d. Integration.
        tmp = integrate_chanel_slant_right(tl, //   -  Index of current time layer.
                //
                trPC, wTrPCI, //   -  double *bv,
                trPN, wTrPNI, //   -  double *uv,
                //
                indCurSqOx, //   -  Indices where trPC and trPN are.
                //
                uv.x, indLB, //   -  double lb  =  Left boundary by Ox.
                //
                indCurSqOy, //   -  Index of current square by Oy axis.
                //
                ox,
                oy,
                density);
        result += tmp;
        //   e. Updating.
        if (isTrDone == false) {
            //   We will compute more. We need to redefine some values.
            wTrPCI = wTrPNI;
            trPC = trPN;
            if (wTrPNI == 1) {
                indCurSqOy[0] += 1;
                indCurSqOy[1] += 1;
            }
            if (wTrPNI == 2) {
                indCurSqOx[0] -= 1;
                indCurSqOx[1] -= 1;
            }
            if (indCurSqOx[0] >= 0) {
                distOx = fabs(trPC.x - ox[ indCurSqOx[0] ]);
            }
            if (indCurSqOx[0] < 0) {
                distOx = fabs(trPC.x - hx * indCurSqOx[0]);
            }
            if (indCurSqOy[1] >= 0) {
                distOy = fabs(oy[ indCurSqOy[1] ] - trPC.y);
            }
            if (indCurSqOy[1] < 0) {
                distOy = fabs(hy * indCurSqOy[1] - trPC.y);
            }
        }
    } while (!isTrDone);
    return result;
}

double integrate_bottom_triangle(int tl,
        point_t &l, //   -  Left vertex of Bottom triangle
        point_t &r, //   -  Right vertex of Bottom triangle
        point_t &m, //   -  Vertex of between L and R
        const double *ox,
        const double *oy,
        double *density) {
    double result = 0.;
    if (m.x == l.x) {
        result = integrate_right_triangle_bottom_right(m, r, tl, ox, oy, density);
    } else if (m.x == r.x) {
        result = integrate_right_triangle_bottom_left(m, l, tl, ox, oy, density);
    } else if (m.x < l.x) {
        result = integrate_right_triangle_bottom_right(m, r, tl, ox, oy, density);
        result -= integrate_right_triangle_bottom_right(m, l, tl, ox, oy, density);
    } else if (m.x > l.x && m.x < r.x) {
        result = integrate_right_triangle_bottom_left(m, l, tl, ox, oy, density);
        result += integrate_right_triangle_bottom_right(m, r, tl, ox, oy, density);
    } else if (m.x > r.x) {
        result = integrate_right_triangle_bottom_left(m, l, tl, ox, oy, density);
        result -= integrate_right_triangle_bottom_left(m, r, tl, ox, oy, density);
    }
    return result;
}

double integrate_upper_triangle(int tl,
        point_t &l, //   -  Left vertex of Upper triangle
        point_t &r, //   -  Right vertex of Upper triangle
        point_t &m, //   -  Vertex of between L and R
        const double *ox, const double *oy,
        double *density) {
    double result = 0.;
    if (m.x == l.x) {
        result = integrate_right_triangle_upper_right(r, m, tl, ox, oy, density);
    } else if (m.x == r.x) {
        result = integrate_right_triangle_upper_left(l, m, tl, ox, oy, density);
    } else if (m.x < l.x) {
        result = integrate_right_triangle_upper_right(r, m, tl, ox, oy, density);
        result -= integrate_right_triangle_upper_right(l, m, tl, ox, oy, density);
    } else if (m.x > l.x && m.x < r.x) {
        result = integrate_right_triangle_upper_left(l, m, tl, ox, oy, density);
        result += integrate_right_triangle_upper_right(r, m, tl, ox, oy, density);
    } else if (m.x > r.x) {
        result = integrate_right_triangle_upper_left(l, m, tl, ox, oy, density);
        result -= integrate_right_triangle_upper_left(r, m, tl, ox, oy, density);
    }
    return result;
}

double integrate_uniform_triangle_wall(int tl,
        const point_t &a,
        const point_t &b,
        const point_t &c,
        const double *ox,
        const double *oy,
        double *density) {
    return 0;
}

double integrate_uniform_triangle(int tl,
        point_t &x,
        point_t &y,
        point_t &z,
        const double *ox,
        const double *oy,
        double *density) {
    //   a * x  +  b * y  = c.
    double a = z.y - x.y;
    if (fabs(a) < _MIN_VALUE) return _MIN_VALUE;
    double b = x.x - z.x;
    double c = b * x.y + a * x.x;
    point_t ip((c - b * y.y) / a, y.y);

    //   Возможны 2 случая расположения точки перечеения относительно средней
    //   слева или справа.
    //   есди средняя точка справа от точки пересечения
    //   обменяем местами  X координаты, чтобы использовать один код для расчета
    if (y.x >= ip.x) {
        double tx = y.x;
        y.x = ip.x;
        ip.x = tx;
    }

    return integrate_upper_triangle(tl, y, ip, z, ox, oy, density)
            + integrate_bottom_triangle(tl, y, ip, x, ox, oy, density);
}

inline double func_u(double x, double y) {
    return B * y * (1. - y) * (_PI / 2. + atan(-x));
}

inline double func_v(double t, double x, double y) {
    return atan((x - LB) * (x - RB) * (1. + t) / 10. * (y - UB) * (y - BB));
}

double func_f(
        double tl_on_tau,
        double x,
        double y) {
    double arg_v = (x - LB) * (x - RB) * (1. + tl_on_tau) / 10. * (y - UB) * (y - BB);
    double rho = analytical_solution(tl_on_tau, x, y);
    double dRhoDT = x * y * cos(tl_on_tau * x * y);
    double dRhoDX = tl_on_tau * y * cos(tl_on_tau * x * y);
    double dRhoDY = tl_on_tau * x * cos(tl_on_tau * x * y);
    double u = func_u(x, y);
    double v = func_v(tl_on_tau, x, y);
    double duDX = -B * y * (1. - y) / (1. + x * x);
    double dvDY = (x - LB) * (x - RB) * (1. + tl_on_tau) / 10. * (y - BB + y - UB);
    dvDY /= (1. + arg_v * arg_v);
    double res = dRhoDT + rho * duDX + u * dRhoDX + rho * dvDY + v * dRhoDY;

    //  printf("x = %f\n", x);
    //  printf("y = %f\n", y);
    //  printf("arg_v = %f\n", arg_v);
    //  printf("rho = %f\n", rho);
    //  printf("dRhoDT = %f\n", dRhoDT);
    //  printf("dRhoDX = %f\n", dRhoDX);
    //  printf("dRhoDY = %f\n", dRhoDY);
    //  printf("u = %f\n", u);
    //  printf("duDX = %f\n", duDX);
    //  printf("v = %f\n", v);
    //  printf("dvDY = %f\n", dvDY);
    //  printf("res = %f\n", res);

    return res;
}

quad_type get_coordinates_on_prev_layer(int cur_tl, int ix, int iy,
        const double *ox,
        const double *oy,
        point_t &alpha, point_t &beta, point_t &gamma, point_t &theta) {
    //   1. First of all let's compute coordinates of square vertexes.
    //  OX:
    if (ix == 0) {
        alpha.x = ox[ ix ];
        beta.x = (ox[ix] + ox[ix + 1]) / 2.;
        gamma.x = (ox[ix] + ox[ix + 1]) / 2.;
        theta.x = ox[ ix ];
    } else if (ix == OX_LEN) {
        alpha.x = (ox[ix - 1] + ox[ix]) / 2.;
        beta.x = ox[ ix ];
        gamma.x = ox[ ix ];
        theta.x = (ox[ix - 1] + ox[ix]) / 2.;
    } else {
        alpha.x = (ox[ix - 1] + ox[ix]) / 2.;
        beta.x = (ox[ix + 1] + ox[ix]) / 2.;
        gamma.x = (ox[ix + 1] + ox[ix]) / 2.;
        theta.x = (ox[ix - 1] + ox[ix]) / 2.;
    }

    //  OY:
    if (iy == 0) {
        alpha.y = oy[ iy ];
        beta.y = oy[ iy ];
        gamma.y = (oy[iy] + oy[ iy + 1]) / 2.;
        theta.y = (oy[iy] + oy[ iy + 1]) / 2.;
    } else if (iy == OY_LEN) {
        alpha.y = (oy[iy] + oy[ iy - 1]) / 2.;
        beta.y = (oy[iy] + oy[ iy - 1]) / 2.;
        gamma.y = oy[ iy ];
        theta.y = oy[ iy ];
    } else {
        alpha.y = (oy[iy] + oy[ iy - 1]) / 2.;
        beta.y = (oy[iy] + oy[ iy - 1]) / 2.;
        gamma.y = (oy[iy] + oy[ iy + 1]) / 2.;
        theta.y = (oy[iy] + oy[ iy + 1]) / 2.;
    }

    double u, v;

    // Now let's compute new coordinates on the previous time level of alpha, beta, gamma, theta points.
    u = func_u(alpha.x, alpha.y);
    v = func_v(TAU*cur_tl, alpha.x, alpha.y);
    alpha.x -= TAU * u;
    alpha.y -= TAU * v;

    u = func_u(beta.x, beta.y);
    v = func_v(TAU*cur_tl, beta.x, beta.y);
    beta.x -= TAU * u;
    beta.y -= TAU * v;

    u = func_u(gamma.x, gamma.y);
    v = func_v(TAU*cur_tl, gamma.x, gamma.y);
    gamma.x -= TAU * u;
    gamma.y -= TAU * v;

    u = func_u(theta.x, theta.y);
    v = func_v(TAU*cur_tl, theta.x, theta.y);
    theta.x -= TAU * u;
    theta.y -= TAU * v;

    point_t intersection = get_intersection_point(alpha, beta, gamma, theta);
    if ((beta.y - intersection.y) * (theta.y - intersection.y) > 0.) return pseudo; // ??
    if ((alpha.x - intersection.x) * (gamma.x - intersection.x) > 0.) return pseudo; // ??
    double product = get_vector_product(alpha, beta, theta); // ?
    if (product < 0.) return pseudo;

    // значит что точка улетела за левую границу
    if (theta.x < 0 ||
            theta.y < 0 ||
            beta.x < 0 ||
            beta.y < 0 ||
            gamma.x < 0 ||
            gamma.y < 0 ||
            alpha.x < 0 ||
            alpha.y < 0) {
        return normal;
        //return wall;
    }
    return normal;
}


// Type of quadrangle: 0 - pseudo; 1 - convex; 2 - concave;

quad_type get_quadrangle_type(int tl, int ix, int iy,
        const double *ox,
        const double *oy,
        point_t &a, //   -  First vertex of first triangle.
        point_t &b, //   -  Second vertex of first triangle.
        point_t &c, //   -  Third vertex of first triangle.
        point_t &k, //   -  First vertex of second triangle.
        point_t &m, //   -  Second vertex of second triangle.
        point_t &n) //   -  Third vertex of second triangle.
{
    point_t alpha, beta, gamma, theta; // coordinates on previous time layer
    quad_type type = get_coordinates_on_prev_layer(tl, ix, iy, ox, oy, alpha, beta, gamma, theta);
    a = alpha;
    b = beta;
    c = gamma;
    k = alpha;
    m = theta;
    n = gamma;
    return type;
}

double integrate(double tl, int ix, int iy,
        const double *ox,
        const double *oy,
        double *density) {
    point_t a1, b1, c1, a2, b2, c2;
    quad_type type = get_quadrangle_type(tl, ix, iy, ox, oy, a1, b1, c1, a2, b2, c2);
    if (type != normal && type != wall) {
        return -1.;
    }

    // чтобы правилно отработала процедура интегрирования
    // точки должны идти в порядке возрастания y координаты
    sort_by_y(a1, b1, c1);
    sort_by_y(a2, b2, c2);

    // check the type of triangle to select appropriate computation method
    double result = 0.;
    switch (type) {
        case wall:
            result += integrate_uniform_triangle_wall(tl, a1, b1, c1, ox, oy,
                    density);
            result += integrate_uniform_triangle_wall(tl, a2, b2, c2, ox, oy,
                    density);
            return result;
        case normal:
            result += integrate_uniform_triangle(tl, a1, b1, c1, ox, oy, density);
            result += integrate_uniform_triangle(tl, a2, b2, c2, ox, oy, density);
            return result;
        case concave:
        case convex:
        case pseudo:
            return 0.;
    }
}

double get_norm_of_error(double* density, int x_length, int y_length, double* ox,
        double* oy,
        double ts_count_mul_steps) {
    double result = 0.;
    for (int k = 1; k < y_length; ++k) {
        for (int j = 1; j < x_length; ++j) {
            result += fabs(analytical_solution(ts_count_mul_steps, ox[j], oy[k])
                    - density[ (x_length + 1) * k + j ]);
        }
    }
    double hx = ox[1] - ox[0];
    double hy = oy[1] - oy[0];
    return hx * hy * result;
}

double solve(const double *ox, const double *oy, double *density) {
    double *prev_density = new double [ XY_LEN ];
    for (int iy = 0; iy < OY_LEN + 1; iy++) {
        for (int ix = 0; ix < OX_LEN + 1; ix++) {
            prev_density[ (OX_LEN + 1) * iy + ix ] = analytical_solution(0., ox[ix], oy[iy]);
        }
    }

    for (int it = 1; it < TIME_STEP_CNT + 1; it++) {
        for (int i = 0; i <= OX_LEN; i++) {
            density[ i ] = init_side(ox[ i ], BB, TAU * it, bottom);
            density[ (OX_LEN + 1) * OY_LEN + i ] = init_side(ox[ i ], UB, TAU * it, up);
        }

        for (int i = 0; i <= OY_LEN; i++) {
            density[ (OX_LEN + 1) * i ] = init_side(LB, oy[ i ], TAU * it, left);
            density[ (OX_LEN + 1) * i + OX_LEN ] = init_side(RB, oy[ i ], TAU * it, right);
        }

        for (int iy = 1; iy < OY_LEN; iy++) {
            for (int ix = 1; ix < OX_LEN; ix++) {
                int index = (OX_LEN + 1) * iy + ix;

                double value = integrate(it, ix, iy, ox, oy, prev_density);

                double h = (ox[ix + 1] - ox[ix - 1]) / 2.;
                value /= h;
                h = (oy[iy + 1] - oy[iy - 1]) / 2.;
                value /= h;

                double rp = func_f(TAU*it, ox[ix], oy[iy]);
                density[ index ] = value + TAU * rp;
            }
        }
        memcpy(prev_density, density, (OX_LEN + 1) * (OY_LEN + 1) * sizeof (double));
    }

    delete[] prev_density;
    return 0;
}

double *compute_density(double b,
        double lb,
        double rb,
        double bb,
        double ub,
        double tau,
        int time_step_count,
        int ox_length,
        int oy_length,
        double* norm) {
    B = b;
    UB = ub;
    BB = bb;
    LB = lb;
    RB = rb;
    TAU = tau;
    OX_LEN = ox_length;
    OY_LEN = oy_length;
    TIME_STEP_CNT = time_step_count;
    XY_LEN = (ox_length + 1)*(oy_length + 1);

    //  point_t *t = new point_t(19,0);
    //  printf("%f \n", (*t)[0]);

    double *density = new double [ XY_LEN ];
    double *ox = new double [ OX_LEN + 1 ];
    double *oy = new double [ OY_LEN + 1 ];

    for (int i = 0; i <= OX_LEN; i++) {
        ox[i] = lb + i * (rb - lb) / OX_LEN;
    }

    for (int i = 0; i <= OY_LEN; i++) {
        oy[i] = bb + i * (ub - bb) / OY_LEN;
    }

    print_params(B, LB, RB, BB, UB, TAU, time_step_count, OX_LEN, OY_LEN);

    solve(ox, oy, density);

    *norm = get_norm_of_error(density, OX_LEN, OY_LEN, ox, oy,
            time_step_count * TAU);
    //  printf("Norm L1 = %f\n", *norm);
    printf("%d x %d wall count = %d\n", ox_length + 1, oy_length + 1, TMP_WALL_CNT);

    delete[] ox;
    delete[] oy;
    return density;
}
