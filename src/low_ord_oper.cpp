#include "common.h"
#include "point.h"

static const double C_pi = 3.14159265358979323846264338327;
static int wall_counter = 0;
static double PAR_B;

double _itemOfInteg_1SpecType(double Py,
                              double Qy,
                              double Gx,
                              double Hx,
                              double a,
                              double b)
{
    double integ = (Hx - a) * (Hx - a) - (Gx - a) * (Gx - a);
    integ *= (Qy - b) * (Qy - b) - (Py - b) * (Py - b);
    return integ / 4;
}

double analytical_solution(double t, double x, double y)
{
    return 1.1 + sin(t * x * y);
}

double _itemOfInteg_2SpecType(double Py,
                              double Qy,
                              double alpha,
                              double a,
                              double b,
                              double betta)
{
    double tmp, integ;
    tmp = (Qy - alpha) * (a * Qy + b - betta) * (a * Qy + b - betta) * (a * Qy + b - betta);
    tmp = tmp - (Py - alpha) * (a * Py + b - betta) * (a * Py + b - betta) * (a * Py + b - betta);
    integ = tmp / (3 * a);
    tmp = (a * Qy + b - betta) * (a * Qy + b - betta) * (a * Qy + b - betta) * (a * Qy + b - betta);
    tmp = tmp - (a * Py + b - betta) * (a * Py + b - betta) * (a * Py + b - betta) * (a * Py + b - betta);
    return integ - tmp / (12 * a * a);
}

double integUnderLeftTr_OneCell(
                                double Py,
                                double Qy,
                                //
                                double a_SL,
                                double b_SL,
                                double Hx,
                                //
                                double tau,
                                int iCurrTL, //   -  Index of current time layer.
                                //
                                int *indCurSqOx, //   -  Index of current square by Ox axis.
                                int *indCurSqOy, //   -  Index of current square by Oy axis.
                                //
                                const double *masOX, //   -  Massive of OX steps. Dimension = numOfOXSt +1.
                                int numOfOXSt, //   -  Number of OX steps.
                                //
                                const double *masOY, //   -  Massive of OY steps. Dimension = numOfOYSt +1.
                                int numOfOYSt, //   -  Number of OY steps.
                                //
                                double *rhoInPrevTL_asV)
{
    double hx = masOX[1] - masOX[0];
    double hy = masOY[1] - masOY[0];
    double integ = 0;
    double buf_D, bufInteg_D;
    double rho[2][2];
    double t = tau * (iCurrTL - 1.);
    double x, y;
    if (indCurSqOx[0] >= 0 && indCurSqOx[1] <= numOfOXSt)
    {
        if (indCurSqOy[0] >= 0 && indCurSqOy[1] <= numOfOYSt)
        {
            rho[0][0] = rhoInPrevTL_asV[ (numOfOXSt + 1) * indCurSqOy[0] + indCurSqOx[0] ];
            rho[0][1] = rhoInPrevTL_asV[ (numOfOXSt + 1) * indCurSqOy[1] + indCurSqOx[0] ];
            rho[1][0] = rhoInPrevTL_asV[ (numOfOXSt + 1) * indCurSqOy[0] + indCurSqOx[1] ];
            rho[1][1] = rhoInPrevTL_asV[ (numOfOXSt + 1) * indCurSqOy[1] + indCurSqOx[1] ];
        }
    }

    // TODO: убрать потому что это неверно (надо расчитывать граничные условия)
    // норма должна уменьшиться
    if (indCurSqOx[0] < 0 || indCurSqOx[1] > numOfOXSt || indCurSqOy[0] < 0 || indCurSqOy[1] > numOfOYSt)
    {
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

        wall_counter++;
    }

    //   1.
    buf_D = (Qy - masOY[ indCurSqOy[1] ]) * (Qy - masOY[ indCurSqOy[1] ]) - (Py - masOY[ indCurSqOy[1] ]) * (Py - masOY[ indCurSqOy[1] ]);
    if ((indCurSqOx[1] >= 0) && (indCurSqOy[1] >= 0))
    {
        buf_D = buf_D * (Hx - masOX[ indCurSqOx[1] ]) * (Hx - masOX[ indCurSqOx[1] ]) / 4.;
        bufInteg_D = _itemOfInteg_2SpecType(Py, Qy, masOY[ indCurSqOy[1] ], a_SL, b_SL, masOX[ indCurSqOx[1] ]);
    }
    else
    {
        buf_D = buf_D * (Hx - hx * indCurSqOx[1]) * (Hx - hx * indCurSqOx[1]) / 4.;
        bufInteg_D = _itemOfInteg_2SpecType(Py, Qy, hy * indCurSqOy[1], a_SL, b_SL, hx * indCurSqOx[1]);
    }
    buf_D = buf_D - bufInteg_D / 2.;
    integ = buf_D * rho[0][0] / hx / hy;

    //   2.
    buf_D = (Qy - masOY[ indCurSqOy[1] ]) * (Qy - masOY[ indCurSqOy[1] ]) - (Py - masOY[ indCurSqOy[1] ]) * (Py - masOY[ indCurSqOy[1] ]);
    if ((indCurSqOx[0] >= 0) && (indCurSqOy[1] >= 0))
    {
        buf_D = -1. * buf_D * (Hx - masOX[ indCurSqOx[0] ]) * (Hx - masOX[ indCurSqOx[0] ]) / 4.;
        bufInteg_D = _itemOfInteg_2SpecType(Py, Qy, masOY[ indCurSqOy[1] ], a_SL, b_SL, masOX[ indCurSqOx[0] ]);
    }
    else
    {
        buf_D = -1. * buf_D * (Hx - hx * indCurSqOx[0]) * (Hx - hx * indCurSqOx[0]) / 4.;
        bufInteg_D = _itemOfInteg_2SpecType(Py, Qy, hy * indCurSqOy[1], a_SL, b_SL, hx * indCurSqOx[0]);
    }
    buf_D = buf_D + bufInteg_D / 2.;
    integ = integ + buf_D * rho[1][0] / hx / hy;

    //   3.
    buf_D = (Qy - masOY[ indCurSqOy[0] ]) * (Qy - masOY[ indCurSqOy[0] ]) - (Py - masOY[ indCurSqOy[0] ]) * (Py - masOY[ indCurSqOy[0] ]);
    if ((indCurSqOx[1] >= 0) && (indCurSqOy[0] >= 0))
    {
        buf_D = -1. * buf_D * (Hx - masOX[ indCurSqOx[1] ]) * (Hx - masOX[ indCurSqOx[1] ]) / 4.;
        bufInteg_D = _itemOfInteg_2SpecType(Py, Qy, masOY[ indCurSqOy[0] ], a_SL, b_SL, masOX[ indCurSqOx[1] ]);
    }
    else
    {
        buf_D = -1. * buf_D * (Hx - hx * indCurSqOx[1]) * (Hx - hx * indCurSqOx[1]) / 4.;
        bufInteg_D = _itemOfInteg_2SpecType(Py, Qy, hy * indCurSqOy[0], a_SL, b_SL, hx * indCurSqOx[1]);
    }
    buf_D = buf_D + bufInteg_D / 2.;
    integ = integ + buf_D * rho[0][1] / hx / hy;

    //   4.
    buf_D = (Qy - masOY[ indCurSqOy[0] ]) * (Qy - masOY[ indCurSqOy[0] ]) - (Py - masOY[ indCurSqOy[0] ]) * (Py - masOY[ indCurSqOy[0] ]);
    if ((indCurSqOx[0] >= 0) && (indCurSqOy[0] >= 0))
    {
        buf_D = buf_D * (Hx - masOX[ indCurSqOx[0] ]) * (Hx - masOX[ indCurSqOx[0] ]) / 4.;
        bufInteg_D = _itemOfInteg_2SpecType(Py, Qy, masOY[ indCurSqOy[0] ], a_SL, b_SL, masOX[ indCurSqOx[0] ]);
    }
    else
    {
        buf_D = buf_D * (Hx - hx * indCurSqOx[0]) * (Hx - hx * indCurSqOx[0]) / 4.;
        bufInteg_D = _itemOfInteg_2SpecType(Py, Qy, hy * indCurSqOy[0], a_SL, b_SL, hx * indCurSqOx[0]);
    }
    buf_D = buf_D - bufInteg_D / 2.;
    integ += buf_D * rho[1][1] / hx / hy;

    return integ;
}

double integUnderRightTr_OneCell(double Py,
                                 double Qy,
                                 double a_SL,
                                 double b_SL,
                                 double Gx,
                                 double tau,
                                 int tl,
                                 int *indCurSqOx, //   -  Index of current square by Ox axis.
                                 int *indCurSqOy, //   -  Index of current square by Oy axis.
                                 const double *ox, //   -  Massive of OX steps. Dimension = numOfOXSt +1.
                                 int numOfOXSt,
                                 const double *oy, //   -  Massive of OY steps. Dimension = numOfOYSt +1.
                                 int numOfOYSt,
                                 double *prev_density)
{
    return -1. * integUnderLeftTr_OneCell(
                                          Py, Qy,
                                          //
                                          a_SL, b_SL,
                                          Gx, //   -  double Hx,
                                          //
                                          tau, tl, //   -  Index of current time layer.
                                          //
                                          indCurSqOx, //   -  Index of current square by Ox axis.
                                          indCurSqOy, //   -  Index of current square by Oy axis.
                                          //
                                          ox, numOfOXSt, //   -  Massive of OX steps. Dimension = numOfOXSt +1. Number of OX steps.
                                          //
                                          oy, numOfOYSt, //   -  Massive of OY steps. Dimension = numOfOYSt +1. Number of OY steps.
                                          //
                                          prev_density);
}

double integUnderRectAng_OneCell(double Py,
                                 double Qy,
                                 double Gx,
                                 double Hx,
                                 double tau,
                                 int iCurrTL,
                                 int *indCurSqOx, //   -  Index of current square by Ox axis.
                                 int *indCurSqOy,
                                 const double *masOX, //   -  Massive of OX steps. Dimension = numOfOXSt +1.
                                 int numOfOXSt,
                                 const double *masOY,
                                 double *rhoInPrevTL_asV)
{
    double hx = masOX[1] - masOX[0];
    double hy = masOY[1] - masOY[0];
    double integ = 0;
    double buf_D;
    double rho[2][2];
    double t = tau * (iCurrTL - 1.);
    double x, y;
    if (indCurSqOx[0] >= 0 && indCurSqOy[0] >= 0)
    {
        rho[0][0] = rhoInPrevTL_asV[ (numOfOXSt + 1) * indCurSqOy[0] + indCurSqOx[0] ];
        rho[0][1] = rhoInPrevTL_asV[ (numOfOXSt + 1) * indCurSqOy[1] + indCurSqOx[0] ];
        rho[1][0] = rhoInPrevTL_asV[ (numOfOXSt + 1) * indCurSqOy[0] + indCurSqOx[1] ];
        rho[1][1] = rhoInPrevTL_asV[ (numOfOXSt + 1) * indCurSqOy[1] + indCurSqOx[1] ];
    }
    else
    {
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

        wall_counter++;

    }

    if (indCurSqOx[1] >= 0 && indCurSqOy[1] >= 0)
    {
        buf_D = _itemOfInteg_1SpecType(Py, Qy, Gx, Hx, masOX[ indCurSqOx[1] ], masOY[ indCurSqOy[1] ]);
    }
    else
    {
        buf_D = _itemOfInteg_1SpecType(Py, Qy, Gx, Hx, hx * indCurSqOx[1], hy * indCurSqOy[1]);
    }
    buf_D = buf_D / hx / hy;
    integ = buf_D * rho[0][0]; //   rhoInPrevTL[ indCurSqOx[0] ][ indCurSqOy[0] ];
    if (indCurSqOx[0] >= 0 && indCurSqOy[1] >= 0)
    {
        buf_D = _itemOfInteg_1SpecType(Py, Qy, Gx, Hx, masOX[ indCurSqOx[0] ], masOY[ indCurSqOy[1] ]);
    }
    else
    {
        buf_D = _itemOfInteg_1SpecType(Py, Qy, Gx, Hx, hx * indCurSqOx[0], hy * indCurSqOy[1]);
    }
    buf_D = buf_D / hx / hy;
    integ = integ - buf_D * rho[1][0]; //   rhoInPrevTL[ indCurSqOx[1] ][ indCurSqOy[0] ];
    if (indCurSqOx[1] >= 0 && indCurSqOy[0] >= 0)
    {
        buf_D = _itemOfInteg_1SpecType(Py, Qy, Gx, Hx, masOX[ indCurSqOx[1] ], masOY[ indCurSqOy[0] ]);
    }
    else
    {
        buf_D = _itemOfInteg_1SpecType(Py, Qy, Gx, Hx, hx * indCurSqOx[1], hy * indCurSqOy[0]);
    }
    buf_D = buf_D / hx / hy;
    integ = integ - buf_D * rho[0][1]; //   rhoInPrevTL[ indCurSqOx[0] ][ indCurSqOy[1] ];
    if (indCurSqOx[0] >= 0 && indCurSqOy[0] >= 0)
    {
        buf_D = _itemOfInteg_1SpecType(Py, Qy, Gx, Hx, masOX[ indCurSqOx[0] ], masOY[ indCurSqOy[0] ]);
    }
    else
    {
        buf_D = _itemOfInteg_1SpecType(Py, Qy, Gx, Hx, hx * indCurSqOx[0], hy * indCurSqOy[0]);
    }
    buf_D = buf_D / hx / hy;
    return integ + buf_D * rho[1][1]; //   rhoInPrevTL[ indCurSqOx[1] ][ indCurSqOy[1] ];
}

double integOfChan_SLRightSd(double tau,
                             int iCurrTL,
                             double *bv, int wTrPCI, //   -  Where travel point current (botton vertex) is.
                             double *uv, int wTrPNI, //   -  Where travel point next (upper vertex) is.
                             //
                             int *indCurSqOx, //   -  Index by OX axis where bv and uv are.
                             //
                             double lb, int *indLB, //   -  Left boundary by Ox. Index by OX axis where lb is.
                             //
                             int *indCurSqOy, //   -  Index of current square by Oy axis.
                             //
                             const double *masOX, //   -  Massive of OX steps. Dimension = numOfOXSt +1.
                             int numOfOXSt, //   -  Number of OX steps.
                             //
                             const double *masOY, //   -  Massive of OY steps. Dimension = numOfOYSt +1.
                             int numOfOYSt, //   -  Number of OY steps.
                             //
                             double *rhoInPrevTL_asV)
{
    double mv[2], rv[2]; //   -  Middle and right vertices.
    int wMvI; //   -  Where middle vertex is.
    int indCurSqOxToCh[2]; //   -  Indices of current square by Ox axis to be changed. Under which we want to integrate.
    double h = masOX[1] - masOX[0];
    double a_SL, b_SL; //   -  Coefficients of slant line: x = a_SL *y  +  b_SL.
    double Gx, Hx; //   -  Left boundary for each integration.
    double integ = 0.;
    double buf_D;

    //   Let's compute helpful values.

    if (uv[0] <= bv[0])
    {
        mv[0] = uv[0];
        mv[1] = uv[1];
        wMvI = wTrPNI;
        rv[0] = bv[0];
        rv[1] = bv[1];
    }

    if (uv[0] > bv[0])
    {
        mv[0] = bv[0];
        mv[1] = bv[1];
        wMvI = wTrPCI;
        rv[0] = uv[0];
        rv[1] = uv[1];
    }

    if ((fabs(uv[1] - bv[1])) <= 1.e-12)
    {
        //   Computation is impossible. Too smale values. Let's return some approximate value.
        //   buf_D  =  (uv[1] - bv[1])  *  ((uv[0] + bv[0]) /2.  -  lb) * rhoInPrevTL[ indCurSqOx[0] ][ indCurSqOy[0] ];
        return fabs(uv[1] - bv[1]); //   fabs(uv[1] - bv[1]);
    }


    //   First step: from "lb" to "masOX[ indCurSqOx[0] ]" by iteration.
    //   integ  += fabs( mv[0] - lb) * fabs(uv[1] - bv[1]);

    indCurSqOxToCh[0] = indLB[0];
    indCurSqOxToCh[1] = indCurSqOxToCh[0] + 1;

    for (int j = indLB[0]; j < indCurSqOx[0]; j++)
    {
        //   If this is first cell we should integrate under rectangle only.
        if (indCurSqOxToCh[0] >= 0)
        {
            Gx = masOX[ indCurSqOxToCh[0] ];
            Hx = masOX[ indCurSqOxToCh[1] ];
        }


        if (indCurSqOxToCh[0] < 0)
        {
            Gx = h * indCurSqOxToCh[0];
            Hx = h * indCurSqOxToCh[1];
        }

        if (j == indLB[0])
        {
            Gx = lb;
        }

        buf_D = integUnderRectAng_OneCell(
                                          bv[1], //   -  double Py,
                                          uv[1], //   -  double Qy,
                                          //
                                          Gx, //   -  double Gx,
                                          Hx, //   -  double Hx,
                                          //
                                          tau, iCurrTL, //   -  Index of current time layer.
                                          //
                                          indCurSqOxToCh, //   -  Index of current square by Ox axis.
                                          indCurSqOy, //   -  Index of current square by Oy axis.
                                          //
                                          masOX, //   -  Massive of OX steps. Dimension = numOfOXSt +1.
                                          numOfOXSt, //   -  Number of OX steps.
                                          //
                                          masOY, //   -  Massive of OY steps. Dimension = numOfOYSt +1.

                                          rhoInPrevTL_asV);

        integ += buf_D;

        indCurSqOxToCh[0] += 1;
        indCurSqOxToCh[1] = indCurSqOxToCh[0] + 1;
    }

    //   Integration. Second step: under [ indCurSqOx[0]; indCurSqOx[1] ] square.

    //   A. Under rectangle.
    if (wMvI == 1)
    {
        if (indCurSqOx[0] == indLB[0])
        {
            Gx = lb;
        }

        if (indCurSqOx[0] > indLB[0])
        {
            if (indCurSqOx[0] >= 0)
            {
                Gx = masOX[ indCurSqOx[0] ];
            }

            if (indCurSqOx[0] < 0)
            {
                Gx = h * indCurSqOx[0];
            }
        }

        buf_D = integUnderRectAng_OneCell(bv[1], //   -  double Py,
                                          uv[1], //   -  double Qy,
                                          //
                                          Gx, //   -  double Gx,
                                          mv[0], //   -  double Hx,
                                          //
                                          tau, iCurrTL, //   -  Index of current time layer.
                                          //
                                          indCurSqOx, //   -  Index of current square by Ox axis.
                                          indCurSqOy, //   -  Index of current square by Oy axis.
                                          //
                                          masOX, //   -  Massive of OX steps. Dimension = numOfOXSt +1.
                                          numOfOXSt, //   -  Number of OX steps.
                                          //
                                          masOY, //   -  Massive of OY steps. Dimension = numOfOYSt +1.

                                          rhoInPrevTL_asV);

        integ += buf_D;

    }

    //   B. Under triangle.

    if (fabs(uv[1] - bv[1]) > 1.e-12)
    {
        //   integ += fabs(uv[1] - bv[1]) * (rv[0] - mv[0]) /2.;
        //   Coefficients of slant line: x = a_SL *y  +  b_SL.
        a_SL = (uv[0] - bv[0]) / (uv[1] - bv[1]);
        b_SL = bv[0] - a_SL * bv[1];


        //   Integration under one cell triangle.

        if (fabs(a_SL) > 1.e-12)
        {
            buf_D = integUnderRightTr_OneCell(
                                              bv[1], //   -  double Py,
                                              uv[1], //   -  double Qy,
                                              //
                                              a_SL,
                                              b_SL,
                                              mv[0], //   -  double Gx,
                                              //
                                              tau, iCurrTL, //   -  Index of current time layer.
                                              //
                                              indCurSqOx, //   -  Index of current square by Ox axis.
                                              indCurSqOy, //   -  Index of current square by Oy axis.
                                              //
                                              masOX, //   -  Massive of OX steps. Dimension = numOfOXSt +1.
                                              numOfOXSt, //   -  Number of OX steps.
                                              //
                                              masOY, //   -  Massive of OY steps. Dimension = numOfOYSt +1.
                                              numOfOYSt, //   -  Number of OY steps.
                                              //
                                              rhoInPrevTL_asV);

            integ += buf_D;
        }
    }

    return integ;
}

double init_bound(double x, double y, double t, bound_side side)
{
    switch (side)
    {
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

double integOfChan_SLLeftSd(
                            double tau,
                            int iCurrTL, //   -  Index of current time layer.
                            //
                            double *bv, int wTrPCI, //   -  Where travel point current (botton vertex) is.
                            double *uv, int wTrPNI, //   -  Where travel point next (upper vertex) is.
                            //
                            int *indCurSqOx, //   -  Index by OX axis where bv and uv are.
                            //
                            double rb, int *indRB, //   -  Right boundary by Ox. Index by OX axis where rb is.
                            //
                            int *indCurSqOy, //   -  Index of current square by Oy axis.
                            //
                            const double *masOX, //   -  Massive of OX steps. Dimension = numOfOXSt +1.
                            int numOfOXSt, //   -  Number of OX steps.
                            //
                            const double *masOY, //   -  Massive of OY steps. Dimension = numOfOYSt +1.
                            int numOfOYSt, //   -  Number of OY steps.
                            //
                            double *rhoInPrevTL_asV)
{
    double lv[2], mv[2]; //   -  Left and middle vertices.
    int wMvI; //   -  Where middle vertex is.
    int indCurSqOxToCh[2]; //   -  Indices of current square by Ox axis to be changed. Under which we want to integrate.
    double h = masOX[1] - masOX[0];
    double a_SL, b_SL; //   -  Coefficients of slant line: x = a_SL *y  +  b_SL.
    double Gx, Hx; //   -  Left and right boundary for each integration.
    double integ = 0.;
    double buf_D;
    int j;

    //   Let's compute helpful values.

    if (uv[0] <= bv[0])
    {
        lv[0] = uv[0];
        lv[1] = uv[1];
        mv[0] = bv[0];
        mv[1] = bv[1];
        wMvI = wTrPCI;
    }

    if (uv[0] > bv[0])
    {
        lv[0] = bv[0];
        lv[1] = bv[1];
        mv[0] = uv[0];
        mv[1] = uv[1];
        wMvI = wTrPNI;
    }

    if ((fabs(uv[1] - bv[1])) <= 1.e-12)
    {
        //   Computation is impossible. Too smale values. Let's return some approximate value.
        //   buf_D  =  (uv[1] - bv[1])  *  (rb  - (uv[0] + bv[0]) /2.) * rhoInPrevTL[ indCurSqOx[0] ][ indCurSqOy[0] ];
        return fabs(uv[1] - bv[1]); //   fabs(uv[1] - bv[1]);
    }

    //   Integration. First step: under [ indCurSqOx[0]; indCurSqOx[1] ] square.

    //   A. Under triangle.

    if (fabs(uv[1] - bv[1]) > 1.e-12)
    {
        //   Coefficients of slant line: x = a_SL *y  +  b_SL.
        a_SL = (uv[0] - bv[0]) / (uv[1] - bv[1]);
        b_SL = bv[0] - a_SL * bv[1];

        //   Integration under one cell triangle.
        if (fabs(a_SL) > 1.e-12)
        {
            buf_D = integUnderLeftTr_OneCell(
                                             bv[1], //   -  double Py,
                                             uv[1], //   -  double Qy,
                                             //
                                             a_SL,
                                             b_SL,
                                             mv[0], //   -  double Hx,
                                             //
                                             tau, iCurrTL, //   -  Index of current time layer.
                                             //
                                             indCurSqOx, //   -  Index of current square by Ox axis.
                                             indCurSqOy, //   -  Index of current square by Oy axis.
                                             //
                                             masOX, //   -  Massive of OX steps. Dimension = numOfOXSt +1.
                                             numOfOXSt, //   -  Number of OX steps.
                                             //
                                             masOY, //   -  Massive of OY steps. Dimension = numOfOYSt +1.
                                             numOfOYSt, //   -  Number of OY steps.
                                             //
                                             rhoInPrevTL_asV);


            integ += buf_D;
        }
    }


    //   B. Under rectangle. Need to be cheking.

    if (wMvI == 1)
    {
        if (indCurSqOx[0] == indRB[0])
        {
            Hx = rb;
        }

        if (indCurSqOx[0] < indRB[0])
        {
            if (indCurSqOx[1] >= 0)
            {
                Hx = masOX[ indCurSqOx[1] ];
            }

            if (indCurSqOx[1] < 0)
            {
                Hx = h * indCurSqOx[1];
            }
        }

        buf_D = integUnderRectAng_OneCell(bv[1], //   -  double Py,
                                          uv[1], //   -  double Qy,
                                          //
                                          mv[0], //   -  double Gx,
                                          Hx, //   -  double Hx,
                                          //
                                          tau, iCurrTL, //   -  Index of current time layer.
                                          //
                                          indCurSqOx, //   -  Index of current square by Ox axis.
                                          indCurSqOy, //   -  Index of current square by Oy axis.
                                          //
                                          masOX, //   -  Massive of OX steps. Dimension = numOfOXSt +1.
                                          numOfOXSt, //   -  Number of OX steps.
                                          //
                                          masOY, //   -  Massive of OY steps. Dimension = numOfOYSt +1.

                                          rhoInPrevTL_asV);

        integ += buf_D;
    }

    //   Second step: from "masOX[ indCurSqOx[1] ]" to "rb" by iteration.


    indCurSqOxToCh[0] = indCurSqOx[0] + 1;
    indCurSqOxToCh[1] = indCurSqOxToCh[0] + 1;

    for (j = indCurSqOx[0] + 1; j < indRB[0] + 1; j++)
    {
        //   If this is first cell we should integrate under triangle only.

        if (indCurSqOxToCh[1] > 0)
        {
            Gx = masOX[ indCurSqOxToCh[0] ];
            Hx = masOX[ indCurSqOxToCh[1] ];
        }


        if (indCurSqOxToCh[1] <= 0)
        {
            Gx = h * indCurSqOxToCh[0];
            Hx = h * indCurSqOxToCh[1];
        }


        if (j == indRB[0])
        {
            Hx = rb;
        }


        buf_D = integUnderRectAng_OneCell(bv[1], //   -  double Py,
                                          uv[1], //   -  double Qy,
                                          //
                                          Gx, //   -  double Gx,
                                          Hx, //   -  double Hx,
                                          //
                                          tau, iCurrTL, //   -  Index of current time layer.
                                          //
                                          indCurSqOxToCh, //   -  Index of current square by Ox axis.
                                          indCurSqOy, //   -  Index of current square by Oy axis.
                                          //
                                          masOX, //   -  Massive of OX steps. Dimension = numOfOXSt +1.
                                          numOfOXSt, //   -  Number of OX steps.
                                          //
                                          masOY, //   -  Massive of OY steps. Dimension = numOfOYSt +1.

                                          rhoInPrevTL_asV);

        integ += buf_D;


        indCurSqOxToCh[0] += 1;
        indCurSqOxToCh[1] = indCurSqOxToCh[0] + 1;
    }

    return integ;
}

double integUnderRigAngTr_BottLeft(
                                   double tau,
                                   int iCurrTL, //   -  Index of current time layer.
                                   //
                                   double *bv,
                                   double *uv,
                                   //
                                   const double *masOX, //   -  Massive of OX steps. Dimension = numOfOXSt +1.
                                   int numOfOXSt, //   -  Number of OX steps.
                                   //
                                   const double *masOY, //   -  Massive of OY steps. Dimension = numOfOYSt +1.
                                   int numOfOYSt, //   -  Number of OY steps.
                                   //
                                   double *rhoInPrevTL_asV)
{
    double trPC[2]; //   -  Travel point current;
    int wTrPCI = 0; //   -  Where travel point current is?
    double trPN[2]; //   -  Travel point next;
    int wTrPNI = 0; //   -  Where travel point next is?
    double ang; //   -  Angle of slant line. Should be greater zero.
    int indCurSqOx[2], indCurSqOy[2]; //   -  Index of current square by Ox and Oy axes.
    int indRB[2]; //   -  Index of right boundary.
    double distOx, distOy; //   -  Distance to near Ox and Oy straight lines.
    bool isTrDone = false; //   -  Is travel done.
    double hx = masOX[1] - masOX[0];
    double hy = masOY[1] - masOY[0];
    double integOfBottTr = 0.; //   -  Value which we are computing.
    double buf_D;
    //   Initial data.
    trPC[0] = bv[0];
    trPC[1] = bv[1];
    if ((fabs(bv[0] - uv[0])) < 1.e-12)
    {
        //   This triangle has very small width. I guess further computation isn't correct.
        return fabs(bv[0] - uv[0]);
    }
    ang = (uv[1] - bv[1]) / (bv[0] - uv[0]);
    if (fabs(ang) < 1.e-12)
    {
        //   This triangle has very small height. I guess further computation isn't correct.
        return fabs(ang);
    }
    indCurSqOx[0] = (int) ((trPC[0] - 1.e-14) / hx); //   -  If trPC[0] is in grid edge I want it will be between in the left side of indCurSqOx[1].
    if ((trPC[0] - 1.e-14) <= 0)
    {
        indCurSqOx[0] -= 1; //   -  The case when "trPC[0]" ia negative.
    }
    indCurSqOx[1] = indCurSqOx[0] + 1; //   -  It's important only in rare case then trPC is in grid edge.
    indRB[0] = indCurSqOx[0];
    indRB[1] = indRB[0] + 1;
    indCurSqOy[0] = (int) ((trPC[1] + 1.e-14) / hy); //   -  If trPC[1] is in grid edge I want it will be between indCurSqOx[0] and indCurSqOx[1].
    if ((trPC[1] + 1.e-14) <= 0)
    {
        indCurSqOy[0] -= 1; //   -  The case when "trPC[0]" ia negative.
    }
    indCurSqOy[1] = indCurSqOy[0] + 1; //   -  It's important only in rare case then trPC is in grid edge.
    if (indCurSqOx[0] >= 0)
    {
        distOx = trPC[0] - masOX[ indCurSqOx[0] ];
    }
    if (indCurSqOx[0] < 0)
    {
        distOx = fabs(trPC[0] - hx * indCurSqOx[0]);
    }
    if (indCurSqOy[1] >= 0)
    {
        distOy = masOY[ indCurSqOy[1] ] - trPC[1];
    }
    if (indCurSqOy[1] < 0)
    {
        distOy = fabs(hy * indCurSqOy[1] - trPC[1]);
    }
    do
    {
        //   a. First case.
        if ((distOy / distOx) <= ang)
        {
            //   Across with straight line parallel Ox axis.
            wTrPNI = 1;
            if (indCurSqOy[1] >= 0)
            {
                trPN[1] = masOY[ indCurSqOy[1] ];
            }
            if (indCurSqOy[1] < 0)
            {
                trPN[1] = hy * indCurSqOy[1];
            }
            trPN[0] = bv[0] - (trPN[1] - bv[1]) / ang;
        }
        //   b. Second case.
        if ((distOy / distOx) > ang)
        {
            //   Across with straight line parallel Oy axis.
            wTrPNI = 2;
            if (indCurSqOx[0] >= 0)
            {
                trPN[0] = masOX[ indCurSqOx[0] ];
            }
            if (indCurSqOx[0] < 0)
            {
                trPN[0] = hx * indCurSqOx[0];
            }
            trPN[1] = bv[1] - ang * (trPN[0] - bv[0]);
        }
        //   c. Cheking.
        if (trPN[0] < (uv[0] + 1.e-14))
        {
            trPN[0] = uv[0];
            trPN[1] = uv[1];
            isTrDone = true;
            wTrPNI = 0;
        }
        //   d. Integration.
        buf_D = integOfChan_SLLeftSd(
                                     tau, iCurrTL, //   -  Index of current time layer.
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
                                     masOX, //   -  Massive of OX steps. Dimension = numOfOXSt +1.
                                     numOfOXSt, //   -  Number of OX steps.
                                     //
                                     masOY, //   -  Massive of OY steps. Dimension = numOfOYSt +1.
                                     numOfOYSt, //   -  Number of OY steps.
                                     //
                                     rhoInPrevTL_asV);
        integOfBottTr = integOfBottTr + buf_D;
        //   e. Updating.
        if (isTrDone == false)
        {
            //   We will compute more. We need to redefine some values.
            wTrPCI = wTrPNI;
            trPC[0] = trPN[0];
            trPC[1] = trPN[1];
            if (wTrPNI == 1)
            {
                indCurSqOy[0] += 1;
                indCurSqOy[1] += 1;
            }
            if (wTrPNI == 2)
            {
                indCurSqOx[0] -= 1;
                indCurSqOx[1] -= 1;
            }
            if (indCurSqOx[0] >= 0)
            {
                distOx = trPC[0] - masOX[ indCurSqOx[0] ];
            }
            if (indCurSqOx[0] < 0)
            {
                distOx = fabs(trPC[0] - hx * indCurSqOx[0]);
            }
            if (indCurSqOy[1] >= 0)
            {
                distOy = masOY[ indCurSqOy[1] ] - trPC[1];
            }
            if (indCurSqOy[1] < 0)
            {
                distOy = fabs(hy * indCurSqOy[1] - trPC[1]);
            }
        }
    }
    while (!isTrDone);
    return integOfBottTr;
}

double integUnderRigAngTr_BottRight(
                                    double tau,
                                    int iCurrTL, //   -  Index of current time layer.
                                    //
                                    double *bv,
                                    double *uv,
                                    //
                                    const double *masOX, //   -  Massive of OX steps. Dimension = numOfOXSt +1.
                                    int numOfOXSt, //   -  Number of OX steps.
                                    //
                                    const double *masOY, //   -  Massive of OY steps. Dimension = numOfOYSt +1.
                                    int numOfOYSt, //   -  Number of OY steps.
                                    //
                                    double *rhoInPrevTL_asV)
{
    double trPC[2]; //   -  Travel point current;
    int wTrPCI = 0; //   -  Where travel point current is?
    double trPN[2]; //   -  Travel point next;
    int wTrPNI = 0; //   -  Where travel point next is?
    double ang; //   -  Angle of slant line. Should be greater zero.
    int indCurSqOx[2], indCurSqOy[2]; //   -  Index of current square by Ox and Oy axes.
    int indLB[2]; //   -  Index of left boundary.
    double distOx, distOy; //   -  Distance to near Ox and Oy straight lines.
    bool isTrDone = false; //   -  Is travel done.
    double hx = masOX[1] - masOX[0];
    double hy = masOY[1] - masOY[0];
    double integOfBottTr = 0.; //   -  Value which we are computing.
    double buf_D;


    trPC[0] = bv[0];
    trPC[1] = bv[1];
    if ((fabs(bv[0] - uv[0])) < 1.e-12) return fabs(bv[0] - uv[0]);

    ang = (uv[1] - bv[1]) / (uv[0] - bv[0]);
    if (fabs(ang) < 1.e-12) return fabs(ang);

    indCurSqOx[0] = (int) ((trPC[0] + 1.e-14) / hx); //   -  If trPC[0] is in grid edge I want it will be between in the right side.

    if ((trPC[0] + 1.e-14) <= 0) indCurSqOx[0] -= 1; //   -  The case when "trPC[0]" is negative.

    indCurSqOx[1] = indCurSqOx[0] + 1; //   -  It's important only in rare case then trPC is in grid edge.
    indLB[0] = indCurSqOx[0];
    indLB[1] = indLB[0] + 1;
    indCurSqOy[0] = (int) ((trPC[1] + 1.e-14) / hy); //   -  If trPC[1] is in grid edge I want it will be in the upper side.
    if ((trPC[1] + 1.e-14) <= 0)
    {
        indCurSqOy[0] -= 1; //   -  The case when "trPC[0]" ia negative.
    }
    indCurSqOy[1] = indCurSqOy[0] + 1; //   -  It's important only in rare case then trPC is in grid edge.

    if (indCurSqOx[1] >= 0)
    {
        distOx = fabs(masOX[ indCurSqOx[1] ] - trPC[0]);
    }
    if (indCurSqOx[1] < 0)
    {
        distOx = fabs(hx * indCurSqOx[1] - trPC[0]);
    }
    if (indCurSqOy[1] >= 0)
    {
        distOy = fabs(masOY[ indCurSqOy[1] ] - trPC[1]);
    }
    if (indCurSqOy[1] < 0)
    {
        distOy = fabs(hy * indCurSqOy[1] - trPC[1]);
    }
    do
    {
        //   a. First case.
        if ((distOy / distOx) <= ang)
        {
            //   Across with straight line parallel Ox axis.
            wTrPNI = 1;
            if (indCurSqOy[1] >= 0)
            {
                trPN[1] = masOY[ indCurSqOy[1] ];
            }
            if (indCurSqOy[1] < 0)
            {
                trPN[1] = hy * indCurSqOy[1];
            }
            trPN[0] = bv[0] + (trPN[1] - bv[1]) / ang;
        }
        //   b. Second case.
        if ((distOy / distOx) > ang)
        {
            //   Across with straight line parallel Oy axis.
            wTrPNI = 2;
            if (indCurSqOx[1] >= 0)
            {
                trPN[0] = masOX[ indCurSqOx[1] ];
            }
            if (indCurSqOx[1] < 0)
            {
                trPN[0] = hx * indCurSqOx[1];
            }
            trPN[1] = bv[1] + ang * (trPN[0] - bv[0]);
        }
        //   c. Cheking.
        if (trPN[0] > (uv[0] - 1.e-14))
        {
            //   -  Without "fabs"!!!
            trPN[0] = uv[0];
            trPN[1] = uv[1];
            isTrDone = true;
            wTrPNI = 0;
        }
        //   d. Integration.
        buf_D = integOfChan_SLRightSd(
                                      tau, iCurrTL, //   -  Index of current time layer.
                                      //
                                      trPC, wTrPCI, //   -  double *bv,
                                      trPN, wTrPNI, //   -  double *uv,
                                      //
                                      indCurSqOx, //   -  Indices where trPC and trPN are.
                                      //
                                      bv[0], indLB, //   -  double lb  =  Left boundary by Ox.
                                      //
                                      indCurSqOy, //   -  Index of current square by Oy axis.
                                      //
                                      masOX, //   -  Massive of OX steps. Dimension = numOfOXSt +1.
                                      numOfOXSt, //   -  Number of OX steps.
                                      //
                                      masOY, //   -  Massive of OY steps. Dimension = numOfOYSt +1.
                                      numOfOYSt, //   -  Number of OY steps.
                                      //
                                      rhoInPrevTL_asV);
        integOfBottTr = integOfBottTr + buf_D;
        //   e. Updating.
        if (isTrDone == false)
        {
            //   We will compute more. We need to redefine some values.
            wTrPCI = wTrPNI;
            trPC[0] = trPN[0];
            trPC[1] = trPN[1];
            if (wTrPNI == 1)
            {
                indCurSqOy[0] += 1;
                indCurSqOy[1] += 1;
            }
            if (wTrPNI == 2)
            {
                indCurSqOx[0] += 1;
                indCurSqOx[1] += 1;
            }
            if (indCurSqOx[1] >= 0)
            {
                distOx = fabs(masOX[ indCurSqOx[1] ] - trPC[0]);
            }
            if (indCurSqOx[1] < 0)
            {
                distOx = fabs(hx * indCurSqOx[1] - trPC[0]);
            }
            if (indCurSqOy[1] >= 0)
            {
                distOy = fabs(masOY[ indCurSqOy[1] ] - trPC[1]);
            }
            if (indCurSqOy[1] < 0)
            {
                distOy = fabs(hy * indCurSqOy[1] - trPC[1]);
            }
        }
    }
    while (!isTrDone);
    return integOfBottTr;
}

double integ_under_bott_triangle(double tau,
                                 int curr_tl,
                                 double *lv, //   -  Left vertex of Bottom triangle.
                                 double *rv, //   -  Right vertex of Bottom triangle.
                                 double *bv, //   -  Bottom vertex of Bottom triangle.
                                 const double *ox,
                                 int ox_length,
                                 const double *oy,
                                 int oy_length,
                                 double *prev_density)
{
    double result = 0.;
    if (bv[0] <= lv[0])
    {
        result = integUnderRigAngTr_BottRight(tau, curr_tl, bv, rv, ox, ox_length, oy, oy_length, prev_density);
        result -= integUnderRigAngTr_BottRight(tau, curr_tl, bv, lv, ox, ox_length, oy, oy_length, prev_density);
        return result;
    }
    if (bv[0] > lv[0] && bv[0] < rv[0])
    {
        result = integUnderRigAngTr_BottLeft(tau, curr_tl, bv, lv, ox, ox_length, oy, oy_length, prev_density);
        result += integUnderRigAngTr_BottRight(tau, curr_tl, bv, rv, ox, ox_length, oy, oy_length, prev_density);
        return result;
    }
    if (bv[0] >= rv[0])
    {
        result = integUnderRigAngTr_BottLeft(tau, curr_tl, bv, lv, ox, ox_length, oy, oy_length, prev_density);
        result -= integUnderRigAngTr_BottLeft(tau, curr_tl, bv, rv, ox, ox_length, oy, oy_length, prev_density);
        return result;
    }
    return result;
}

double integUnderRigAngTr_UppLeft(double tau,
                                  int iCurrTL,
                                  double *bv,
                                  double *uv,
                                  const double *masOX,
                                  int numOfOXSt,
                                  const double *masOY,
                                  int numOfOYSt,
                                  double *rhoInPrevTL_asV)
{
    double trPC[2]; //   -  Travel point current;
    int wTrPCI = 0; //   -  Where travel point current is?
    double trPN[2]; //   -  Travel point next;
    int wTrPNI = 0; //   -  Where travel point next is?
    double ang; //   -  Angle of slant line. Should be greater zero.
    int indCurSqOx[2], indCurSqOy[2]; //   -  Index of current square by Ox and Oy axes.
    int indRB[2]; //   -  Index of right boundary.
    double distOx, distOy; //   -  Distance to near Ox and Oy straight lines.
    bool isTrDone = false; //   -  Is travel done.
    double hx = masOX[1] - masOX[0];
    double hy = masOY[1] - masOY[0];
    double integOfUppTr = 0.; //   -  Value which we are computing.
    double buf_D;
    //   Initial data.
    trPC[0] = bv[0];
    trPC[1] = bv[1];
    if ((fabs(bv[0] - uv[0])) < 1.e-12) return fabs(bv[0] - uv[0]);

    ang = (uv[1] - bv[1]) / (uv[0] - bv[0]);
    if (fabs(ang) < 1.e-12) return fabs(ang);

    //   The follow equations are quite important.
    indCurSqOx[0] = (int) ((trPC[0] + 1.e-14) / hx); //   -  If trPC[0] is in grid edge I want it will be in the right side.
    if ((trPC[0] + 1.e-14) <= 0)
    {
        indCurSqOx[0] -= 1; //   -  The case when "trPC[0]" ia negative.
    }
    indCurSqOx[1] = indCurSqOx[0] + 1; //   -  It's important only in rare case then trPC is in grid edge.
    indCurSqOy[0] = (int) ((trPC[1] + 1.e-14) / hy); //   -  If trPC[1] is in grid edge I want it will be in the upper square.
    if ((trPC[1] + 1.e-14) <= 0)
    {
        indCurSqOy[0] -= 1; //   -  The case when "trPC[0]" ia negative.
    }
    indCurSqOy[1] = indCurSqOy[0] + 1;
    indRB[0] = (int) ((uv[0] - 1.e-14) / hy); //   -  If uv[0] is in grid edge I want it will be in the left side.
    if ((uv[0] - 1.e-14) <= 0)
    {
        indRB[0] -= 1; //   -  The case when "trPC[0]" ia negative.
    }
    indRB[1] = indRB[0] + 1;
    if (indCurSqOx[1] >= 0)
    {
        distOx = masOX[ indCurSqOx[1] ] - trPC[0];
    }
    if (indCurSqOx[1] < 0)
    {
        distOx = fabs(hx * indCurSqOx[1] - trPC[0]);
    }
    if (indCurSqOy[1] >= 0)
    {
        distOy = masOY[ indCurSqOy[1] ] - trPC[1];
    }
    if (indCurSqOy[1] < 0)
    {
        distOy = fabs(hy * indCurSqOy[1] - trPC[1]);
    }
    do
    {
        //   a. First case.
        if ((distOy / distOx) <= ang)
        {
            //   Across with straight line parallel Ox axis.
            wTrPNI = 1;
            if (indCurSqOy[1] >= 0)
            {
                trPN[1] = masOY[ indCurSqOy[1] ];
            }
            if (indCurSqOy[1] < 0)
            {
                trPN[1] = hy * indCurSqOy[1];
            }
            trPN[0] = bv[0] + (trPN[1] - bv[1]) / ang;
        }
        //   b. Second case.
        if ((distOy / distOx) > ang)
        {
            //   Across with straight line parallel Oy axis.
            wTrPNI = 2;
            if (indCurSqOx[1] >= 0)
            {
                trPN[0] = masOX[ indCurSqOx[1] ];
            }
            if (indCurSqOx[1] < 0)
            {
                trPN[0] = hx * indCurSqOx[1];
            }
            trPN[1] = bv[1] + ang * (trPN[0] - bv[0]);
        }
        //   c. Cheking.
        if (trPN[0] > (uv[0] - 1.e-14))
        {
            trPN[0] = uv[0];
            trPN[1] = uv[1];
            isTrDone = true;
            wTrPNI = 0;
        }
        //   d. Integration.
        buf_D = integOfChan_SLLeftSd(
                                     tau, iCurrTL, //   -  Index of current time layer.
                                     //
                                     trPC, wTrPCI, //   -  double *bv,
                                     trPN, wTrPNI, //   -  double *uv,
                                     //
                                     indCurSqOx, //   -  Indices where trPC and trPN are.
                                     //
                                     uv[0], indRB, //   -  double rb  =  Right boundary by Ox.
                                     //
                                     indCurSqOy, //   -  Index of current square by Oy axis.
                                     //
                                     masOX, //   -  Massive of OX steps. Dimension = numOfOXSt +1.
                                     numOfOXSt, //   -  Number of OX steps.
                                     //
                                     masOY, //   -  Massive of OY steps. Dimension = numOfOYSt +1.
                                     numOfOYSt, //   -  Number of OY steps.
                                     //
                                     rhoInPrevTL_asV);
        integOfUppTr = integOfUppTr + buf_D;
        //   e. Updating.
        if (isTrDone == false)
        {
            //   We will compute more. We need to redefine some values.
            wTrPCI = wTrPNI;
            trPC[0] = trPN[0];
            trPC[1] = trPN[1];
            if (wTrPNI == 1)
            {
                indCurSqOy[0] += 1;
                indCurSqOy[1] += 1;
            }
            if (wTrPNI == 2)
            {
                indCurSqOx[0] += 1;
                indCurSqOx[1] += 1;
            }
            if (indCurSqOx[1] >= 0)
            {
                distOx = fabs(masOX[ indCurSqOx[1] ] - trPC[0]);
            }
            if (indCurSqOx[1] < 0)
            {
                distOx = fabs(hx * indCurSqOx[1] - trPC[0]);
            }
            if (indCurSqOy[1] >= 0)
            {
                distOy = fabs(masOY[ indCurSqOy[1] ] - trPC[1]);
            }
            if (indCurSqOy[1] < 0)
            {
                distOy = fabs(hy * indCurSqOy[1] - trPC[1]);
            }
        }
    }
    while (!isTrDone);
    return integOfUppTr;
}

double integUnderRigAngTr_UppRight(double tau,
                                   int iCurrTL,
                                   double *bv,
                                   double *uv,
                                   const double *masOX,
                                   int numOfOXSt,
                                   const double *masOY,
                                   int numOfOYSt,
                                   double *rhoInPrevTL_asV)
{
    double trPC[2]; //   -  Travel point current;
    int wTrPCI = 0; //   -  Where travel point current is?
    double trPN[2]; //   -  Travel point next;
    int wTrPNI = 0; //   -  Where travel point next is?
    double ang; //   -  Angle of slant line. Should be greater zero.
    int indCurSqOx[2], indCurSqOy[2]; //   -  Index of current square by Ox and Oy axes.
    int indLB[2]; //   -  Index of left boundary.
    double distOx, distOy; //   -  Distance to near Ox and Oy straight lines.
    bool isTrDone = false; //   -  Is travel done.
    double hx = masOX[1] - masOX[0];
    double hy = masOY[1] - masOY[0];
    double integOfUppTr = 0.; //   -  Value which we are computing.
    double buf_D;
    //   Initial data.
    trPC[0] = bv[0];
    trPC[1] = bv[1];
    if ((fabs(bv[0] - uv[0])) < 1.e-12)
    {
        //   This triangle has very small width. I guess further computation isn't correct.
        return fabs(bv[0] - uv[0]);
    }
    ang = (uv[1] - bv[1]) / (bv[0] - uv[0]);
    if (fabs(ang) < 1.e-12)
    {
        //   This triangle has very small height. I guess further computation isn't correct.
        return fabs(ang);
    }
    indCurSqOx[0] = (int) ((trPC[0] - 1.e-14) / hx); //   -  If trPC[0] is in grid edge I want it will be between in the left side.
    if ((trPC[0] - 1.e-14) <= 0)
    {
        indCurSqOx[0] -= 1; //   -  The case when "trPC[0]" ia negative.
    }
    indCurSqOx[1] = indCurSqOx[0] + 1; //   -  It's important only in rare case then trPC is in grid edge.
    indLB[0] = (int) ((uv[0] + 1.e-14) / hx);
    if ((uv[0] + 1.e-14) <= 0)
    {
        indLB[0] -= 1; //   -  The case when "trPC[0]" ia negative.
    }
    indLB[1] = indLB[0] + 1;
    indCurSqOy[0] = (int) ((trPC[1] + 1.e-14) / hy); //   -  If trPC[1] is in grid edge I want it will be in the upper side.
    if ((trPC[1] + 1.e-14) <= 0)
    {
        indCurSqOy[0] -= 1; //   -  The case when "trPC[0]" ia negative.
    }
    indCurSqOy[1] = indCurSqOy[0] + 1; //   -  It's important only in rare case then trPC is in grid edge.
    if (indCurSqOx[0] >= 0)
    {
        distOx = fabs(trPC[0] - masOX[ indCurSqOx[0] ]);
    }
    if (indCurSqOx[0] < 0)
    {
        distOx = fabs(trPC[0] - hx * indCurSqOx[0]);
    }
    if (indCurSqOy[1] >= 0)
    {
        distOy = fabs(masOY[ indCurSqOy[1] ] - trPC[1]);
    }
    if (indCurSqOy[1] < 0)
    {
        distOy = fabs(hy * indCurSqOy[1] - trPC[1]);
    }
    do
    {
        //   a. First case.
        if ((distOy / distOx) <= ang)
        {
            //   Across with straight line parallel Ox axis.
            wTrPNI = 1;
            if (indCurSqOy[1] >= 0)
            {
                trPN[1] = masOY[ indCurSqOy[1] ];
            }
            if (indCurSqOy[1] < 0)
            {
                trPN[1] = hy * indCurSqOy[1];
            }
            trPN[0] = bv[0] - (trPN[1] - bv[1]) / ang;
        }
        //   b. Second case.
        if ((distOy / distOx) > ang)
        {
            //   Across with straight line parallel Oy axis.
            wTrPNI = 2;
            if (indCurSqOx[0] >= 0)
            {
                trPN[0] = masOX[ indCurSqOx[0] ];
            }
            if (indCurSqOx[0] < 0)
            {
                trPN[0] = hx * indCurSqOx[0];
            }
            trPN[1] = bv[1] - ang * (trPN[0] - bv[0]);
        }
        //   c. Cheking.
        if (trPN[0] < (uv[0] + 1.e-14))
        {
            trPN[0] = uv[0];
            trPN[1] = uv[1];
            isTrDone = true;
            wTrPNI = 0;
        }
        //   d. Integration.
        buf_D = integOfChan_SLRightSd(
                                      tau, iCurrTL, //   -  Index of current time layer.
                                      //
                                      trPC, wTrPCI, //   -  double *bv,
                                      trPN, wTrPNI, //   -  double *uv,
                                      //
                                      indCurSqOx, //   -  Indices where trPC and trPN are.
                                      //
                                      uv[0], indLB, //   -  double lb  =  Left boundary by Ox.
                                      //
                                      indCurSqOy, //   -  Index of current square by Oy axis.
                                      //
                                      masOX, //   -  Massive of OX steps. Dimension = numOfOXSt +1.
                                      numOfOXSt, //   -  Number of OX steps.
                                      //
                                      masOY, //   -  Massive of OY steps. Dimension = numOfOYSt +1.
                                      numOfOYSt, //   -  Number of OY steps.
                                      //
                                      rhoInPrevTL_asV);
        integOfUppTr = integOfUppTr + buf_D;
        //   e. Updating.
        if (isTrDone == false)
        {
            //   We will compute more. We need to redefine some values.
            wTrPCI = wTrPNI;
            trPC[0] = trPN[0];
            trPC[1] = trPN[1];
            if (wTrPNI == 1)
            {
                indCurSqOy[0] += 1;
                indCurSqOy[1] += 1;
            }
            if (wTrPNI == 2)
            {
                indCurSqOx[0] -= 1;
                indCurSqOx[1] -= 1;
            }
            if (indCurSqOx[0] >= 0)
            {
                distOx = fabs(trPC[0] - masOX[ indCurSqOx[0] ]);
            }
            if (indCurSqOx[0] < 0)
            {
                distOx = fabs(trPC[0] - hx * indCurSqOx[0]);
            }
            if (indCurSqOy[1] >= 0)
            {
                distOy = fabs(masOY[ indCurSqOy[1] ] - trPC[1]);
            }
            if (indCurSqOy[1] < 0)
            {
                distOy = fabs(hy * indCurSqOy[1] - trPC[1]);
            }
        }
    }
    while (!isTrDone);
    return integOfUppTr;
}

double integ_under_upper_triangle(double tau,
                                  int curr_tl,
                                  double *lv, //   -  Left vertex of Upper triangle.
                                  double *rv, //   -  Right vertex of Upper triangle.
                                  double *uv, //   -  Upper vertex of Upper triangle.
                                  const double *ox,
                                  int ox_length,
                                  const double *oy,
                                  int oy_length,
                                  double *prev_density)
{
    double result = 0.;
    if (uv[0] <= lv[0])
    {
        result = integUnderRigAngTr_UppRight(tau, curr_tl, rv, uv, ox, ox_length, oy, oy_length, prev_density);
        result -= integUnderRigAngTr_UppRight(tau, curr_tl, lv, uv, ox, ox_length, oy, oy_length, prev_density);
        return result;
    }
    if (uv[0] > lv[0] && uv[0] < rv[0])
    {
        result = integUnderRigAngTr_UppLeft(tau, curr_tl, lv, uv, ox, ox_length, oy, oy_length, prev_density);
        result += integUnderRigAngTr_UppRight(tau, curr_tl, rv, uv, ox, ox_length, oy, oy_length, prev_density);
        return result;
    }
    if (uv[0] >= rv[0])
    {
        result = integUnderRigAngTr_UppLeft(tau, curr_tl, lv, uv, ox, ox_length, oy, oy_length, prev_density);
        result -= integUnderRigAngTr_UppLeft(tau, curr_tl, rv, uv, ox, ox_length, oy, oy_length, prev_density);
        return result;
    }
    return result;
}

double wall_integ_under_uniform_triangle(double tau,
                                         int curr_tl,
                                         point_t *a,
                                         point_t *b,
                                         point_t *c,
                                         const double *ox,
                                         int ox_length,
                                         const double *oy,
                                         int oy_length,
                                         double *prev_density)
{
    return 0;
}

double integ_under_uniform_triangle(double tau,
                                    int tl,
                                    point_t *x,
                                    point_t *y,
                                    point_t *z,
                                    const double *ox,
                                    int ox_length,
                                    const double *oy,
                                    int oy_length,
                                    double *prev_density)
{
    const double MIN_VALUE = 1.e-12;
    double bv[2], mv[2], uv[2], ip[2]; //   -  Bottom, middle and upper vertices + intersection point
    
    bv[0] = x->x;
    bv[1] = x->y;
    mv[0] = y->x;
    mv[1] = y->y;
    uv[0] = z->x;
    uv[1] = z->y;

    //   2. I want to compute intersection point.
    //   2.a Let's compute line coefficients between "bv" and "uv" vertices.
    //   a * x  +  b * y  = c.
    double a = z->y - x->y;
    if (fabs(a) < MIN_VALUE) return MIN_VALUE;
    double b = x->x - z->x;
    double c = b * x->y + a * x->x;
    ip[1] = mv[1];
    ip[0] = (c - b * mv[1]) / a;

    //   3. Возможны 2 случая расположения точки перечеения относительно средней
    //   слева или справа.

    if (mv[0] >= ip[0]) // средняя точка справа от точки пересечения
    { // обменяем местами  X координаты, чтобы использовать один код для расчета
        double tx = mv[0];
        mv[0] = ip[0];
        ip[0] = tx;
    }

    return integ_under_bott_triangle(tau, tl,
                                     mv, ip, bv,
                                     ox, ox_length,
                                     oy, oy_length,
                                     prev_density)
            + integ_under_upper_triangle(tau, tl,
                                         mv, ip, uv,
                                         ox, ox_length,
                                         oy, oy_length,
                                         prev_density);
}

inline double u_function(double x, double y)
{
    return PAR_B * y * (1. - y) * (C_pi / 2. + atan(-x));
}

inline double v_function(
                         double lbDom,
                         double rbDom,
                         double bbDom,
                         double ubDom,
                         double t, double x, double y)
{
    return atan((x - lbDom) * (x - rbDom) * (1. + t) / 10. * (y - ubDom) * (y - bbDom));
}

double f_function(
                  double lbDom, //   -  Left and right boundaries of rectangular domain.
                  double rbDom,
                  //
                  double bbDom, //   -  Bottom and upper boundaries of rectangular domain.
                  double ubDom,
                  //
                  double tau,
                  int iCurrTL, //   -  Index of current time layer.
                  //
                  int iOfOXN, //   -  Index of current OX node.
                  const double *masOX, //   -  Massive of OX steps. Dimension = numOfOXSt

                  int iOfOYN, //   -  Index of current OY node.
                  const double *masOY) //   -  Number of OY steps (segments).
{
    //printf("\ncpu f function \n");
    double t = tau * iCurrTL;
    //printf("cpu t = %f\n", t);
    double x = masOX[ iOfOXN ];
    //printf("cpu x = %f\n", x);
    double y = masOY[ iOfOYN ];
    //printf("cpu y = %f\n", y);
    double arg_v = (x - lbDom) * (x - rbDom) * (1. + t) / 10. * (y - ubDom) * (y - bbDom);
    //printf("cpu arg_v = %f\n", arg_v);
    double rho, dRhoDT, dRhoDX, dRhoDY;
    double u, duDX;
    double v, dvDY;
    rho = analytical_solution(t, x, y);
    //  printf("cpu rho = %f\n", rho);
    dRhoDT = x * y * cos(t * x * y);
    //  printf("cpu dRhoDT = %f\n", dRhoDT);
    dRhoDX = t * y * cos(t * x * y);
    //  printf("cpu dRhoDX = %f\n", dRhoDX);
    dRhoDY = t * x * cos(t * x * y);
    //  printf("cpu dRhoDY = %f\n", dRhoDY);
    u = u_function(x, y);
    //  printf("cpu u = %f\n", u);
    duDX = -PAR_B * y * (1. - y) / (1. + x * x);
    //  printf("cpu duDX = %f\n", duDX);
    v = v_function(lbDom, rbDom, bbDom, ubDom, t, x, y);
    //  printf("cpu v = %f\n", v);
    dvDY = (x - lbDom) * (x - rbDom) * (1. + t) / 10. * (y - bbDom + y - ubDom);
    //  printf("cpu dvDY 1 = %f\n", dvDY);
    dvDY = dvDY / (1. + arg_v * arg_v);
    //  printf("cpu dvDY 2 = %f\n", dvDY);
    double res = dRhoDT + rho * duDX + u * dRhoDX + rho * dvDY + v * dRhoDY;
    //  printf("cpu res = %f\n", res);
    return res;
}

quad_type compute_coordinate_on_prev_layer(double lb,
                                           double rb,
                                           double bb,
                                           double ub,
                                           double tau,
                                           int cur_tl,
                                           int i_ox,
                                           const double *masOX,
                                           int ox_length,
                                           int i_oy,
                                           const double *masOY,
                                           int oy_length,
                                           point_t *alpha, point_t *beta, point_t *gamma, point_t *theta)
{
    //   1. First of all let's compute coordinates of square vertexes.
    //  OX:
    if (i_ox == 0)
    {
        alpha->x = masOX[ i_ox ];
        beta->x = (masOX[i_ox] + masOX[i_ox + 1]) / 2.;
        gamma->x = (masOX[i_ox] + masOX[i_ox + 1]) / 2.;
        theta->x = masOX[ i_ox ];
    }
    else if (i_ox == ox_length)
    {
        alpha->x = (masOX[i_ox - 1] + masOX[i_ox]) / 2.;
        beta->x = masOX[ i_ox ];
        gamma->x = masOX[ i_ox ];
        theta->x = (masOX[i_ox - 1] + masOX[i_ox]) / 2.;
    }
    else
    {
        alpha->x = (masOX[i_ox - 1] + masOX[i_ox]) / 2.;
        beta->x = (masOX[i_ox + 1] + masOX[i_ox]) / 2.;
        gamma->x = (masOX[i_ox + 1] + masOX[i_ox]) / 2.;
        theta->x = (masOX[i_ox - 1] + masOX[i_ox]) / 2.;
    }

    //  OY:
    if (i_oy == 0)
    {
        alpha->y = masOY[ i_oy ];
        beta->y = masOY[ i_oy ];
        gamma->y = (masOY[i_oy] + masOY[ i_oy + 1]) / 2.;
        theta->y = (masOY[i_oy] + masOY[ i_oy + 1]) / 2.;
    }
    else if (i_oy == oy_length)
    {
        alpha->y = (masOY[i_oy] + masOY[ i_oy - 1]) / 2.;
        beta->y = (masOY[i_oy] + masOY[ i_oy - 1]) / 2.;
        gamma->y = masOY[ i_oy ];
        theta->y = masOY[ i_oy ];
    }
    else
    {
        alpha->y = (masOY[i_oy] + masOY[ i_oy - 1]) / 2.;
        beta->y = (masOY[i_oy] + masOY[ i_oy - 1]) / 2.;
        gamma->y = (masOY[i_oy] + masOY[ i_oy + 1]) / 2.;
        theta->y = (masOY[i_oy] + masOY[ i_oy + 1]) / 2.;
    }

    double u, v;

    // Now let's compute new coordinates on the previous time level of alpha, beta, gamma, theta points.
    u = u_function(alpha->x, alpha->y);
    v = v_function(lb, rb, bb, ub, tau*cur_tl, alpha->x, alpha->y);
    alpha->x -= tau * u;
    alpha->y -= tau * v;

    u = u_function(beta->x, beta->y);
    v = v_function(lb, rb, bb, ub, tau*cur_tl, beta->x, beta->y);
    beta->x -= tau * u;
    beta->y -= tau * v;

    u = u_function(gamma->x, gamma->y);
    v = v_function(lb, rb, bb, ub, tau*cur_tl, gamma->x, gamma->y);
    gamma->x -= tau * u;
    gamma->y -= tau * v;

    u = u_function(theta->x, theta->y);
    v = v_function(lb, rb, bb, ub, tau*cur_tl, theta->x, theta->y);
    theta->x -= tau * u;
    theta->y -= tau * v;

    point_t intersection = get_intersection_point(alpha, beta, gamma, theta);
    if ((beta->y - intersection.y) * (theta->y - intersection.y) > 0.) return pseudo; // ??
    if ((alpha->x - intersection.x) * (gamma->x - intersection.x) > 0.) return pseudo; // ??
    double product = get_vector_product(alpha, beta, theta); // ?
    if (product < 0.) return pseudo;

    // значит что точка удетела за левую границу
    if (theta->x < 0 ||
            theta->y < 0 ||
            beta->x < 0 ||
            beta->y < 0 ||
            gamma->x < 0 ||
            gamma->y < 0 ||
            alpha->x < 0 ||
            alpha->y < 0)
    {
        return normal;
        //return wall;
    }
    return normal;
}


// Type of quadrangle: 0 - pseudo; 1 - convex; 2 - concave;

quad_type get_quadrangle_type(double lb,
                              double rb,
                              double bb,
                              double ub,
                              double tau,
                              int curr_tl,
                              int i_ox,
                              const double *ox,
                              int ox_length,
                              int i_oy,
                              const double *oy,
                              int oy_length,
                              point_t *t_1_a, //   -  First vertex of first triangle.
                              point_t *t_1_b, //   -  Second vertex of first triangle.
                              point_t *t_1_c, //   -  Third vertex of first triangle.
                              //
                              point_t *t_2_a, //   -  First vertex of second triangle.
                              point_t *t_2_b, //   -  Second vertex of second triangle.
                              point_t *t_2_c) //   -  Third vertex of second triangle.
{
    point_t alpha, beta, gamma, theta; // coordinates on previous time layer

    quad_type type = compute_coordinate_on_prev_layer(lb, rb,
                                                      bb, ub,
                                                      tau, curr_tl,
                                                      i_ox, ox, ox_length,
                                                      i_oy, oy, oy_length, &alpha, &beta, &gamma, &theta);

    // Convex quadrangle DO HAS WRITE anticlockwise vertices sequence order. 
    // It's convex.

    ptcpy(t_1_a, &alpha);
    ptcpy(t_1_b, &beta);
    ptcpy(t_1_c, &gamma);
    ptcpy(t_2_a, &alpha);
    ptcpy(t_2_b, &theta);
    ptcpy(t_2_c, &gamma);

    return type;
}

double compute_value(double lb,
                     double rb,
                     double bb,
                     double ub,
                     double tau,
                     double curr_tl,
                     int i_ox,
                     const double *ox,
                     int ox_length,
                     int i_oy,
                     const double *oy,
                     int oy_length,
                     double *prev_density)
{
    point_t t_1_a, t_1_b, t_1_c, t_2_a, t_2_b, t_2_c;

    quad_type type = get_quadrangle_type(lb, rb, bb, ub, tau, curr_tl,
                                         i_ox, ox, ox_length,
                                         i_oy, oy, oy_length,
                                         &t_1_a, &t_1_b, &t_1_c,
                                         &t_2_a, &t_2_b, &t_2_c);

    if (type != normal && type != wall)
    {
        return -1.;
    }

    // чтобы правилно отработала процедура интегрирования
    // точки должны идти в порядке возрастания y координаты
    sort_by_y(t_1_a, t_1_b, t_1_c);
    sort_by_y(t_2_a, t_2_b, t_2_c);


    // check the type of triangle to select appropriate computation method
    double result = 0.;
    switch (type)
    {
    case wall:
        result += wall_integ_under_uniform_triangle(tau, curr_tl,
                                                    &t_1_a, &t_1_b, &t_1_c,
                                                    ox, ox_length,
                                                    oy, oy_length,
                                                    prev_density);
        result += wall_integ_under_uniform_triangle(tau, curr_tl,
                                                    &t_2_a, &t_2_b, &t_2_c,
                                                    ox, ox_length,
                                                    oy, oy_length,
                                                    prev_density);
        return result;
    case normal:
        result += integ_under_uniform_triangle(tau, curr_tl,
                                               &t_1_a, &t_1_b, &t_1_c,
                                               ox, ox_length,
                                               oy, oy_length,
                                               prev_density);
        result += integ_under_uniform_triangle(tau, curr_tl,
                                               &t_2_a, &t_2_b, &t_2_c,
                                               ox, ox_length,
                                               oy, oy_length,
                                               prev_density);
        return result;
    case concave:
    case convex:
    case pseudo:
        return 0.;
    }
}

void print_params(double b,
                  double lb,
                  double rb,
                  double bb,
                  double ub,
                  double tau,
                  int tl_count,
                  int ox_length,
                  int oy_length)
{
    printf("b = %f\n", b);
    printf("lbDom = %f\n", lb);
    printf("rbDom = %f\n", rb);
    printf("bbDom = %f\n", bb);
    printf("ubDom = %f\n", ub);
    printf("tau = %f\n", tau);
    printf("Time level count = %d\n", tl_count);
    printf("ox length = %d\n", ox_length + 1);
    printf("oy length = %d\n", oy_length + 1);
}

void print_params(int index, int needed_index,
                  double b,
                  double lb,
                  double rb,
                  double bb,
                  double ub,
                  double tau,
                  int tl,
                  int tl_count,
                  int ox_length,
                  int oy_length,
                  int cur_x,
                  int cur_y,
                  double value)
{
    if (index == needed_index)
    {
        printf("index = %d\n", index);
        printf("b = %f\n", b);
        printf("lbDom = %f\n", lb);
        printf("rbDom = %f\n", rb);
        printf("bbDom = %f\n", bb);
        printf("ubDom = %f\n", ub);
        printf("tau = %f\n", tau);
        printf("Time level count = %d\n", tl_count);
        printf("current time level = %d\n", tl);
        printf("current x = %d\n", cur_x);
        printf("current y = %d\n", cur_y);
        printf("ox length = %d\n", ox_length + 1);
        printf("oy length = %d\n", oy_length + 1);
        printf("value = %f\n", value);
    }
}

double get_norm_of_error(double* density, int x_length, int y_length, double* ox,
                         double* oy,
                         double ts_count_mul_steps)
{
    double result = 0.;
    for (int k = 1; k < y_length; ++k)
    {
        for (int j = 1; j < x_length; ++j)
        {
            result += fabs(analytical_solution(ts_count_mul_steps, ox[j], oy[k])
                           - density[ (x_length + 1) * k + j ]);
        }
    }
    double hx = ox[1] - ox[0];
    double hy = oy[1] - oy[0];
    return hx * hy * result;
}

double solve(double lb,
             double rb,
             double bb,
             double ub,
             double tau,
             int time_step_count,
             double *ox,
             int ox_length,
             double *oy,
             int oy_length,
             double *density)
{
    double *prev_density = new double [ (ox_length + 1) * (oy_length + 1) ];
    for (int i_oy = 0; i_oy < oy_length + 1; i_oy++)
    {
        for (int i_ox = 0; i_ox < ox_length + 1; i_ox++)
        {
            prev_density[ (ox_length + 1) * i_oy + i_ox ] = analytical_solution(0., ox[i_ox], oy[i_oy]);
        }
    }

    for (int i_tl = 1; i_tl < time_step_count + 1; i_tl++)
    {
        for (int i = 0; i <= ox_length; i++)
        {
            density[ i ] = init_bound(ox[ i ], bb, tau * i_tl, bottom);
            density[ (ox_length + 1) * oy_length + i ] = init_bound(ox[ i ], ub, tau * i_tl, up);
        }

        for (int i = 0; i <= oy_length; i++)
        {
            density[ (ox_length + 1) * i ] = init_bound(lb, oy[ i ], tau * i_tl, left);
            density[ (ox_length + 1) * i + ox_length ] = init_bound(rb, oy[ i ], tau * i_tl, right);
        }

        for (int i_oy = 1; i_oy < oy_length; i_oy++)
        {
            for (int i_ox = 1; i_ox < ox_length; i_ox++)
            {
                int index = (ox_length + 1) * i_oy + i_ox;

                double value = compute_value(lb, rb,
                                             bb, ub,
                                             tau, i_tl,
                                             i_ox, ox, ox_length,
                                             i_oy, oy, oy_length,
                                             prev_density);
                /*print_params(index, 12,
                             b,
                             lb,
                             rb,
                             bb,
                             ub,
                             tau,
                             i_tl,
                             time_step_count,
                             ox_length,
                             oy_length,
                             i_ox,
                             i_oy,
                             value);*/

                double h = (ox[i_ox + 1] - ox[i_ox - 1]) / 2.;
                value /= h;
                h = (oy[i_oy + 1] - oy[i_oy - 1]) / 2.;
                value /= h;

                double rp = f_function(lb,
                                       rb,
                                       bb,
                                       ub,
                                       tau,
                                       i_tl,
                                       i_ox,
                                       ox,
                                       i_oy,
                                       oy);
                density[ index ] = value + tau * rp;
            }
        }
        memcpy(prev_density, density, (ox_length + 1) * (oy_length + 1) * sizeof (double));
    }

    delete[] prev_density;
    return 0;
}

double *solve(double b,
              double lb,
              double rb,
              double bb,
              double ub,
              double time_step,
              int time_step_count,
              int ox_length,
              int oy_length,
              double* norm)
{
    double *density = new double [ (ox_length + 1) * (oy_length + 1) ];
    double *ox = new double [ ox_length + 1 ];
    double *oy = new double [ oy_length + 1 ];

    for (int i = 0; i <= ox_length; i++)
    {
        ox[i] = lb + i * (rb - lb) / ox_length;
    }

    for (int i = 0; i <= oy_length; i++)
    {
        oy[i] = bb + i * (ub - bb) / oy_length;
    }

    PAR_B = b;
    
    print_params(PAR_B,
                 lb,
                 rb,
                 bb,
                 ub,
                 time_step,
                 time_step_count,
                 ox_length,
                 oy_length);

    solve(lb, rb,
          bb, ub,
          time_step,
          time_step_count,
          ox,
          ox_length,
          oy,
          oy_length,
          density);

    *norm = get_norm_of_error(density, ox_length, oy_length, ox, oy,
                              time_step_count * time_step);
    //  printf("Norm L1 = %f\n", *norm);
    printf("%d x %d wall count = %d\n", ox_length + 1, oy_length + 1, wall_counter);

    delete[] ox;
    delete[] oy;
    return density;
}
