


#include <math.h>

#include "InitAndBoundaryData.h"                      //   -  Initial and Boundary data.



extern double itemOfInteg_2SpecType(
              double Py,
              double Qy,
              //
              double alpha,
              //
              double a,
              double b,
              double betta );





extern double itemOfInteg_1SpecType(
              double Py,
              double Qy,
              //
              double Gx,
              double Hx,
              //
              double a,
              double b );





extern double integUnderRectAng_OneCell(
              double par_a,                           //   -  Solution parameter.
              //
              double lbDom,                           //   -  Left and right boundaries of rectangular domain.
              double rbDom,
              //
              double bbDom,                           //   -  Botton and upper boundaries of rectangular domain.
              double ubDom,
              //
              double Py,
              double Qy,
              //
              double Gx,
              double Hx,
              //
              double tau,
              int iCurrTL,                            //   -  Index of current time layer.
              //
              int * indCurSqOx,                       //   -  Index of current square by Ox axis.
              int * indCurSqOy,                       //   -  Index of current square by Oy axis.
              //
              double * masOX,                         //   -  Massive of OX steps. Dimension = numOfOXSt +1.
              int numOfOXSt,                          //   -  Number of OX steps.
              //
              double * masOY,                         //   -  Massive of OY steps. Dimension = numOfOYSt +1.
              int numOfOYSt,                          //   -  Number of OY steps.
              //
              double * rhoInPrevTL_asV );





extern double integUnderLeftTr_OneCell(
              double par_a,                           //   -  Solution parameter.
              //
              double lbDom,                           //   -  Left and right boundaries of rectangular domain.
              double rbDom,
              //
              double bbDom,                           //   -  Botton and upper boundaries of rectangular domain.
              double ubDom,
              //
              double Py,
              double Qy,
              //
              double a_SL,
              double b_SL,
              double Hx,
              //
              double tau,
              int iCurrTL,                            //   -  Index of current time layer.
              //
              int * indCurSqOx,                       //   -  Index of current square by Ox axis.
              int * indCurSqOy,                       //   -  Index of current square by Oy axis.
              //
              double * masOX,                         //   -  Massive of OX steps. Dimension = numOfOXSt +1.
              int numOfOXSt,                          //   -  Number of OX steps.
              //
              double * masOY,                         //   -  Massive of OY steps. Dimension = numOfOYSt +1.
              int numOfOYSt,                          //   -  Number of OY steps.
              //
              double * rhoInPrevTL_asV );





extern double integUnderRightTr_OneCell(
              double par_a,                           //   -  Solution parameter.
              //
              double lbDom,                           //   -  Left and right boundaries of rectangular domain.
              double rbDom,
              //
              double bbDom,                           //   -  Botton and upper boundaries of rectangular domain.
              double ubDom,
              //
              double Py,
              double Qy,
              //
              double a_SL,
              double b_SL,
              double Gx,
              //
              double tau,
              int iCurrTL,                            //   -  Index of current time layer.
              //
              int * indCurSqOx,                       //   -  Index of current square by Ox axis.
              int * indCurSqOy,                       //   -  Index of current square by Oy axis.
              //
              double * masOX,                         //   -  Massive of OX steps. Dimension = numOfOXSt +1.
              int numOfOXSt,                          //   -  Number of OX steps.
              //
              double * masOY,                         //   -  Massive of OY steps. Dimension = numOfOYSt +1.
              int numOfOYSt,                          //   -  Number of OY steps.
              //
              double * rhoInPrevTL_asV );


