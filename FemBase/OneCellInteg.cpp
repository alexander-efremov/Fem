


#include <math.h>

#include "InitAndBoundaryData.h"                      //   -  Initial and Boundary data.



double itemOfInteg_2SpecType(
              double Py,
              double Qy,
              //
              double alpha,
              //
              double a,
              double b,
              double betta )
{



double buf_D, integ;



//   Computing...

buf_D = (Qy - alpha) * (a*Qy + b - betta) * (a*Qy + b - betta) * (a*Qy + b - betta);

buf_D = buf_D  -  (Py - alpha) * (a*Py + b - betta) * (a*Py + b - betta) * (a*Py + b - betta);

integ = buf_D / (3. * a);



buf_D = (a*Qy + b - betta) * (a*Qy + b - betta) * (a*Qy + b - betta) * (a*Qy + b - betta);

buf_D = buf_D - (a*Py + b - betta) * (a*Py + b - betta) * (a*Py + b - betta) * (a*Py + b - betta);

integ = integ  -  buf_D / (12. *a *a);



return integ;
}







double itemOfInteg_1SpecType(
              double Py,
              double Qy,
              //
              double Gx,
              double Hx,
              //
              double a,
              double b )
{



double integ;



//   Computing...

integ = (Hx - a)*(Hx - a)  -  (Gx - a)*(Gx - a);

integ = integ * (  (Qy - b)*(Qy - b)  -  (Py - b)*(Py - b)  );

integ = integ / 4.;



return integ;
}







double integUnderRectAng_OneCell(
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
              double * rhoInPrevTL_asV )
{



//   return ( fabs( (Qy - Py) * (Hx - Gx) ) );



double hx = masOX[1] - masOX[0];

double hy = masOY[1] - masOY[0];

double integ = 0;

double buf_D;



double rho[2][2];

double t = tau * (iCurrTL -1.);

double x, y;



if(   (indCurSqOx[0] >=0) && (indCurSqOy[0] >=0)  )
{
   rho[0][0] = rhoInPrevTL_asV[ (numOfOXSt +1)*indCurSqOy[0] + indCurSqOx[0] ];

   rho[0][1] = rhoInPrevTL_asV[ (numOfOXSt +1)*indCurSqOy[1] + indCurSqOx[0] ];

   rho[1][0] = rhoInPrevTL_asV[ (numOfOXSt +1)*indCurSqOy[0] + indCurSqOx[1] ];

   rho[1][1] = rhoInPrevTL_asV[ (numOfOXSt +1)*indCurSqOy[1] + indCurSqOx[1] ];

   //   anSol  =  analytSolut( par_a,   lbDom, rbDom,   bbDom, ubDom,   iCurrTL *tau, masOX[ iOfOXN ], masOY[ iOfOYN ] );

   //   rhoInPrevTL_asV[ (numOfOXSt +1)*iOfOYN + iOfOXN ]  =  anSol  -  rhoInPrevTL_asV[ (numOfOXSt +1)*iOfOYN + iOfOXN ];

   //   rho[0][0] = rhoInPrevTL[ indCurSqOx[0] ][ indCurSqOy[0] ];

   //   rho[0][1] = rhoInPrevTL[ indCurSqOx[0] ][ indCurSqOy[1] ];

   //   rho[1][0] = rhoInPrevTL[ indCurSqOx[1] ][ indCurSqOy[0] ];

   //   rho[1][1] = rhoInPrevTL[ indCurSqOx[1] ][ indCurSqOy[1] ];
}
else
{
   x = indCurSqOx[0] * hx;   y = indCurSqOy[0] * hy;

   rho[0][0]  =  analytSolut( par_a,   lbDom, rbDom, bbDom, ubDom,   t, x, y );


   x = indCurSqOx[0] * hx;   y = indCurSqOy[1] * hy;

   rho[0][1]  =  analytSolut( par_a,   lbDom, rbDom, bbDom, ubDom,   t, x, y );


   x = indCurSqOx[1] * hx;   y = indCurSqOy[0] * hy;

   rho[1][0]  =  analytSolut( par_a,   lbDom, rbDom, bbDom, ubDom,   t, x, y );


   x = indCurSqOx[1] * hx;   y = indCurSqOy[1] * hy;

   rho[1][1]  =  analytSolut( par_a,   lbDom, rbDom, bbDom, ubDom,   t, x, y );
}




/*
rhoInPrevTL[ indCurSqOx[0] ][ indCurSqOy[0] ] = 1.;

rhoInPrevTL[ indCurSqOx[0] ][ indCurSqOy[1] ] = 1.;

rhoInPrevTL[ indCurSqOx[1] ][ indCurSqOy[0] ] = 1.;

rhoInPrevTL[ indCurSqOx[1] ][ indCurSqOy[1] ] = 1.;



Py = masOY[ indCurSqOy[0] ];

Qy = masOY[ indCurSqOy[1] ];

Gx = masOX[ indCurSqOx[0] ];

Hx = masOX[ indCurSqOx[1] ];
*/




if(   (indCurSqOx[1] >= 0) && (indCurSqOy[1] >= 0)   )
{
   buf_D = itemOfInteg_1SpecType( Py, Qy, Gx, Hx, masOX[ indCurSqOx[1] ], masOY[ indCurSqOy[1] ] );
}
else
{
   buf_D = itemOfInteg_1SpecType( Py, Qy, Gx, Hx,   hx *indCurSqOx[1]   , hy * indCurSqOy[1] );
}


buf_D = buf_D  /hx /hy;

integ = buf_D * rho[0][0];                            //   rhoInPrevTL[ indCurSqOx[0] ][ indCurSqOy[0] ];



if(   (indCurSqOx[0] >= 0)  &&   (indCurSqOy[1] >= 0)   )
{
   buf_D = itemOfInteg_1SpecType( Py, Qy, Gx, Hx, masOX[ indCurSqOx[0] ], masOY[ indCurSqOy[1] ] );
}
else
{
   buf_D = itemOfInteg_1SpecType( Py, Qy, Gx, Hx,   hx * indCurSqOx[0]  , hy * indCurSqOy[1] );
}


buf_D = buf_D  /hx /hy;

integ = integ - buf_D * rho[1][0];                    //   rhoInPrevTL[ indCurSqOx[1] ][ indCurSqOy[0] ];



if(   (indCurSqOx[1] >= 0)  &&  (indCurSqOy[0] >= 0)   )
{
   buf_D = itemOfInteg_1SpecType( Py, Qy, Gx, Hx, masOX[ indCurSqOx[1] ], masOY[ indCurSqOy[0] ] );
}
else
{
   buf_D = itemOfInteg_1SpecType( Py, Qy, Gx, Hx,   hx * indCurSqOx[1]  , hy * indCurSqOy[0] );
}


buf_D = buf_D  /hx /hy;

integ = integ - buf_D * rho[0][1];                    //   rhoInPrevTL[ indCurSqOx[0] ][ indCurSqOy[1] ];



if(   (indCurSqOx[0] >= 0)  &&  (indCurSqOy[0] >= 0)   )
{
   buf_D = itemOfInteg_1SpecType( Py, Qy, Gx, Hx, masOX[ indCurSqOx[0] ], masOY[ indCurSqOy[0] ] );
}
else
{
   buf_D = itemOfInteg_1SpecType( Py, Qy, Gx, Hx,   hx * indCurSqOx[0]  , hy * indCurSqOy[0] );
}


buf_D = buf_D  /hx /hy;

integ = integ + buf_D * rho[1][1];                    //   rhoInPrevTL[ indCurSqOx[1] ][ indCurSqOy[1] ];





return integ;
}







double integUnderLeftTr_OneCell(
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
              double * rhoInPrevTL_asV )
{



//   return ( fabs( (Qy-Py) * a_SL * (Qy-Py) /2.) );



double hx = masOX[1] - masOX[0];

double hy = masOY[1] - masOY[0];

double integ = 0;

double buf_D, bufInteg_D;



double rho[2][2];

double t = tau * (iCurrTL - 1.);

double x, y;



if(  (indCurSqOx[0] >=0)  &&  (indCurSqOx[1] <=numOfOXSt)  )
{
   if(  (indCurSqOy[0] >=0)  &&  (indCurSqOy[1] <=numOfOYSt)  )
   {
      rho[0][0] = rhoInPrevTL_asV[ (numOfOXSt +1)*indCurSqOy[0] + indCurSqOx[0] ];

      rho[0][1] = rhoInPrevTL_asV[ (numOfOXSt +1)*indCurSqOy[1] + indCurSqOx[0] ];

      rho[1][0] = rhoInPrevTL_asV[ (numOfOXSt +1)*indCurSqOy[0] + indCurSqOx[1] ];

      rho[1][1] = rhoInPrevTL_asV[ (numOfOXSt +1)*indCurSqOy[1] + indCurSqOx[1] ];

      //   anSol  =  analytSolut( par_a,   lbDom, rbDom,   bbDom, ubDom,   iCurrTL *tau, masOX[ iOfOXN ], masOY[ iOfOYN ] );

      //   rhoInPrevTL_asV[ (numOfOXSt +1)*iOfOYN + iOfOXN ]  =  anSol  -  rhoInPrevTL_asV[ (numOfOXSt +1)*iOfOYN + iOfOXN ];

      //   rho[0][0] = rhoInPrevTL[ indCurSqOx[0] ][ indCurSqOy[0] ];

      //   rho[0][1] = rhoInPrevTL[ indCurSqOx[0] ][ indCurSqOy[1] ];

      //   rho[1][0] = rhoInPrevTL[ indCurSqOx[1] ][ indCurSqOy[0] ];

      //   rho[1][1] = rhoInPrevTL[ indCurSqOx[1] ][ indCurSqOy[1] ];
   }
}



if(  (indCurSqOx[0] < 0)  ||  (indCurSqOx[1] > numOfOXSt)  ||  (indCurSqOy[0] < 0)  ||  (indCurSqOy[1] > numOfOYSt)  )
{
   x = indCurSqOx[0] * hx;   y = indCurSqOy[0] * hy;

   rho[0][0]  =  analytSolut( par_a,   lbDom, rbDom, bbDom, ubDom,   t, x, y );


   x = indCurSqOx[0] * hx;   y = indCurSqOy[1] * hy;

   rho[0][1]  =  analytSolut( par_a,   lbDom, rbDom, bbDom, ubDom,   t, x, y );


   x = indCurSqOx[1] * hx;   y = indCurSqOy[0] * hy;

   rho[1][0]  =  analytSolut( par_a,   lbDom, rbDom, bbDom, ubDom,   t, x, y );


   x = indCurSqOx[1] * hx;   y = indCurSqOy[1] * hy;

   rho[1][1]  =  analytSolut( par_a,   lbDom, rbDom, bbDom, ubDom,   t, x, y );
}




/*
Py = masOY[ indCurSqOy[0] ];

Qy = masOY[ indCurSqOy[1] ];


Hx = masOX[ indCurSqOx[1] ];

a_SL =-1.;

b_SL = Hx + Py;   //   Hx - Qy;   //   Hx + Py;



rhoInPrevTL[ indCurSqOx[0] ][ indCurSqOy[0] ] = 1.;

rhoInPrevTL[ indCurSqOx[0] ][ indCurSqOy[1] ] = 1.;

rhoInPrevTL[ indCurSqOx[1] ][ indCurSqOy[0] ] = 1.;

rhoInPrevTL[ indCurSqOx[1] ][ indCurSqOy[1] ] = 1.;
*/




//   1.

buf_D = (Qy - masOY[ indCurSqOy[1] ]) * (Qy - masOY[ indCurSqOy[1] ])  -  (Py - masOY[ indCurSqOy[1] ]) * (Py - masOY[ indCurSqOy[1] ]);


if(  (indCurSqOx[1] >= 0)  &&  (indCurSqOy[1] >= 0)  )
{
   buf_D = buf_D  *  (Hx - masOX[ indCurSqOx[1] ])  *  (Hx - masOX[ indCurSqOx[1] ]) /4.;

   bufInteg_D = itemOfInteg_2SpecType( Py, Qy, masOY[ indCurSqOy[1] ], a_SL, b_SL, masOX[ indCurSqOx[1] ] );
}
else
{
   buf_D = buf_D  *  (Hx -   hx * indCurSqOx[1]  )  *  (Hx -   hx * indCurSqOx[1]  ) /4.;

   bufInteg_D = itemOfInteg_2SpecType( Py, Qy, hy * indCurSqOy[1], a_SL, b_SL,    hx * indCurSqOx[1]  );
}


buf_D = buf_D  -  bufInteg_D /2.;

integ = buf_D * rho[0][0] /hx /hy;



//   2.

buf_D = (Qy - masOY[ indCurSqOy[1] ]) * (Qy - masOY[ indCurSqOy[1] ])  -  (Py - masOY[ indCurSqOy[1] ]) * (Py - masOY[ indCurSqOy[1] ]);


if(  (indCurSqOx[0] >= 0)  &&  (indCurSqOy[1] >= 0)  )
{
   buf_D = -1. * buf_D  *  (Hx - masOX[ indCurSqOx[0] ])  *  (Hx - masOX[ indCurSqOx[0] ]) /4.;

   bufInteg_D = itemOfInteg_2SpecType( Py, Qy, masOY[ indCurSqOy[1] ], a_SL, b_SL, masOX[ indCurSqOx[0] ] );
}
else
{
   buf_D = -1. * buf_D  *  (Hx -   hx * indCurSqOx[0]  )  *  (Hx -   hx * indCurSqOx[0]  ) /4.;

   bufInteg_D = itemOfInteg_2SpecType( Py, Qy, hy * indCurSqOy[1], a_SL, b_SL,   hx * indCurSqOx[0]   );
}


buf_D = buf_D  +  bufInteg_D /2.;

integ = integ  +  buf_D * rho[1][0] /hx /hy;



//   3.

buf_D = (Qy - masOY[ indCurSqOy[0] ]) * (Qy - masOY[ indCurSqOy[0] ])  -  (Py - masOY[ indCurSqOy[0] ]) * (Py - masOY[ indCurSqOy[0] ]);


if(  (indCurSqOx[1] >= 0)  &&  (indCurSqOy[0] >= 0)  )
{
   buf_D = -1. * buf_D  *  (Hx - masOX[ indCurSqOx[1] ])  *  (Hx - masOX[ indCurSqOx[1] ]) /4.;

   bufInteg_D = itemOfInteg_2SpecType( Py, Qy, masOY[ indCurSqOy[0] ], a_SL, b_SL, masOX[ indCurSqOx[1] ] );
}
else
{
   buf_D = -1. * buf_D  *  (Hx -   hx * indCurSqOx[1]  )  *  (Hx -   hx * indCurSqOx[1]  ) /4.;

   bufInteg_D = itemOfInteg_2SpecType( Py, Qy, hy * indCurSqOy[0], a_SL, b_SL,   hx * indCurSqOx[1]   );
}


buf_D = buf_D  +  bufInteg_D /2.;

integ = integ  +  buf_D * rho[0][1] /hx /hy;



//   4.

buf_D = (Qy - masOY[ indCurSqOy[0] ]) * (Qy - masOY[ indCurSqOy[0] ])  -  (Py - masOY[ indCurSqOy[0] ]) * (Py - masOY[ indCurSqOy[0] ]);


if(  (indCurSqOx[0] >= 0)  &&  (indCurSqOy[0] >= 0)  )
{
   buf_D = buf_D  *  (Hx - masOX[ indCurSqOx[0] ])  *  (Hx - masOX[ indCurSqOx[0] ]) /4.;

   bufInteg_D = itemOfInteg_2SpecType( Py, Qy, masOY[ indCurSqOy[0] ], a_SL, b_SL, masOX[ indCurSqOx[0] ] );
}
else
{
   buf_D = buf_D  *  (Hx -   hx * indCurSqOx[0]  )  *  (Hx -   hx * indCurSqOx[0]  ) /4.;

   bufInteg_D = itemOfInteg_2SpecType( Py, Qy, hy * indCurSqOy[0], a_SL, b_SL,   hx * indCurSqOx[0]   );
}


buf_D = buf_D  -  bufInteg_D /2.;

integ = integ  +  buf_D * rho[1][1] /hx /hy;





return integ;
}







double integUnderRightTr_OneCell(
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
              double * rhoInPrevTL_asV )
{



double bufInteg_D, integ;



bufInteg_D = integUnderLeftTr_OneCell(
              par_a,                                  //   -  Solution parameter.
              //
              lbDom, rbDom,                           //   -  Left and right boundaries of rectangular domain.
              //
              bbDom, ubDom,                           //   -  Botton and upper boundaries of rectangular domain.
              //
              Py, Qy,
              //
              a_SL, b_SL,
              Gx,                                     //   -  double Hx,
              //
              tau, iCurrTL,                           //   -  Index of current time layer.
              //
              indCurSqOx,                             //   -  Index of current square by Ox axis.
              indCurSqOy,                             //   -  Index of current square by Oy axis.
              //
              masOX, numOfOXSt,                       //   -  Massive of OX steps. Dimension = numOfOXSt +1. Number of OX steps.
              //
              masOY, numOfOYSt,                       //   -  Massive of OY steps. Dimension = numOfOYSt +1. Number of OY steps.
              //
              rhoInPrevTL_asV );



integ = -1. * bufInteg_D;





return integ;
}
