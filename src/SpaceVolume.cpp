


#include <iostream>

using namespace std;

#include <math.h>

#include "InitAndBoundaryData.h"                      //   -  Initial and Boundary data.

#include "integUndRigAngTr.h"





double integUnderBottTr(
              double par_a,                           //   -  Item of left and right setback (parameter "a" in test).
              double par_b,                           //   -  Item of second parameter from "u_funcion".
              //
              double lbDom,                           //   -  Left and right boundaries of rectangular domain.
              double rbDom,
              //
              double bbDom,                           //   -  Botton and upper boundaries of rectangular domain.
              double ubDom,
              //
              double tau,
              int iCurrTL,                            //   -  Index of current time layer.
              //
              double * LvBt,                          //   -  Left, Right and Botton vertices of Botton triangle.
              double * RvBt,                          //   -  Left, Right and Botton vertices of Botton triangle.
              double * BvBt,                          //   -  Left, Right and Botton vertices of Botton triangle.
              //
              double * masOX,                         //   -  Massive of OX steps. Dimension = numOfOXSt +1.
              int numOfOXSt,                          //   -  Number of OX steps.
              //
              double * masOY,                         //   -  Massive of OY steps. Dimension = numOfOYSt +1.
              int numOfOYSt,                          //   -  Number of OY steps.
              //
              double * rhoInPrevTL_asV )
{



double integOfBottTr;

double buf_D;



//   Three ways are possible.

//   1.

if(  BvBt[0] <= LvBt[0]  )
{

   buf_D = integUnderRigAngTr_BottRight(
              par_a,   lbDom, rbDom,   bbDom, ubDom,   tau, iCurrTL,
              //
              BvBt, RvBt,   masOX, numOfOXSt, masOY, numOfOYSt,   rhoInPrevTL_asV );

   integOfBottTr = buf_D;


   buf_D = integUnderRigAngTr_BottRight(
              par_a,   lbDom, rbDom,   bbDom, ubDom,   tau, iCurrTL,
              //
              BvBt, LvBt,   masOX, numOfOXSt, masOY, numOfOYSt,   rhoInPrevTL_asV );

   integOfBottTr = integOfBottTr - buf_D;


   return integOfBottTr;
}



//   2.

if(  (BvBt[0] > LvBt[0]) && (BvBt[0] < RvBt[0]) )
{
/*
   BvBt[0] = masOX[ 40 ];   BvBt[1] = masOY[ 43 ];

   LvBt[0] = masOX[ 37 ];   LvBt[1] = masOY[ 45 ] + 0.01;
*/
   buf_D = integUnderRigAngTr_BottLeft(
              par_a,   lbDom, rbDom,   bbDom, ubDom,   tau, iCurrTL,
              //
              BvBt, LvBt,   masOX, numOfOXSt, masOY, numOfOYSt,   rhoInPrevTL_asV );

   integOfBottTr = buf_D;


/*
   BvBt[0] = 0.25;       BvBt[1] = 0.1;

   RvBt[0] = 0.41;       RvBt[1] = 0.6;

   buf_D = integUnderRigAngTr_BottRight(
              par_a,   lbDom, rbDom,   bbDom, ubDom,   tau, iCurrTL,
              //
              BvBt, RvBt,   masOX, numOfOXSt, masOY, numOfOYSt,   rhoInPrevTL );


   BvBt[0] =-0.25;       BvBt[1] = 0.1;

   LvBt[0] =-0.41;       LvBt[1] = 0.6;

   RvBt[0] = integUnderRigAngTr_BottLeft(
              par_a,   lbDom, rbDom,   bbDom, ubDom,   tau, iCurrTL,
              //
              BvBt, LvBt,   masOX, numOfOXSt, masOY, numOfOYSt,   rhoInPrevTL );
*/


   buf_D = integUnderRigAngTr_BottRight(
              par_a,   lbDom, rbDom,   bbDom, ubDom,   tau, iCurrTL,
              //
              BvBt, RvBt,   masOX, numOfOXSt, masOY, numOfOYSt,   rhoInPrevTL_asV );

   integOfBottTr = integOfBottTr + buf_D;





   return integOfBottTr;
}



//   3.

if(  BvBt[0] >= RvBt[0]  )
{
/*
   BvBt[0] = masOX[ 20 ];   BvBt[1] = masOY[ 43 ];

   LvBt[0] = masOX[ 18 ];   LvBt[1] = masOY[ 44 ];
*/
   buf_D = integUnderRigAngTr_BottLeft(
              par_a,   lbDom, rbDom,   bbDom, ubDom,   tau, iCurrTL,
              //
              BvBt, LvBt,   masOX, numOfOXSt, masOY, numOfOYSt,   rhoInPrevTL_asV );

   integOfBottTr = buf_D;


   buf_D = integUnderRigAngTr_BottLeft(
              par_a,   lbDom, rbDom,   bbDom, ubDom,   tau, iCurrTL,
              //
              BvBt, RvBt,   masOX, numOfOXSt, masOY, numOfOYSt,   rhoInPrevTL_asV );

   integOfBottTr = integOfBottTr - buf_D;


   return integOfBottTr;
}





return integOfBottTr;
}







double integUnderUpperTr(
              double par_a,                           //   -  Item of left and right setback (parameter "a" in test).
              double par_b,                           //   -  Item of second parameter from "u_funcion".
              //
              double lbDom,                           //   -  Left and right boundaries of rectangular domain.
              double rbDom,
              //
              double bbDom,                           //   -  Botton and upper boundaries of rectangular domain.
              double ubDom,
              //
              double tau,
              int iCurrTL,                            //   -  Index of current time layer.
              //
              double * LvUt,                          //   -  Left, Right and Upper vertices of Upper triangle.
              double * RvUt,                          //   -  Left, Right and Upper vertices of Upper triangle.
              double * UvUt,                          //   -  Left, Right and Upper vertices of Upper triangle.
              //
              double * masOX,                         //   -  Massive of OX steps. Dimension = numOfOXSt +1.
              int numOfOXSt,                          //   -  Number of OX steps.
              //
              double * masOY,                         //   -  Massive of OY steps. Dimension = numOfOYSt +1.
              int numOfOYSt,                          //   -  Number of OY steps.
              //
              double * rhoInPrevTL_asV )
{



double integOfUppTr;

double buf_D;



//   Three ways are possible.

//   1.

if(  UvUt[0] <= LvUt[0]  )
{
   buf_D = integUnderRigAngTr_UppRight(
              par_a,   lbDom, rbDom,   bbDom, ubDom,   tau, iCurrTL,
              //
              RvUt, UvUt,   masOX, numOfOXSt, masOY, numOfOYSt,   rhoInPrevTL_asV );

   integOfUppTr = buf_D;


   buf_D = integUnderRigAngTr_UppRight(
              par_a,   lbDom, rbDom,   bbDom, ubDom,   tau, iCurrTL,
              //
              LvUt, UvUt,   masOX, numOfOXSt, masOY, numOfOYSt,   rhoInPrevTL_asV );

   integOfUppTr = integOfUppTr - buf_D;


   return integOfUppTr;
}



//   2.

if(  (UvUt[0] > LvUt[0]) && (UvUt[0] < RvUt[0]) )
{
   buf_D = integUnderRigAngTr_UppLeft(
              par_a,   lbDom, rbDom,   bbDom, ubDom,   tau, iCurrTL,
              //
              LvUt, UvUt,   masOX, numOfOXSt, masOY, numOfOYSt,   rhoInPrevTL_asV );

   integOfUppTr = buf_D;






/*
   RvUt[0] = 0.05;         RvUt[1] = 0.;

   UvUt[0] = 0.025;        UvUt[1] = 0.05;

   buf_D = integUnderRigAngTr_UppRight(
              par_a,   lbDom, rbDom,   bbDom, ubDom,   tau, iCurrTL,
              //
              RvUt, UvUt,   masOX, numOfOXSt, masOY, numOfOYSt,   rhoInPrevTL );


   LvUt[0] =-0.05;         LvUt[1] = 0.;

   UvUt[0] =-0.025;        UvUt[1] = 0.05;

   RvUt[0] = integUnderRigAngTr_UppLeft(
              par_a,   lbDom, rbDom,   bbDom, ubDom,   tau, iCurrTL,
              //
              LvUt, UvUt,   masOX, numOfOXSt, masOY, numOfOYSt,   rhoInPrevTL );
*/







   buf_D = integUnderRigAngTr_UppRight(
              par_a,   lbDom, rbDom,   bbDom, ubDom,   tau, iCurrTL,
              //
              RvUt, UvUt,   masOX, numOfOXSt, masOY, numOfOYSt,   rhoInPrevTL_asV );

   integOfUppTr = integOfUppTr + buf_D;


   return integOfUppTr;
}



//   3.

if(  UvUt[0] >= RvUt[0]  )
{
   buf_D = integUnderRigAngTr_UppLeft(
              par_a,   lbDom, rbDom,   bbDom, ubDom,   tau, iCurrTL,
              //
              LvUt, UvUt,   masOX, numOfOXSt, masOY, numOfOYSt,   rhoInPrevTL_asV );

   integOfUppTr = buf_D;


   buf_D = integUnderRigAngTr_UppLeft(
              par_a,   lbDom, rbDom,   bbDom, ubDom,   tau, iCurrTL,
              //
              RvUt, UvUt,   masOX, numOfOXSt, masOY, numOfOYSt,   rhoInPrevTL_asV );

   integOfUppTr = integOfUppTr - buf_D;


   return integOfUppTr;
}





return integOfUppTr;
}







double integUnderUnunifTr(
              double par_a,                           //   -  Item of left and right setback (parameter "a" in test).
              double par_b,                           //   -  Item of second parameter from "u_funcion".
              //
              double lbDom,                           //   -  Left and right boundaries of rectangular domain.
              double rbDom,
              //
              double bbDom,                           //   -  Botton and upper boundaries of rectangular domain.
              double ubDom,
              //
              double tau,
              int iCurrTL,                            //   -  Index of current time layer.
              //
              double * firVer,                        //   -  First vertex of triangle.
              double * secVer,                        //   -  Second vertex of triangle.
              double * thiVer,                        //   -  Third vertex of triangle.
              //
              double * masOX,                         //   -  Massive of OX steps. Dimension = numOfOXSt +1.
              int numOfOXSt,                          //   -  Number of OX steps.
              //
              double * masOY,                         //   -  Massive of OY steps. Dimension = numOfOYSt +1.
              int numOfOYSt,                          //   -  Number of OY steps.
              //
              double * rhoInPrevTL_asV )
{



double bv[2], mv[2], uv[2];                           //   -  Botton, middle and upper vertices of triangle.


bool isFirVUsed = false;

bool isSecVUsed = false;

bool isThiVUsed = false;


bool is1VUsed, is2VUsed, is3VUsed;


double a_LC, b_LC, c_LC;                              //   -  Coefficients of line betweeen "bv" and "uv" vertices.

double ap[2];                                         //   -  Across point of line through "bv" to "uv" and "y == mv[1]"


double LvBt[2], RvBt[2], BvBt[2];                     //   -  Left, Right and Botton vertices of Botton triangle.

double integOfBottTr;                                 //   -  Item of integral under Botton triangle.

double LvUt[2], RvUt[2], UvUt[2];                     //   -  Left, Right and Upper vertices of Upper triangle.

double integOfUppTr;                                  //   -  Item of integral under Upper triangle.


double integ = 0.;                                    //   -  Item which I'm computing.





//   1. I need to understand which vertex is botton, middle and upper.

bv[1] = firVer[1];   bv[0] = firVer[0];   isFirVUsed = true;

if( bv[1] > secVer[1] ){ bv[1] = secVer[1];   bv[0] = secVer[0];   isFirVUsed = false;   isSecVUsed = true; }

if( bv[1] > thiVer[1] ){ bv[1] = thiVer[1];   bv[0] = thiVer[0];   isFirVUsed = false;   isSecVUsed = false;  isThiVUsed = true; }



uv[1] = masOY[0];                                     //   -  The minimum possible value.

is1VUsed = false;   is2VUsed = false;   is3VUsed = false;

if(  (uv[1] < firVer[1])  &&  (isFirVUsed == false)  ){  uv[1] = firVer[1];   uv[0] = firVer[0];   is1VUsed = true; }

if(  (uv[1] < secVer[1])  &&  (isSecVUsed == false)  ){  uv[1] = secVer[1];   uv[0] = secVer[0];   is2VUsed = true;   is1VUsed = false; }

if(  (uv[1] < thiVer[1])  &&  (isThiVUsed == false)  ){  uv[1] = thiVer[1];   uv[0] = thiVer[0];   is3VUsed = true;   is2VUsed = false;   is1VUsed = false; }



//   Dangerous.

if(  (isFirVUsed == false) &&  (is1VUsed == false)  ){  mv[1] = firVer[1];   mv[0] = firVer[0];  }

if(  (isSecVUsed == false) &&  (is2VUsed == false)  ){  mv[1] = secVer[1];   mv[0] = secVer[0];  }

if(  (isThiVUsed == false) &&  (is3VUsed == false)  ){  mv[1] = thiVer[1];   mv[0] = thiVer[0];  }



//   2. I want to compute across point.

//   2.a Let's compute line coefficients betweeen "bv" and "uv" vertices.
//   a_LC * x  +  b_LC * y  = c_LC.

a_LC  =  uv[1] - bv[1];

b_LC  =  bv[0] - uv[0];

c_LC  =  (bv[0] - uv[0])*bv[1]  +  (uv[1] - bv[1])*bv[0];



//   2.b Across point.

ap[1]  =  mv[1];

if( fabs(a_LC) < 1.e-12 )
{
   //   This triangle has very small height. I guess further computation isn't correct.

   return 1.e-12;
}

ap[0]  =  (c_LC  -  b_LC * ap[1])  /a_LC;



//   3. There the middle vertex relativly straight line is? Two ways are possible.

if( mv[0] < ap[0] )
{
   //   Left, Right and Botton vertices of Botton triangle.

   LvBt[0]  =  mv[0];   LvBt[1]  =  mv[1];

   RvBt[0]  =  ap[0];   RvBt[1]  =  ap[1];

   BvBt[0]  =  bv[0];   BvBt[1]  =  bv[1];

   integOfBottTr = integUnderBottTr(
              par_a, par_b,
              //
              lbDom, rbDom,                           //   -  Left and right boundaries of rectangular domain.
              //
              bbDom, ubDom,                           //   -  Botton and upper boundaries of rectangular domain.
              //
              tau, iCurrTL,                           //   -  Index of current time layer.
              //
              LvBt, RvBt, BvBt,                       //   -  Left, Right and Botton vertices of Botton triangle.
              //
              masOX, numOfOXSt,                       //   -  Number of OX steps.
              //
              masOY, numOfOYSt,                       //   -  Number of OY steps.
              //
              rhoInPrevTL_asV );

   integ = integOfBottTr;


   //   Left, Right and Upper vertices of Upper triangle.

   LvUt[0]  =  mv[0];   LvUt[1]  =  mv[1];

   RvUt[0]  =  ap[0];   RvUt[1]  =  ap[1];

   UvUt[0]  =  uv[0];   UvUt[1]  =  uv[1];

   integOfUppTr = integUnderUpperTr(
              par_a, par_b,
              //
              lbDom, rbDom,                           //   -  Left and right boundaries of rectangular domain.
              //
              bbDom, ubDom,                           //   -  Botton and upper boundaries of rectangular domain.
              //
              tau, iCurrTL,                           //   -  Index of current time layer.
              //
              LvUt, RvUt, UvUt,                       //   -  Left, Right and Botton vertices of Upper triangle.
              //
              masOX, numOfOXSt,                       //   -  Number of OX steps.
              //
              masOY, numOfOYSt,                       //   -  Number of OY steps.
              //
              rhoInPrevTL_asV );

   integ = integ + integOfUppTr;

   return integ;
}



if( mv[0] >= ap[0] )
{
   //   Left, Right and Botton vertices of Botton triangle.

   LvBt[0]  =  ap[0];   LvBt[1]  =  ap[1];

   RvBt[0]  =  mv[0];   RvBt[1]  =  mv[1];

   BvBt[0]  =  bv[0];   BvBt[1]  =  bv[1];

   integOfBottTr = integUnderBottTr(
              par_a, par_b,
              //
              lbDom, rbDom,                           //   -  Left and right boundaries of rectangular domain.
              //
              bbDom, ubDom,                           //   -  Botton and upper boundaries of rectangular domain.
              //
              tau, iCurrTL,                           //   -  Index of current time layer.
              //
              LvBt, RvBt, BvBt,                       //   -  Left, Right and Botton vertices of Botton triangle.
              //
              masOX, numOfOXSt,                       //   -  Number of OX steps.
              //
              masOY, numOfOYSt,                       //   -  Number of OY steps.
              //
              rhoInPrevTL_asV );

   integ = integOfBottTr;


   //   Left, Right and Upper vertices of Upper triangle.

   LvUt[0]  =  ap[0];   LvUt[1]  =  ap[1];

   RvUt[0]  =  mv[0];   RvUt[1]  =  mv[1];

   UvUt[0]  =  uv[0];   UvUt[1]  =  uv[1];

   integOfUppTr = integUnderUpperTr(
              par_a, par_b,
              //
              lbDom, rbDom,                           //   -  Left and right boundaries of rectangular domain.
              //
              bbDom, ubDom,                           //   -  Botton and upper boundaries of rectangular domain.
              //
              tau, iCurrTL,                           //   -  Index of current time layer.
              //
              LvUt, RvUt, UvUt,                       //   -  Left, Right and Botton vertices of Upper triangle.
              //
              masOX, numOfOXSt,                       //   -  Number of OX steps.
              //
              masOY, numOfOYSt,                       //   -  Number of OY steps.
              //
              rhoInPrevTL_asV );

   integ = integ + integOfUppTr;

   return integ;
}



return integ;
}







int quadrAngleType(
              double par_a,                           //   -  Item of left and right setback (parameter "a" in test).
              double par_b,                           //   -  Item of second parameter from "u_funcion".
              //
              double lbDom,                           //   -  Left and right boundaries of rectangular domain.
              double rbDom,
              //
              double bbDom,                           //   -  Botton and upper boundaries of rectangular domain.
              double ubDom,
              //
              double tau,
              double iCurrTL,                         //   -  Index of current time layer. Necessary for velocity.
              //
              int iOfOXN,                             //   -  Index of current OX node.
              double * masOX,                         //   -  Massive of OX steps. Dimension = numOfOXSt +1.
              int numOfOXSt,                          //   -  Number of OX steps.
              //
              int iOfOYN,                             //   -  Index of current OY node.
              double * masOY,                         //   -  Massive of OY steps. Dimension = numOfOYSt +1.
              int numOfOYSt,                          //   -  Number of OY steps.
              //
              double * firVfirT,                      //   -  First vertex of first triangle.
              double * secVfirT,                      //   -  Second vertex of first triangle.
              double * thiVfirT,                      //   -  Third vertex of first triangle.
              //
              double * firVsecT,                      //   -  First vertex of second triangle.
              double * secVsecT,                      //   -  Second vertex of second triangle.
              double * thiVsecT )                     //   -  Third vertex of second triangle.

{



int qAngType = 0;   // Be default.                    //   -   if "qAngType == 0"  pseudo case. COMPUTATION IS IMPOSSIBLE.
                                                      //   -   if "qAngType == 1"  the quadrangle is convex. Computation possible.
                                                      //   -   if "qAngType == 2"  the quadrangle is concave. Computation possible.

double alpha[2], betta[2], gamma[2], theta[2];        //   -  Vertexes of square. Anticlockwise order from left botton vertex.

double u, v;                                          //   -  Items of velocity components.

double alNew[2], beNew[2], gaNew[2], thNew[2];        //   -  New positions of vertexes. Vertexes of quadrangle.



double vectAlGa[2], vectBeTh[2];                      //   -  Vectors: 1) from "alNew" to "gaNew" 2) from "beNew" to "thNew".

double a_1LC, b_1LC, c_1LC;                           //   -  a_1LC * x  +  b_1LC * y  = c_1LC. Line equation through "alNew" and "gaNew".

double a_2LC, b_2LC, c_2LC;                           //   -  a_2LC * x  +  b_2LC * y  = c_2LC. Line equation through "beNew" and "thNew".

double AcrP[2];                                       //   -  Across point of two lines.



double vectAlBe[2];                                   //   -  Vectors for computing vertexes sequence order by vector production.

double vectAlTh[2];                                   //   -  Vectors for computing vertexes sequence order by vector production.

double vectBeGa[2];                                   //   -  Vectors for computing vertexes sequence order by vector production.

double vectBeAl[2];                                   //   -  Vectors for computing vertexes sequence order by vector production.


double vectProdOz;                                    //   -  Z-axis of vector production.

double scalProd;                                      //   -  Scalar production of two vectors.



bool bul;





//   1. First of all let's compute coordinates of square vertexes.

//  OX:

if( iOfOXN == 0 )
{
   alpha[0]  =  masOX[ iOfOXN ];

   betta[0]  =  ( masOX[iOfOXN]  +  masOX[iOfOXN +1] ) /2.;

   gamma[0]  =  ( masOX[iOfOXN]  +  masOX[iOfOXN +1] ) /2.;

   theta[0]  =  masOX[ iOfOXN ];
}


if( iOfOXN == numOfOXSt )
{
   alpha[0]  =  ( masOX[iOfOXN -1]  +  masOX[iOfOXN] ) /2.;

   betta[0]  =  masOX[ iOfOXN ];

   gamma[0]  =  masOX[ iOfOXN ];

   theta[0]  =  ( masOX[iOfOXN -1]  +  masOX[iOfOXN] ) /2.;
}


if( (iOfOXN > 0)  &&  (iOfOXN < numOfOXSt) )
{
   alpha[0]  =  ( masOX[iOfOXN -1]  +  masOX[iOfOXN] ) /2.;

   betta[0]  =  ( masOX[iOfOXN +1]  +  masOX[iOfOXN] ) /2.;

   gamma[0]  =  ( masOX[iOfOXN +1]  +  masOX[iOfOXN] ) /2.;

   theta[0]  =  ( masOX[iOfOXN -1]  +  masOX[iOfOXN] ) /2.;
}



//  OY:

if( iOfOYN == 0 )
{
   alpha[1]  =  masOY[ iOfOYN ];

   betta[1]  =  masOY[ iOfOYN ];

   gamma[1]  =  ( masOY[iOfOYN]  +  masOY[ iOfOYN +1] ) /2.;

   theta[1]  =  ( masOY[iOfOYN]  +  masOY[ iOfOYN +1] ) /2.;
}


if( iOfOYN == numOfOYSt )
{
   alpha[1]  =  ( masOY[iOfOYN]  +  masOY[ iOfOYN -1] ) /2.;

   betta[1]  =  ( masOY[iOfOYN]  +  masOY[ iOfOYN -1] ) /2.;

   gamma[1]  =  masOY[ iOfOYN ];

   theta[1]  =  masOY[ iOfOYN ];
}


if( (iOfOYN > 0) && (iOfOYN < numOfOYSt) )
{
   alpha[1]  =  ( masOY[iOfOYN]  +  masOY[ iOfOYN -1] ) /2.;

   betta[1]  =  ( masOY[iOfOYN]  +  masOY[ iOfOYN -1] ) /2.;

   gamma[1]  =  ( masOY[iOfOYN]  +  masOY[ iOfOYN +1] ) /2.;

   theta[1]  =  ( masOY[iOfOYN]  +  masOY[ iOfOYN +1] ) /2.;
}



//   2. Now let's compute new coordinates on the previous time level of alpha, betta, gamma, theta points.

//  alNew.

u = u_function( par_b,   lbDom, rbDom,   bbDom, ubDom,   iCurrTL * tau, alpha[0], alpha[1] );

v = v_function( par_b,   lbDom, rbDom,   bbDom, ubDom,   iCurrTL * tau, alpha[0], alpha[1] );


alNew[0]  =  alpha[0]  -  tau * u;

alNew[1]  =  alpha[1]  -  tau * v;


//  beNew.

u = u_function( par_b,   lbDom, rbDom,   bbDom, ubDom,   iCurrTL * tau, betta[0], betta[1] );

v = v_function( par_b,   lbDom, rbDom,   bbDom, ubDom,   iCurrTL * tau, betta[0], betta[1] );

beNew[0]  =  betta[0]  -  tau * u;

beNew[1]  =  betta[1]  -  tau * v;


//  gaNew.

u = u_function( par_b,   lbDom, rbDom,   bbDom, ubDom,   iCurrTL * tau, gamma[0], gamma[1] );

v = v_function( par_b,   lbDom, rbDom,   bbDom, ubDom,   iCurrTL * tau, gamma[0], gamma[1] );

gaNew[0]  =  gamma[0]  -  tau * u;

gaNew[1]  =  gamma[1]  -  tau * v;


//  thNew.

u = u_function( par_b,   lbDom, rbDom,   bbDom, ubDom,   iCurrTL * tau, theta[0], theta[1] );

v = v_function( par_b,   lbDom, rbDom,   bbDom, ubDom,   iCurrTL * tau, theta[0], theta[1] );

thNew[0]  =  theta[0]  -  tau * u;

thNew[1]  =  theta[1]  -  tau * v;




/*
alNew[0] =-0.3;    alNew[1] = -0.1;

beNew[0] = 0.;     beNew[1] = 0.;

gaNew[0] = 0.1;    gaNew[1] = 0.1;

thNew[0] =-0.2;    thNew[1] = 0.1;
*/




//   3.a Let's compute coefficients of first line betweeen "alNew" and "gaNew" points.
//   a_1LC * x  +  b_1LC * y  = c_1LC.

vectAlGa[0] = gaNew[0] - alNew[0];

vectAlGa[1] = gaNew[1] - alNew[1];


a_1LC  =  vectAlGa[1];

b_1LC  = -vectAlGa[0];

c_1LC  =  vectAlGa[1] * alNew[0]  -  vectAlGa[0] * alNew[1];



//   3.b Let's compute coefficients of second line betweeen "beNew" and "thNew" points.
//   a_2LC * x  +  b_2LC * y  = c_2LC.

vectBeTh[0] = thNew[0] - beNew[0];

vectBeTh[1] = thNew[1] - beNew[1];


a_2LC  =  vectBeTh[1];

b_2LC  = -vectBeTh[0];

c_2LC  =  vectBeTh[1] * beNew[0]  -  vectBeTh[0] * beNew[1];



//   4. Let's compute coordinates of across point of this two lines.

//   Are lines parallel?

if( fabs(b_1LC*a_2LC - b_2LC*a_1LC)  <  1.e-14 )
{
   //   Not checked.

   qAngType = 0.;

   //   Pseudo case. Anyway I need to compute some values.

   //   First triangle.

   firVfirT[0] = alNew[0];   firVfirT[1] = alNew[1];

   secVfirT[0] = beNew[0];   secVfirT[1] = beNew[1];

   thiVfirT[0] = gaNew[0];   thiVfirT[1] = gaNew[1];


   //   Vertices of second triagle depends on scalar production.

   vectAlGa[0] = gaNew[0] - alNew[0];

   vectAlGa[1] = gaNew[1] - alNew[1];


   vectBeTh[0] = thNew[0] - beNew[0];

   vectBeTh[1] = thNew[1] - beNew[1];


   scalProd = vectAlGa[0] * vectBeTh[0]  +  vectAlGa[1] * vectBeTh[1];


   firVsecT[0] = beNew[0];   firVsecT[1] = beNew[1];

   secVsecT[0] = thNew[0];   secVsecT[1] = thNew[1];


   if( scalProd >= 0. )
   {
      thiVsecT[0] = gaNew[0];   thiVsecT[1] = gaNew[1];
   }


   if( scalProd < 0. )
   {
      thiVsecT[0] = alNew[0];   thiVsecT[1] = alNew[1];
   }

   return qAngType;
}


AcrP[0]  =  ( b_1LC*c_2LC - b_2LC*c_1LC ) / ( b_1LC*a_2LC - b_2LC*a_1LC );

AcrP[1]  =  ( a_1LC*c_2LC - a_2LC*c_1LC ) / (-b_1LC*a_2LC + b_2LC*a_1LC );



//   5. Checking "beNew" and "thNew" points position. It's a main criterion.

//   There are two ways:

//  5.a "(  (beNew[1] - AcrP[1])*(thNew[1] - AcrP[1])  )  > 0."

//  5.b "(  (beNew[1] - AcrP[1])*(thNew[1] - AcrP[1])  )  <=0."

//  Now let's consider first case 5.a "(  (beNew[1] - AcrP[1])*(thNew[1] - AcrP[1])  )  > 0."

if (  (  (beNew[1] - AcrP[1])*(thNew[1] - AcrP[1])  )  >  0.  )
{
   //   Vertices "beNew" and "thNew" are on the same side of across point by Oy-axis.

   //   Second criterion. Is across point between "alNew" and "gaNew" vertices by Ox-axis?

   //   Second criterion. First case.

   if(  (  (alNew[0] - AcrP[0])*(gaNew[0] - AcrP[0])  )  >  0.  )
   {
      //   The across point is NOT between "alNew" and "gaNew" vertices by Ox-axis!

      qAngType = 0;

      //   Pseudo case.  Anyway I need to find some solution. So

      firVfirT[0] = alNew[0];   firVfirT[1] = alNew[1];

      secVfirT[0] = beNew[0];   secVfirT[1] = beNew[1];

      thiVfirT[0] = gaNew[0];   thiVfirT[1] = gaNew[1];

      //   Second triangle.

      firVsecT[0] = beNew[0];   firVsecT[1] = beNew[1];

      secVsecT[0] = thNew[0];   secVsecT[1] = thNew[1];


      //   Third vertex computing...

      vectAlGa[0] = gaNew[0] - alNew[0];

      vectAlGa[1] = gaNew[1] - alNew[1];


      vectBeTh[0] = thNew[0] - beNew[0];

      vectBeTh[1] = thNew[1] - beNew[1];


      scalProd = vectAlGa[0] * vectBeTh[0]  +  vectAlGa[1] * vectBeTh[1];

      if( scalProd >= 0. )
      {
         thiVsecT[0] = gaNew[0];   thiVsecT[1] = gaNew[1];
      }


      if( scalProd < 0. )
      {
         thiVsecT[0] = alNew[0];   thiVsecT[1] = alNew[1];
      }

      return qAngType;

   }   //   "if(  (  (alNew[0] - AcrP[0])*(gaNew[0] - AcsP[0])  )  >  0.  )".


   //   Second criterion. Second case.

   if(  (  (alNew[0] - AcrP[0])*(gaNew[0] - AcrP[0])  )  <=  0.  )
   {
      //   The across point IS BETWEEN "alNew" and "gaNew" vertices by Ox-axis! It's good.

      //   THIRD criterion: is vertex "beNew" in triangle "alNew - gaNew - thNew".

      //   This is the same criterion for question: is vertex "thNew" in triangle "alNew - gaNew - beNew".

      //   We can check it by vector production. This is a second criterion.

      vectAlBe[0] = beNew[0] - alNew[0];

      vectAlBe[1] = beNew[1] - alNew[1];


      vectAlTh[0] = thNew[0] - alNew[0];

      vectAlTh[1] = thNew[1] - alNew[1];


      vectProdOz = vectAlBe[0] * vectAlTh[1]  -  vectAlBe[1] * vectAlTh[0];

      if( vectProdOz < 0. )
      {
         //   The vertex "beNew" is NOT in triangle "alNew - gaNew - thNew".

         qAngType = 0;

         //   Pseudo case. Anyway I need to find some solution. So

         firVfirT[0] = alNew[0];   firVfirT[1] = alNew[1];

         secVfirT[0] = beNew[0];   secVfirT[1] = beNew[1];

         thiVfirT[0] = thNew[0];   thiVfirT[1] = thNew[1];

         //   Second triangle.

         firVsecT[0] = beNew[0];   firVsecT[1] = beNew[1];

         secVsecT[0] = thNew[0];   secVsecT[1] = thNew[1];

         thiVsecT[0] = gaNew[0];   thiVsecT[1] = gaNew[1];

         return qAngType;
      }

      if( vectProdOz >= 0. )
      {
         //  It's all write. We have a good concave quadrangle.

         //   Now let's compute all vertices which I need.

         qAngType = 2;

         //   First triangle.

         firVfirT[0] = alNew[0];   firVfirT[1] = alNew[1];

         secVfirT[0] = beNew[0];   secVfirT[1] = beNew[1];

         thiVfirT[0] = thNew[0];   thiVfirT[1] = thNew[1];

         //   Second triangle.

         firVsecT[0] = beNew[0];   firVsecT[1] = beNew[1];

         secVsecT[0] = thNew[0];   secVsecT[1] = thNew[1];

         thiVsecT[0] = gaNew[0];   thiVsecT[1] = gaNew[1];

         return qAngType;
      }

   }   //   "if(  (  (alNew[0] - AcsP[0])*(gaNew[0] - AcsP[0])  )  <=  0.  )".   //   Last second case of second criterion.

}   //   end of "if (  (  (beNew[1] - AcrP[1])*(thNew[1] - AcrP[1])  )  >  0.  )"



//  Now let's consider SECOND case 5.b "(  (beNew[1] - AcrP[1])*(thNew[1] - AcrP[1])  )  <= 0."

if (  (  (beNew[1] - AcrP[1])*(thNew[1] - AcrP[1])  )  <= 0.  )
{
   //   Vertices "beNew" and "thNew" are on the different sides of across point.

   //   Second criterion. Is across point between "alNew" and "gaNew" vertices by Ox-axis?

   //   Second criterion. First case.

   if(  (  (alNew[0] - AcrP[0])*(gaNew[0] - AcrP[0])  )  >  0.  )
   {
      //  It means the across point IS NOT between "alNew" and "gaNew" vertices by Ox-axis?

      //   O.K. the quadrangle IS NOT CONVEX. Is it concave or pseudo? Third criterion.

      vectBeGa[0]  =  gaNew[0] - beNew[0];

      vectBeGa[1]  =  gaNew[1] - beNew[1];

      vectBeAl[0]  =  alNew[0] - beNew[0];

      vectBeAl[1]  =  alNew[1] - beNew[1];


      vectProdOz  =  vectBeGa[0] * vectBeAl[1]  -  vectBeGa[1] * vectBeAl[0];

      if( vectProdOz >= 0. )
      {
         qAngType = 2;

         //   The quadrangle is concave. First triangle.

         firVfirT[0] = alNew[0];   firVfirT[1] = alNew[1];

         secVfirT[0] = beNew[0];   secVfirT[1] = beNew[1];

         thiVfirT[0] = gaNew[0];   thiVfirT[1] = gaNew[1];

         //   Second triangle.

         firVsecT[0] = alNew[0];   firVsecT[1] = alNew[1];

         secVsecT[0] = thNew[0];   secVsecT[1] = thNew[1];

         thiVsecT[0] = gaNew[0];   thiVsecT[1] = gaNew[1];

         return qAngType;
      }

      if( vectProdOz < 0. )
      {
         qAngType = 0;

         //   This concave quadrangle do has NO write anticlockwise vertices sequence order. It's pseudo.

         firVfirT[0] = alNew[0];   firVfirT[1] = alNew[1];

         secVfirT[0] = beNew[0];   secVfirT[1] = beNew[1];

         thiVfirT[0] = gaNew[0];   thiVfirT[1] = gaNew[1];

         //   Second triangle.

         firVsecT[0] = alNew[0];   firVsecT[1] = alNew[1];

         secVsecT[0] = thNew[0];   secVsecT[1] = thNew[1];

         thiVsecT[0] = gaNew[0];   thiVsecT[1] = gaNew[1];

         return qAngType;
      }
   }   //   end of "if(  (  (alNew[0] - AcrP[0])*(gaNew[0] - AcsP[0])  )  >  0.  )". First case of second criterion.


   //   Second criterion. Second case.

   if(  (  (alNew[0] - AcrP[0])*(gaNew[0] - AcrP[0])  )  <=  0.  )
   {
      //   O.K. the quadrangle is convex. Is it has the same vertices sequence order.

      vectAlBe[0]  =  beNew[0] - alNew[0];

      vectAlBe[1]  =  beNew[1] - alNew[1];

      vectAlTh[0]  =  thNew[0] - alNew[0];

      vectAlTh[1]  =  thNew[1] - alNew[1];


      vectProdOz  =  vectAlBe[0] * vectAlTh[1]  -  vectAlBe[1] * vectAlTh[0];

      if( vectProdOz >= 0. )
      {
         qAngType = 1;

         //   Convex quadrangle DO HAS WRITE anticlockwise vertices sequence order. It's convex.

         firVfirT[0] = alNew[0];   firVfirT[1] = alNew[1];

         secVfirT[0] = beNew[0];   secVfirT[1] = beNew[1];

         thiVfirT[0] = gaNew[0];   thiVfirT[1] = gaNew[1];

         //   Second triangle.

         firVsecT[0] = alNew[0];   firVsecT[1] = alNew[1];

         secVsecT[0] = thNew[0];   secVsecT[1] = thNew[1];

         thiVsecT[0] = gaNew[0];   thiVsecT[1] = gaNew[1];

         return qAngType;
      }

      if( vectProdOz < 0. )
      {
         qAngType = 0;

         //   Convex quadrangle do has NO write anticlockwise vertices sequence order. But it's convex.

         firVfirT[0] = alNew[0];   firVfirT[1] = alNew[1];

         secVfirT[0] = beNew[0];   secVfirT[1] = beNew[1];

         thiVfirT[0] = gaNew[0];   thiVfirT[1] = gaNew[1];

         //   Second triangle.

         firVsecT[0] = alNew[0];   firVsecT[1] = alNew[1];

         secVsecT[0] = thNew[0];   secVsecT[1] = thNew[1];

         thiVsecT[0] = gaNew[0];   thiVsecT[1] = gaNew[1];

         return qAngType;
      }

   }   //   end of "if(  (  (alNew[0] - AcrP[0])*(gaNew[0] - AcsP[0])  )  <=  0.  )". //   Second case of second criterion.

}



// That's all.

return qAngType;
}







double spaceVolumeInPrevTL(
              double par_a,                           //   -  Item of left and right setback (parameter "a" in test).
              double par_b,                           //   -  Item of second parameter from "u_funcion".
              //
              double lbDom,                           //   -  Left and right boundaries of rectangular domain.
              double rbDom,
              //
              double bbDom,                           //   -  Botton and upper boundaries of rectangular domain.
              double ubDom,
              //
              double tau,
              double iCurrTL,                         //   -  Index of current time layer. Necessary for velocity.
              //
              int iOfOXN,                             //   -  Index of current OX node.
              double * masOX,                         //   -  Massive of OX steps. Dimension = numOfOXSt +1.
              int numOfOXSt,                          //   -  Number of OX steps.
              //
              int iOfOYN,                             //   -  Index of current OY node.
              double * masOY,                         //   -  Massive of OY steps. Dimension = numOfOYSt +1.
              int numOfOYSt,                          //   -  Number of OY steps.
              //
              bool & isSpVolComp,                     //   -  Is space volume computed?
              double * rhoInPrevTL_asV )
{


/*
double ** rhoInPrevTL;

rhoInPrevTL  =  new double * [ numOfOXSt +1 ];

for( int j=0; j< numOfOXSt +1; j++ )
{
   rhoInPrevTL[j] = new double [numOfOYSt +1];
}



for( int k=0; k< numOfOYSt +1; k++ )
{
   for( int j=0; j< numOfOXSt +1; j++ )
   {
      rhoInPrevTL[ j ][ k ]  =  rhoInPrevTL_asV[ (numOfOXSt +1)*k + j ];
   }
}
*/




double firVfirT[2], secVfirT[2], thiVfirT[2];         //   -  First, second and third vertices of first triangle.

double firVsecT[2], secVsecT[2], thiVsecT[2];         //   -  First, second and third vertices of second triangle.


double spVolInPrevTL = 0.;                            //   -  Space value in previous time level which we are computing.

int qAngType;                                         //   -  Type of quadrangle: 0 - pseudo; 1 - convex; 2 - concave;

double buf_D;



//   Let's understand what type of quadrangle we have.

qAngType = quadrAngleType(
              par_a, par_b,                           //   -  Analitycal solution parameters.
              //
              lbDom, rbDom,                           //   -  Left and right boundaries of rectangular domain.
              //
              bbDom, ubDom,                           //   -  Botton and upper boundaries of rectangular domain.
              //
              tau, iCurrTL,                           //   -  Time data. Necessary for velocity.
              //
              iOfOXN, masOX, numOfOXSt,               //   -  OX data.
              //
              iOfOYN, masOY, numOfOYSt,               //   -  OY data.
              //
              firVfirT, secVfirT, thiVfirT,           //   -  Vertices of first triangle.
              //
              firVsecT, secVsecT, thiVsecT );         //   -  Vertices of second triangle.



//   We consider only two cases.

if( qAngType != 1 )
{
   cout<<"\n\nqAngType != 1 "<<endl;

   isSpVolComp = false;

/*
   for( int j=0; j< numOfOXSt +1; j++ )
   {
      delete rhoInPrevTL[j];;
   }

   delete rhoInPrevTL;
*/

   return -1.;
}



//

if( qAngType == 1 )
{

   buf_D = integUnderUnunifTr(
              par_a, par_b,                           //   -  Analitycal solution parameters.
              //
              lbDom, rbDom,                           //   -  Left and right boundaries of rectangular domain.
              //
              bbDom, ubDom,                           //   -  Botton and upper boundaries of rectangular domain.
              //
              tau, iCurrTL,                           //   -  Index of current time layer.
              //
              firVfirT, secVfirT, thiVfirT,           //   -  Vertices of first triangle.
              //
              masOX, numOfOXSt,                       //   -  Number of OX steps.
              //
              masOY, numOfOYSt,                       //   -  Number of OY steps.
              //
              rhoInPrevTL_asV );




/*
   double Sq, Per, a, b, c;

   a = sqrt(  ( secVfirT[0] - firVfirT[0] )*( secVfirT[0] - firVfirT[0] )  +  ( secVfirT[1] - firVfirT[1] )*( secVfirT[1] - firVfirT[1] )   );

   b = sqrt(  ( thiVfirT[0] - secVfirT[0] )*( thiVfirT[0] - secVfirT[0] )  +  ( thiVfirT[1] - secVfirT[1] )*( thiVfirT[1] - secVfirT[1] )   );

   c = sqrt(  ( firVfirT[0] - thiVfirT[0] )*( firVfirT[0] - thiVfirT[0] )  +  ( firVfirT[1] - thiVfirT[1] )*( firVfirT[1] - thiVfirT[1] )   );

   Per = (a + b + c) /2.;

   Sq = sqrt( fabs(  Per * (Per - a) * (Per - b) * (Per - c) )  );
*/




   spVolInPrevTL = buf_D;



   buf_D = integUnderUnunifTr(
              par_a, par_b,                           //   -  Analitycal solution parameters.
              //
              lbDom, rbDom,                           //   -  Left and right boundaries of rectangular domain.
              //
              bbDom, ubDom,                           //   -  Botton and upper boundaries of rectangular domain.
              //
              tau, iCurrTL,                           //   -  Index of current time layer.
              //
              firVsecT, secVsecT, thiVsecT,           //   -  Vertices of first triangle.
              //
              masOX, numOfOXSt,                       //   -  Number of OX steps.
              //
              masOY, numOfOYSt,                       //   -  Number of OY steps.
              //
              rhoInPrevTL_asV );




/*
   a = sqrt(  ( secVsecT[0] - firVsecT[0] )*( secVsecT[0] - firVsecT[0] )  +  ( secVsecT[1] - firVsecT[1] )*( secVsecT[1] - firVsecT[1] )   );

   b = sqrt(  ( thiVsecT[0] - secVsecT[0] )*( thiVsecT[0] - secVsecT[0] )  +  ( thiVsecT[1] - secVsecT[1] )*( thiVsecT[1] - secVsecT[1] )   );

   c = sqrt(  ( firVsecT[0] - thiVsecT[0] )*( firVsecT[0] - thiVsecT[0] )  +  ( firVsecT[1] - thiVsecT[1] )*( firVsecT[1] - thiVsecT[1] )   );

   Per = (a + b + c) /2.;

   Sq = sqrt( fabs(  Per * (Per - a) * (Per - b) * (Per - c) )  );
*/




   spVolInPrevTL = spVolInPrevTL + buf_D ;

   isSpVolComp = true;


/*
   for( int j=0; j< numOfOXSt +1; j++ )
   {
      delete rhoInPrevTL[j];;
   }

   delete rhoInPrevTL;
*/


   return spVolInPrevTL;
}




/*
for( int j=0; j< numOfOXSt +1; j++ )
{
   delete rhoInPrevTL[j];;
}

delete rhoInPrevTL;
*/


return spVolInPrevTL;
}




