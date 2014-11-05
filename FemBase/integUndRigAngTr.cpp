


#include <math.h>

//   #include "InitAndBoundaryData.h"                      //   -  Initial and Boundary data.

#include "IntegUndGorChan.h"



double integUnderRigAngTr_BottLeft(
              double par_a,                           //   -  Solution parameter.
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
              double *bv,
              double *uv,
              //
              double * masOX,                         //   -  Massive of OX steps. Dimension = numOfOXSt +1.
              int numOfOXSt,                          //   -  Number of OX steps.
              //
              double * masOY,                         //   -  Massive of OY steps. Dimension = numOfOYSt +1.
              int numOfOYSt,                          //   -  Number of OY steps.
              //
              double * rhoInPrevTL_asV )
{



//   return ( fabs( (uv[1] - bv[1]) * (bv[0] - uv[0]) /2.) );



double trPC[2];                                       //   -  Travel point current;

int wTrPCI = 0;                                       //   -  Where travel point current is?


double trPN[2];                                       //   -  Travel point next;

int wTrPNI = 0;                                       //   -  Where travel point next is?


double ang;                                           //   -  Angle of slant line. Should be greater zero.


int indCurSqOx[2], indCurSqOy[2];                     //   -  Index of current square by Ox and Oy axes.

int indRB[2];                                         //   -  Index of right boundary.

double distOx, distOy;                                //   -  Distance to near Ox and Oy straight lines.


bool isTrDone = false;                                //   -  Is travel done.


double hx = masOX[1] - masOX[0];

double hy = masOY[1] - masOY[0];


double integOfBottTr = 0.;                            //   -  Value which we are computing.

double buf_D;





//   Initial data.

trPC[0] = bv[0];   trPC[1] = bv[1];

if(  ( fabs(bv[0] - uv[0]) )  <  1.e-12  )
{
   //   This triangle has very small width. I guess further computation isn't correct.

   return fabs(bv[0] - uv[0]);
}



ang = (uv[1] - bv[1]) / (bv[0] - uv[0]);

if(  fabs(ang)  <  1.e-12  )
{
   //   This triangle has very small height. I guess further computation isn't correct.

   return fabs(ang);
}





indCurSqOx[0] = (int)(  (trPC[0] - 1.e-14) /hx);      //   -  If trPC[0] is in grid edge I want it will be between in the left side of indCurSqOx[1].

if( (trPC[0] - 1.e-14) <= 0 ){ indCurSqOx[0] -= 1; }  //   -  The case when "trPC[0]" ia negative.

indCurSqOx[1] = indCurSqOx[0] +1;                     //   -  It's important only in rare case then trPC is in grid edge.


indRB[0] = indCurSqOx[0];

indRB[1] = indRB[0] +1;


indCurSqOy[0] = (int)(  (trPC[1] + 1.e-14) /hy);      //   -  If trPC[1] is in grid edge I want it will be between indCurSqOx[0] and indCurSqOx[1].

if( (trPC[1] + 1.e-14) <= 0 ){ indCurSqOy[0] -= 1; }  //   -  The case when "trPC[0]" ia negative.

indCurSqOy[1] = indCurSqOy[0] +1;                     //   -  It's important only in rare case then trPC is in grid edge.



if( indCurSqOx[0] >= 0)
{
   distOx = trPC[0]  -  masOX[ indCurSqOx[0] ];
}

if( indCurSqOx[0] < 0 )
{
   distOx = fabs( trPC[0]  -  hx * indCurSqOx[0] );
}


if( indCurSqOy[1] >= 0 )
{
   distOy = masOY[ indCurSqOy[1] ]  -  trPC[1];
}


if( indCurSqOy[1] < 0 )
{
   distOy = fabs( hy * indCurSqOy[1]  -  trPC[1] );
}



do{
   //   a. First case.

   if( (distOy /distOx) <= ang )
   {
      //   Across with straight line parallel Ox axis.

      wTrPNI = 1;


      if( indCurSqOy[1] >= 0)
      {
         trPN[1] = masOY[ indCurSqOy[1] ];
      }

      if( indCurSqOy[1] < 0)
      {
         trPN[1] = hy * indCurSqOy[1];
      }


      trPN[0] = bv[0] - (trPN[1] - bv[1]) /ang;
   }


   //   b. Second case.

   if( (distOy /distOx) > ang )
   {
      //   Across with straight line parallel Oy axis.

      wTrPNI = 2;


      if( indCurSqOx[0] >= 0 )
      {
         trPN[0]  =  masOX[ indCurSqOx[0] ];
      }

      if( indCurSqOx[0] < 0 )
      {
         trPN[0]  =  hx * indCurSqOx[0];
      }


      trPN[1]  =  bv[1]  -  ang * (trPN[0] - bv[0]);
   }



   //   c. Cheking.

   if(  trPN[0]  <  (uv[0] + 1.e-14)  )
   {
      trPN[0] = uv[0];

      trPN[1] = uv[1];

      isTrDone = true;

      wTrPNI = 0;
   }



   //   d. Integration.

   buf_D = integOfChan_SLLeftSd(                      //   -  The domain is Channel with Slant Line on the left side.
              par_a,                                  //   -  Solution parameter.
              //
              lbDom, rbDom,
              //
              bbDom, ubDom,
              //
              tau, iCurrTL,                           //   -  Index of current time layer.
              //
              trPC,  wTrPCI,                          //   -  double *bv,
              trPN,  wTrPNI,                          //   -  double *uv,
              //
              indCurSqOx,                             //   -  Indices where trPC and trPN are.
              //
              bv[0], indRB,                           //   -  double rb  =  Right boundary by Ox.
              //
              indCurSqOy,                             //   -  Index of current square by Oy axis.
              //
              masOX,                                  //   -  Massive of OX steps. Dimension = numOfOXSt +1.
              numOfOXSt,                              //   -  Number of OX steps.
              //
              masOY,                                  //   -  Massive of OY steps. Dimension = numOfOYSt +1.
              numOfOYSt,                              //   -  Number of OY steps.
              //
              rhoInPrevTL_asV );

   integOfBottTr = integOfBottTr + buf_D;



   //   e. Updating.

   if( isTrDone == false )
   {
      //   We will compute more. We need to redefine some values.

      wTrPCI = wTrPNI;

      trPC[0] = trPN[0];   trPC[1] = trPN[1];


      if( wTrPNI == 1)
      {
         indCurSqOy[0] += 1;

         indCurSqOy[1] += 1;
      }

      if( wTrPNI == 2)
      {
         indCurSqOx[0] -= 1;

         indCurSqOx[1] -= 1;
      }



      if( indCurSqOx[0] >= 0)
      {
         distOx = trPC[0]  -  masOX[ indCurSqOx[0] ];
      }

      if( indCurSqOx[0] < 0)
      {
         distOx = fabs( trPC[0]  -  hx * indCurSqOx[0] );
      }


      if( indCurSqOy[1] >= 0 )
      {
         distOy = masOY[ indCurSqOy[1] ]  -  trPC[1];
      }


      if( indCurSqOy[1] < 0 )
      {
         distOy = fabs( hy * indCurSqOy[1]  -  trPC[1] );
      }
   }

}
while( isTrDone == false );





return integOfBottTr;
}







double integUnderRigAngTr_UppLeft(
              double par_a,                           //   -  Solution parameter.
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
              double *bv,
              double *uv,
              //
              double * masOX,                         //   -  Massive of OX steps. Dimension = numOfOXSt +1.
              int numOfOXSt,                          //   -  Number of OX steps.
              //
              double * masOY,                         //   -  Massive of OY steps. Dimension = numOfOYSt +1.
              int numOfOYSt,                          //   -  Number of OY steps.
              //
              double * rhoInPrevTL_asV )
{



//   return ( fabs( (uv[1] - bv[1]) * (bv[0] - uv[0]) /2.) );



double trPC[2];                                       //   -  Travel point current;

int wTrPCI = 0;                                       //   -  Where travel point current is?


double trPN[2];                                       //   -  Travel point next;

int wTrPNI = 0;                                       //   -  Where travel point next is?


double ang;                                           //   -  Angle of slant line. Should be greater zero.


int indCurSqOx[2], indCurSqOy[2];                     //   -  Index of current square by Ox and Oy axes.

int indRB[2];                                         //   -  Index of right boundary.

double distOx, distOy;                                //   -  Distance to near Ox and Oy straight lines.


bool isTrDone = false;                                //   -  Is travel done.


double hx = masOX[1] - masOX[0];

double hy = masOY[1] - masOY[0];


double integOfUppTr = 0.;                             //   -  Value which we are computing.

double buf_D;





//   Initial data.

trPC[0] = bv[0];   trPC[1] = bv[1];

if(  ( fabs(bv[0] - uv[0]) )  <  1.e-12  )
{
   //   This triangle has very small width. I guess further computation isn't correct.

   return fabs(bv[0] - uv[0]);
}



ang = (uv[1] - bv[1]) / (uv[0] - bv[0]);

if(  fabs(ang)  <  1.e-12  )
{
   //   This triangle has very small height. I guess further computation isn't correct.

   return fabs(ang);
}




//   The follow equations are quite important.

indCurSqOx[0] = (int)(  (trPC[0] + 1.e-14) /hx);      //   -  If trPC[0] is in grid edge I want it will be in the right side.

if( (trPC[0] + 1.e-14) <= 0 ){ indCurSqOx[0] -= 1; }  //   -  The case when "trPC[0]" ia negative.

indCurSqOx[1] = indCurSqOx[0] +1;                     //   -  It's important only in rare case then trPC is in grid edge.


indCurSqOy[0] = (int)(  (trPC[1] + 1.e-14) /hy);      //   -  If trPC[1] is in grid edge I want it will be in the upper square.

if( (trPC[1] + 1.e-14) <= 0 ){ indCurSqOy[0] -= 1; }  //   -  The case when "trPC[0]" ia negative.

indCurSqOy[1] = indCurSqOy[0] +1;


indRB[0] = (int)(  (uv[0] - 1.e-14) /hy);             //   -  If uv[0] is in grid edge I want it will be in the left side.

if( (uv[0] - 1.e-14) <= 0 ){  indRB[0] -= 1;  }       //   -  The case when "trPC[0]" ia negative.

indRB[1] = indRB[0] +1;



if( indCurSqOx[1] >= 0)
{
   distOx = masOX[ indCurSqOx[1] ]  -  trPC[0];
}

if( indCurSqOx[1] < 0)
{
   distOx = fabs( hx * indCurSqOx[1]  -  trPC[0] );
}


if( indCurSqOy[1] >= 0 )
{
   distOy = masOY[ indCurSqOy[1] ]  -  trPC[1];
}

if( indCurSqOy[1] < 0 )
{
   distOy = fabs( hy * indCurSqOy[1]  -  trPC[1] );
}




do{
   //   a. First case.

   if( (distOy /distOx) <= ang )
   {
      //   Across with straight line parallel Ox axis.

      wTrPNI = 1;


      if( indCurSqOy[1] >= 0 )
      {
         trPN[1] = masOY[ indCurSqOy[1] ];
      }

      if( indCurSqOy[1] < 0 )
      {
         trPN[1] = hy * indCurSqOy[1];
      }


      trPN[0] = bv[0] + (trPN[1] - bv[1]) /ang;
   }

   //   b. Second case.

   if( (distOy /distOx) > ang )
   {
      //   Across with straight line parallel Oy axis.

      wTrPNI = 2;


      if( indCurSqOx[1] >= 0 )
      {
         trPN[0]  =  masOX[ indCurSqOx[1] ];
      }

      if( indCurSqOx[1] < 0 )
      {
         trPN[0]  =  hx * indCurSqOx[1];
      }


      trPN[1]  =  bv[1]  +  ang * (trPN[0] - bv[0]);
   }



   //   c. Cheking.

   if(  trPN[0]  >  (uv[0] - 1.e-14)  )
   {
      trPN[0] = uv[0];

      trPN[1] = uv[1];

      isTrDone = true;

      wTrPNI = 0;
   }



   //   d. Integration.

   buf_D = integOfChan_SLLeftSd(                      //   -  The domain is Channel with Slant Line on the left side.
              par_a,                                  //   -  Solution parameter.
              //
              lbDom, rbDom,
              //
              bbDom, ubDom,
              //
              tau, iCurrTL,                           //   -  Index of current time layer.
              //
              trPC,  wTrPCI,                          //   -  double *bv,
              trPN,  wTrPNI,                          //   -  double *uv,
              //
              indCurSqOx,                             //   -  Indices where trPC and trPN are.
              //
              uv[0], indRB,                           //   -  double rb  =  Right boundary by Ox.
              //
              indCurSqOy,                             //   -  Index of current square by Oy axis.
              //
              masOX,                                  //   -  Massive of OX steps. Dimension = numOfOXSt +1.
              numOfOXSt,                              //   -  Number of OX steps.
              //
              masOY,                                  //   -  Massive of OY steps. Dimension = numOfOYSt +1.
              numOfOYSt,                              //   -  Number of OY steps.
              //
              rhoInPrevTL_asV );

   integOfUppTr = integOfUppTr + buf_D;



   //   e. Updating.

   if( isTrDone == false )
   {
      //   We will compute more. We need to redefine some values.

      wTrPCI = wTrPNI;

      trPC[0] = trPN[0];   trPC[1] = trPN[1];


      if( wTrPNI == 1)
      {
         indCurSqOy[0] += 1;

         indCurSqOy[1] += 1;
      }

      if( wTrPNI == 2)
      {
         indCurSqOx[0] += 1;

         indCurSqOx[1] += 1;
      }


      if( indCurSqOx[1] >= 0)
      {
         distOx = fabs( masOX[ indCurSqOx[1] ]  -  trPC[0] );
      }

      if( indCurSqOx[1] < 0)
      {
         distOx = fabs( hx * indCurSqOx[1]  -  trPC[0] );
      }


      if( indCurSqOy[1] >= 0 )
      {
         distOy = fabs( masOY[ indCurSqOy[1] ]  -  trPC[1] );
      }

      if( indCurSqOy[1] < 0 )
      {
         distOy = fabs( hy * indCurSqOy[1]  -  trPC[1] );
      }
   }

}
while( isTrDone == false );





return integOfUppTr;
}







double integUnderRigAngTr_BottRight(
              double par_a,                           //   -  Solution parameter.
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
              double *bv,
              double *uv,
              //
              double * masOX,                         //   -  Massive of OX steps. Dimension = numOfOXSt +1.
              int numOfOXSt,                          //   -  Number of OX steps.
              //
              double * masOY,                         //   -  Massive of OY steps. Dimension = numOfOYSt +1.
              int numOfOYSt,                          //   -  Number of OY steps.
              //
              double * rhoInPrevTL_asV )
{



//   return ( fabs( (uv[1] - bv[1]) * (bv[0] - uv[0]) /2.) );



double trPC[2];                                       //   -  Travel point current;

int wTrPCI = 0;                                       //   -  Where travel point current is?


double trPN[2];                                       //   -  Travel point next;

int wTrPNI = 0;                                       //   -  Where travel point next is?


double ang;                                           //   -  Angle of slant line. Should be greater zero.


int indCurSqOx[2], indCurSqOy[2];                     //   -  Index of current square by Ox and Oy axes.

int indLB[2];                                         //   -  Index of left boundary.

double distOx, distOy;                                //   -  Distance to near Ox and Oy straight lines.


bool isTrDone = false;                                //   -  Is travel done.


double hx = masOX[1] - masOX[0];

double hy = masOY[1] - masOY[0];


double integOfBottTr = 0.;                            //   -  Value which we are computing.

double buf_D;





//   Initial data.

trPC[0] = bv[0];   trPC[1] = bv[1];

if(  ( fabs(bv[0] - uv[0]) )  <  1.e-12  )
{
   //   This triangle has very small width. I guess further computation isn't correct.

   return fabs(bv[0] - uv[0]);
}



ang = (uv[1] - bv[1]) / (uv[0] - bv[0]);

if(  fabs(ang)  <  1.e-12  )
{
   //   This triangle has very small height. I guess further computation isn't correct.

   return fabs(ang);
}





indCurSqOx[0] = (int)(  (trPC[0] + 1.e-14) /hx);      //   -  If trPC[0] is in grid edge I want it will be between in the right side.

if( (trPC[0] + 1.e-14) <= 0 ){ indCurSqOx[0] -= 1; }  //   -  The case when "trPC[0]" ia negative.

indCurSqOx[1] = indCurSqOx[0] +1;                     //   -  It's important only in rare case then trPC is in grid edge.



indLB[0] = indCurSqOx[0];

indLB[1] = indLB[0] +1;



indCurSqOy[0] = (int)(  (trPC[1] + 1.e-14) /hy);      //   -  If trPC[1] is in grid edge I want it will be in the upper side.

if( (trPC[1] + 1.e-14) <= 0 ){ indCurSqOy[0] -= 1; }  //   -  The case when "trPC[0]" ia negative.

indCurSqOy[1] = indCurSqOy[0] +1;                     //   -  It's important only in rare case then trPC is in grid edge.



if( indCurSqOx[1] >=0 )
{
   distOx = fabs( masOX[ indCurSqOx[1] ]  -  trPC[0] );
}

if( indCurSqOx[1] < 0 )
{
   distOx = fabs( hx * indCurSqOx[1]  -  trPC[0] );
}



if( indCurSqOy[1] >=0 )
{
   distOy = fabs( masOY[ indCurSqOy[1] ]  -  trPC[1] );
}

if( indCurSqOy[1] < 0 )
{
   distOy = fabs( hy * indCurSqOy[1]  -  trPC[1] );
}



do{
   //   a. First case.

   if( (distOy /distOx) <= ang )
   {
      //   Across with straight line parallel Ox axis.

      wTrPNI = 1;


      if( indCurSqOy[1] >=0 )
      {
         trPN[1] = masOY[ indCurSqOy[1] ];
      }

      if( indCurSqOy[1] < 0 )
      {
         trPN[1] = hy * indCurSqOy[1];
      }


      trPN[0] = bv[0] + (trPN[1] - bv[1]) /ang;
   }



   //   b. Second case.

   if( (distOy /distOx) > ang )
   {
      //   Across with straight line parallel Oy axis.

      wTrPNI = 2;


      if( indCurSqOx[1] >= 0 )
      {
         trPN[0]  =  masOX[ indCurSqOx[1] ];
      }

      if( indCurSqOx[1]  < 0 )
      {
         trPN[0]  =  hx * indCurSqOx[1];
      }


      trPN[1]  =  bv[1]  +  ang * (trPN[0] - bv[0]);
   }



   //   c. Cheking.

   if(  trPN[0]  >  (uv[0] - 1.e-14)  )               //   -  Without "fabs"!!!
   {
      trPN[0] = uv[0];

      trPN[1] = uv[1];

      isTrDone = true;

      wTrPNI = 0;
   }



   //   d. Integration.

   buf_D = integOfChan_SLRightSd(                     //   -  The domain is Channel with Slant Line on the Right side.
              par_a,                                  //   -  Solution parameter.
              //
              lbDom, rbDom,
              //
              bbDom, ubDom,
              //
              tau, iCurrTL,                           //   -  Index of current time layer.
              //
              trPC,  wTrPCI,                          //   -  double *bv,
              trPN,  wTrPNI,                          //   -  double *uv,
              //
              indCurSqOx,                             //   -  Indices where trPC and trPN are.
              //
              bv[0], indLB,                           //   -  double lb  =  Left boundary by Ox.
              //
              indCurSqOy,                             //   -  Index of current square by Oy axis.
              //
              masOX,                                  //   -  Massive of OX steps. Dimension = numOfOXSt +1.
              numOfOXSt,                              //   -  Number of OX steps.
              //
              masOY,                                  //   -  Massive of OY steps. Dimension = numOfOYSt +1.
              numOfOYSt,                              //   -  Number of OY steps.
              //
              rhoInPrevTL_asV );

   integOfBottTr = integOfBottTr + buf_D;



   //   e. Updating.

   if( isTrDone == false )
   {
      //   We will compute more. We need to redefine some values.

      wTrPCI = wTrPNI;

      trPC[0] = trPN[0];   trPC[1] = trPN[1];


      if( wTrPNI == 1)
      {
         indCurSqOy[0] += 1;

         indCurSqOy[1] += 1;
      }

      if( wTrPNI == 2)
      {
         indCurSqOx[0] += 1;

         indCurSqOx[1] += 1;
      }



      if( indCurSqOx[1] >=0 )
      {
         distOx = fabs( masOX[ indCurSqOx[1] ]  -  trPC[0] );
      }

      if( indCurSqOx[1] < 0 )
      {
         distOx = fabs( hx * indCurSqOx[1]  -  trPC[0] );
      }



      if( indCurSqOy[1] >=0 )
      {
         distOy = fabs( masOY[ indCurSqOy[1] ]  -  trPC[1] );
      }

      if( indCurSqOy[1] < 0 )
      {
         distOy = fabs( hy * indCurSqOy[1]  -  trPC[1] );
      }
   }

}
while( isTrDone == false );




return integOfBottTr;
}







double integUnderRigAngTr_UppRight(
              double par_a,                           //   -  Solution parameter.
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
              double *bv,
              double *uv,
              //
              double * masOX,                         //   -  Massive of OX steps. Dimension = numOfOXSt +1.
              int numOfOXSt,                          //   -  Number of OX steps.
              //
              double * masOY,                         //   -  Massive of OY steps. Dimension = numOfOYSt +1.
              int numOfOYSt,                          //   -  Number of OY steps.
              //
              double * rhoInPrevTL_asV )
{



//   return ( fabs( (uv[1] - bv[1]) * (bv[0] - uv[0]) /2.) );



double trPC[2];                                       //   -  Travel point current;

int wTrPCI = 0;                                       //   -  Where travel point current is?


double trPN[2];                                       //   -  Travel point next;

int wTrPNI = 0;                                       //   -  Where travel point next is?


double ang;                                           //   -  Angle of slant line. Should be greater zero.


int indCurSqOx[2], indCurSqOy[2];                     //   -  Index of current square by Ox and Oy axes.

int indLB[2];                                         //   -  Index of left boundary.

double distOx, distOy;                                //   -  Distance to near Ox and Oy straight lines.


bool isTrDone = false;                                //   -  Is travel done.


double hx = masOX[1] - masOX[0];

double hy = masOY[1] - masOY[0];


double integOfUppTr = 0.;                             //   -  Value which we are computing.

double buf_D;





//   Initial data.

trPC[0] = bv[0];   trPC[1] = bv[1];

if(  ( fabs(bv[0] - uv[0]) )  <  1.e-12  )
{
   //   This triangle has very small width. I guess further computation isn't correct.

   return fabs(bv[0] - uv[0]);
}



ang = (uv[1] - bv[1]) / (bv[0] - uv[0]);

if(  fabs(ang)  <  1.e-12  )
{
   //   This triangle has very small height. I guess further computation isn't correct.

   return fabs(ang);
}





indCurSqOx[0] = (int)(  (trPC[0] - 1.e-14) /hx);      //   -  If trPC[0] is in grid edge I want it will be between in the left side.

if( (trPC[0] - 1.e-14) <= 0 ){ indCurSqOx[0] -= 1; }  //   -  The case when "trPC[0]" ia negative.

indCurSqOx[1] = indCurSqOx[0] +1;                     //   -  It's important only in rare case then trPC is in grid edge.



indLB[0] = (int)( (uv[0] + 1.e-14) /hx);

if( (uv[0] + 1.e-14) <=0 ){   indLB[0] -= 1;  }       //   -  The case when "trPC[0]" ia negative.

indLB[1] = indLB[0] +1;



indCurSqOy[0] = (int)(  (trPC[1] + 1.e-14) /hy);      //   -  If trPC[1] is in grid edge I want it will be in the upper side.

if( (trPC[1] + 1.e-14) <= 0 ){ indCurSqOy[0] -= 1; }  //   -  The case when "trPC[0]" ia negative.

indCurSqOy[1] = indCurSqOy[0] +1;                     //   -  It's important only in rare case then trPC is in grid edge.



if( indCurSqOx[0] >= 0 )
{
   distOx = fabs( trPC[0]  -  masOX[ indCurSqOx[0] ] );
}

if( indCurSqOx[0] < 0 )
{
   distOx = fabs( trPC[0]  -  hx * indCurSqOx[0] );
}


if( indCurSqOy[1] >= 0 )
{
   distOy = fabs( masOY[ indCurSqOy[1] ]  -  trPC[1] );
}

if( indCurSqOy[1] < 0 )
{
   distOy = fabs( hy * indCurSqOy[1]  -  trPC[1] );
}



do{
   //   a. First case.

   if( (distOy /distOx) <= ang )
   {
      //   Across with straight line parallel Ox axis.

      wTrPNI = 1;


      if( indCurSqOy[1] >= 0 )
      {
         trPN[1] = masOY[ indCurSqOy[1] ];
      }

      if( indCurSqOy[1] < 0 )
      {
         trPN[1] = hy * indCurSqOy[1];
      }


      trPN[0] = bv[0] - (trPN[1] - bv[1]) /ang;
   }

   //   b. Second case.

   if( (distOy /distOx) > ang )
   {
      //   Across with straight line parallel Oy axis.

      wTrPNI = 2;


      if( indCurSqOx[0] >= 0 )
      {
         trPN[0]  =  masOX[ indCurSqOx[0] ];
      }

      if( indCurSqOx[0] < 0 )
      {
         trPN[0]  =  hx * indCurSqOx[0];
      }


      trPN[1]  =  bv[1]  -  ang * (trPN[0] - bv[0]);
   }



   //   c. Cheking.

   if(  trPN[0]  <  (uv[0] + 1.e-14)  )
   {
      trPN[0] = uv[0];

      trPN[1] = uv[1];

      isTrDone = true;

      wTrPNI = 0;
   }



   //   d. Integration.

   buf_D = integOfChan_SLRightSd(                     //   -  The domain is Channel with Slant Line on the Right side.
              par_a,                                  //   -  Solution parameter.
              //
              lbDom, rbDom,
              //
              bbDom, ubDom,
              //
              tau, iCurrTL,                           //   -  Index of current time layer.
              //
              trPC,  wTrPCI,                          //   -  double *bv,
              trPN,  wTrPNI,                          //   -  double *uv,
              //
              indCurSqOx,                             //   -  Indices where trPC and trPN are.
              //
              uv[0], indLB,                           //   -  double lb  =  Left boundary by Ox.
              //
              indCurSqOy,                             //   -  Index of current square by Oy axis.
              //
              masOX,                                  //   -  Massive of OX steps. Dimension = numOfOXSt +1.
              numOfOXSt,                              //   -  Number of OX steps.
              //
              masOY,                                  //   -  Massive of OY steps. Dimension = numOfOYSt +1.
              numOfOYSt,                              //   -  Number of OY steps.
              //
              rhoInPrevTL_asV );

   integOfUppTr = integOfUppTr + buf_D;



   //   e. Updating.

   if( isTrDone == false )
   {
      //   We will compute more. We need to redefine some values.

      wTrPCI = wTrPNI;

      trPC[0] = trPN[0];   trPC[1] = trPN[1];


      if( wTrPNI == 1)
      {
         indCurSqOy[0] += 1;

         indCurSqOy[1] += 1;
      }

      if( wTrPNI == 2)
      {
         indCurSqOx[0] -= 1;

         indCurSqOx[1] -= 1;
      }


      if( indCurSqOx[0] >= 0 )
      {
         distOx = fabs( trPC[0]  -  masOX[ indCurSqOx[0] ] );
      }

      if( indCurSqOx[0] < 0 )
      {
         distOx = fabs( trPC[0]  -  hx * indCurSqOx[0] );
      }


      if( indCurSqOy[1] >= 0 )
      {
         distOy = fabs( masOY[ indCurSqOy[1] ]  -  trPC[1] );
      }

      if( indCurSqOy[1] < 0 )
      {
         distOy = fabs( hy * indCurSqOy[1]  -  trPC[1] );
      }
   }

}
while( isTrDone == false );





return integOfUppTr;
}
