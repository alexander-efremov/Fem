


#include <math.h>

const double C_pi = 3.14159265358979323846264338327;





double u_function(
              double par_b,                           //   -  Item of second parameter from "u_funcion" or "v_funcion".
              //
              double lbDom,                           //   -  Left and right boundaries of rectangular domain.
              double rbDom,
              //
              double bbDom,                           //   -  Botton and upper boundaries of rectangular domain.
              double ubDom,
              //
              double t, double x, double y )
{



double  u;



//  Test 1.0

//   u  =  x * (1. + t) /10. * y * (1.-y);

//   u  =  u  *  (1.+t) /10. * (y - ubDom) * (y - bbDom);

//   u  =  par_b * atan ( u );



//   Test 2.0

u  =  par_b * y * (1.-y) * (  C_pi /2. + atan( -x )  );



//   Test 3.0

//   u  =  par_b * y * (1.-y) * atan( x );





return u;
}







double v_function(
              double par_b,                           //   -  Item of second parameter from "u_funcion" or "v_funcion".
              //
              double lbDom,                           //   -  Left and right boundaries of rectangular domain.
              double rbDom,
              //
              double bbDom,                           //   -  Botton and upper boundaries of rectangular domain.
              double ubDom,
              //
              double t, double x, double y )
{



double v  =  (x - lbDom) * (x - rbDom) * (1.+t) /10. * (y - ubDom) * (y - bbDom);

v  =  atan( v );





return v;
}







double analytSolut(
              double par_a,
              //
              double lbDom,                           //   -  Left and right boundaries of rectangular domain.
              double rbDom,
              //
              double bbDom,                           //   -  Botton and upper boundaries of rectangular domain.
              double ubDom,
              //
              double t, double x, double y )
{



double anSol = 1.1  +  sin( t *x *y);                 //   5. * t  +  cos( x *y ) + 10. * y;   // 1.1  +  sin( t *x *y);





return anSol;
}







double initDataOfSol(
              double par_a,                           //   -  Item of left and right setback (parameter "a" in test).
              //
              double lbDom,                           //   -  Left and right boundaries of rectangular domain.
              double rbDom,
              //
              double bbDom,                           //   -  Botton and upper boundaries of rectangular domain.
              double ubDom,
              //
              int iOfOXN,                             //   -  Index of current OX node.
              double *masOX,                          //   -  Massive of abscissa grid steps. Dimension = numOfOxSt +1.
              //
              int iOfOYN,                             //   -  Index of current OY node.
              double *masOY )                         //   -  Massive of ordinate grid steps. Dimension = numOfOySt +1.
{



double rho;

double x = masOX[ iOfOXN ];

double y = masOY[ iOfOYN ];



rho  =  analytSolut( par_a,
              //
              lbDom, rbDom,
              //
              bbDom, ubDom,
              //
              0., x, y );





return rho;
}







double leftBound(
              double par_a,                           //   -  Item of left and right setback (parameter "a" in test).
              //
              double lbDom,                           //   -  Left and right boundaries of rectangular domain.
              double rbDom,
              //
              double bbDom,                           //   -  Botton and upper boundaries of rectangular domain.
              double ubDom,
              //
              double t,
              double y )
{



double lb;

lb  =  analytSolut( par_a,
              //
              lbDom, rbDom,
              //
              bbDom, ubDom,
              //
              t, lbDom, y );





return lb;
}







double rightBound(
              double par_a,                           //   -  Item of left and right setback (parameter "a" in test).
              //
              double lbDom,                           //   -  Left and right boundaries of rectangular domain.
              double rbDom,
              //
              double bbDom,                           //   -  Botton and upper boundaries of rectangular domain.
              double ubDom,
              //
              double t,
              double y )
{



double rb;

rb  =  analytSolut( par_a,
              //
              lbDom, rbDom,
              //
              bbDom, ubDom,
              //
              t, rbDom, y );





return rb;
}







double bottonBound(
              double par_a,                           //   -  Item of left and right setback (parameter "a" in test).
              //
              double lbDom,                           //   -  Left and right boundaries of rectangular domain.
              double rbDom,
              //
              double bbDom,                           //   -  Botton and upper boundaries of rectangular domain.
              double ubDom,
              //
              double t,
              double x )
{



double bb;

bb  =  analytSolut( par_a,
              //
              lbDom, rbDom,
              //
              bbDom, ubDom,
              //
              t, x, bbDom );





return bb;
}







double upperBound(
              double par_a,                           //   -  Item of left and right setback (parameter "a" in test).
              //
              double lbDom,                           //   -  Left and right boundaries of rectangular domain.
              double rbDom,
              //
              double bbDom,                           //   -  Botton and upper boundaries of rectangular domain.
              double ubDom,
              //
              double t,
              double x )
{



double ub;

ub  =  analytSolut( par_a,
              //
              lbDom, rbDom,
              //
              bbDom, ubDom,
              //
              t, x, ubDom );





return ub;
}







double f_function(                                    //   -  It's item of right part of differential equation.
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
              int iOfOXN,                             //   -  Index of current OX node.
              double *masOX,                          //   -  Massive of OX steps. Dimension = numOfOXSt +1.
              int numOfOXSt,                          //   -  Number of OX steps (segments).
               //
              int iOfOYN,                             //   -  Index of current OY node.
              double *masOY,                          //   -  Massive of OY steps. Dimension = numOfOYSt +1.
              int numOfOYSt )                         //   -  Number of OY steps (segments).
{



double t  =  tau * iCurrTL;

double x  =  masOX[ iOfOXN ];

double y  =  masOY[ iOfOYN ];


double arg_u  =   x * (1.+t) /10. * y * (1. - y);

double arg_v  =  (x - lbDom) * (x - rbDom) * (1.+t) /10. * (y - ubDom) * (y - bbDom);


double rho, dRhoDT, dRhoDX, dRhoDY;

double u, duDX;

double v, dvDY;

double f;



rho  =  analytSolut(
              par_a,
              //
              lbDom, rbDom,
              //
              bbDom, ubDom,
              //
              t, x, y );

dRhoDT  =  x * y * cos( t*x*y );

dRhoDX  =  t * y * cos( t*x*y );

dRhoDY  =  t * x * cos( t*x*y );



u  =  u_function(
              par_b,
              //
              lbDom, rbDom,
              //
              bbDom, ubDom,
              //
              t, x, y );

//   Test 1.0:   duDX  =  par_b * (1.+t) /10. * y * (1.-y)  /  ( 1.  +  arg_u * arg_u );

//   Test 2.0:   duDX  = -par_b * y * (1.-y)  /  ( 1.  +  x * x );

//   Test 3.0    duDX  =  par_b * y * (1.-y)  /  (1.  +  x*x );

duDX  = -par_b * y * (1.-y)  /  ( 1.  +  x * x );



v  =  v_function(
              par_b,
              //
              lbDom, rbDom,
              //
              bbDom, ubDom,
              //
              t, x, y );

dvDY  =  (x - lbDom) * (x - rbDom) * (1.+t) /10. * (y - bbDom + y - ubDom);

dvDY  =  dvDY  /  ( 1.  +  arg_v * arg_v );



f  =  dRhoDT   +   rho * duDX   +   u * dRhoDX   +   rho * dvDY   +   v * dRhoDY;





return f;
}







double ErrorWithAnalSol(
              double par_a,
              //
              double lbDom,                           //   -  Left and right boundaries of rectangular domain.
              double rbDom,
              //
              double bbDom,                           //   -  Botton and upper boundaries of rectangular domain.
              double ubDom,
              //
              double t,
              double x,
              double y,
              //
              double rho )                            //   -  Numerical solution.
{



double anSol;

//   Analytical solution.

anSol  =  analytSolut(
              par_a,
              //
              lbDom, rbDom,
              //
              bbDom, ubDom,
              //
              t, x, y );


anSol  =  anSol  - rho;





return anSol;
}







double MaxModItemOfVect( double *vect, int dim )
{



double maxItemOfComp;



maxItemOfComp = fabs( vect[0] );

for( int j=0; j< dim; j++ )
{
	if( maxItemOfComp < fabs(vect[j]) ){  maxItemOfComp = fabs( vect[j] ); }
}





return maxItemOfComp;
}







double MaxModItemOfMatr( double **mat, int dim1, int dim2 )
{



double maxItemOfComp;



maxItemOfComp = fabs( mat[0][0] );

for( int j=0; j< dim1; j++ )
{
   for( int k=0; k < dim2; k++ )
   {
      if( maxItemOfComp < fabs(mat[j][k]) )
      {
         maxItemOfComp = fabs( mat[j][k] );
      }
   }
}





return maxItemOfComp;
}







double normOfMatrAtL2(
              double *masOX,                          //   -  Massive of OX grid nodes. Dimension = dimOX.
              int dimOX,
              //
              double *masOY,                          //   -  Massive of OY grid nodes. Dimension = dimOY.
              int dimOY,
              //
              double ** mat )
{



double norm = 0.;

double h  =  masOX[1]  -  masOX[0];



for( int j=1; j< dimOX -1; j++ )
{
   for( int k=1; k< dimOY -1; k++ )
   {
      norm  +=  mat[j][k] * mat[j][k];
   }
	//   norm  +=  ( fabs(vect[j])  +  fabs( vect[j+1] ) )/2. * h;
}



norm = h * sqrt( norm );





return norm;
}







double normOfMatrAtL1(
              double *masOX,                          //   -  Massive of OX grid nodes. Dimension = dimOX.
              int dimOX,
              //
              double *masOY,                          //   -  Massive of OY grid nodes. Dimension = dimOY.
              int dimOY,
              //
              double ** mat )
{



double norm = 0.;

double h  =  masOX[1]  -  masOX[0];



for( int j=1; j< dimOX -1; j++ )
{
   for( int k=1; k< dimOY -1; k++ )
   {
      norm  +=  fabs( mat[j][k] );
   }
}



norm = h * h * norm ;





return norm;
}







double normOfMatrAtL1_asV(
              double *masOX,                          //   -  Massive of OX grid nodes. Dimension = dimOX.
              int dimOX,
              //
              double *masOY,                          //   -  Massive of OY grid nodes. Dimension = dimOY.
              int dimOY,
              //
              double * mat_asV )
{



double norm = 0.;

double hx  =  masOX[1]  -  masOX[0];

double hy  =  masOY[1]  -  masOY[0];



for( int k=1; k< dimOY -1; k++ )
{
   for( int j=1; j< dimOX -1; j++ )
   {
      norm  +=  fabs( mat_asV[ dimOX*k + j ] );
   }
}



norm = hx * hy * norm ;





return norm;
}






/*
double normOfVectAtL2(
              double *masOx,                          //   -  Massive of grid steps. Dimension = numOfGridSteps +1.
              int numOfGridSteps,
              //
              double *vect )
{



double norm = 0.;

double h  =  masOx[1]  -  masOx[0];



for( int j=0; j< numOfGridSteps; j++ )
{
	//   norm  +=  ( fabs(vect[j])  +  fabs( vect[j+1] ) )/2. * h;

	norm  +=  ( vect[j]*vect[j]  +  vect[j+1] * vect[j+1] )/2. * h;
}



norm = sqrt( norm );





return norm;
}
*/
