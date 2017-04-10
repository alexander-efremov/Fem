


extern double u_function(
              double par_b,                           //   -  Item of second parameter from "u_funcion" or "v_funcion".
              //
              double lbDom,                           //   -  Left and right boundaries of rectangular domain.
              double rbDom,
              //
              double bbDom,                           //   -  Botton and upper boundaries of rectangular domain.
              double ubDom,
              //
              double t, double x, double y );





extern double v_function(
              double par_b,                           //   -  Item of second parameter from "u_funcion" or "v_funcion".
              //
              double lbDom,                           //   -  Left and right boundaries of rectangular domain.
              double rbDom,
              //
              double bbDom,                           //   -  Botton and upper boundaries of rectangular domain.
              double ubDom,
              //
              double t, double x, double y );




extern double analytSolut(
              double par_a,
              //
              double lbDom,                           //   -  Left and right boundaries of rectangular domain.
              double rbDom,
              //
              double bbDom,                           //   -  Botton and upper boundaries of rectangular domain.
              double ubDom,
              //
              double t, double x, double y );





extern double initDataOfSol(
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
              double *masOY );                        //   -  Massive of ordinate grid steps. Dimension = numOfOySt +1.





extern double leftBound(
              double par_a,                           //   -  Item of left and right setback (parameter "a" in test).
              //
              double lbDom,                           //   -  Left and right boundaries of rectangular domain.
              double rbDom,
              //
              double bbDom,                           //   -  Botton and upper boundaries of rectangular domain.
              double ubDom,
              //
              double t,
              double y );





extern double rightBound(
              double par_a,                           //   -  Item of left and right setback (parameter "a" in test).
              //
              double lbDom,                           //   -  Left and right boundaries of rectangular domain.
              double rbDom,
              //
              double bbDom,                           //   -  Botton and upper boundaries of rectangular domain.
              double ubDom,
              //
              double t,
              double y );




extern double bottonBound(
              double par_a,                           //   -  Item of left and right setback (parameter "a" in test).
              //
              double lbDom,                           //   -  Left and right boundaries of rectangular domain.
              double rbDom,
              //
              double bbDom,                           //   -  Botton and upper boundaries of rectangular domain.
              double ubDom,
              //
              double t,
              double x );





extern double upperBound(
              double par_a,                           //   -  Item of left and right setback (parameter "a" in test).
              //
              double lbDom,                           //   -  Left and right boundaries of rectangular domain.
              double rbDom,
              //
              double bbDom,                           //   -  Botton and upper boundaries of rectangular domain.
              double ubDom,
              //
              double t,
              double x );





extern double f_function(                             //   -  It's item of right part of differential equation.
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
              int numOfOYSt );                        //   -  Number of OY steps (segments).





extern double ErrorWithAnalSol(
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
              double rho );                           //   -  Numerical solution.





extern double MaxModItemOfVect( double *vect, int dim );





extern double MaxModItemOfMatr( double **mat, int dim1, int dim2 );





extern double normOfMatrAtL2(
              double *masOX,                          //   -  Massive of OX grid nodes. Dimension = dimOX.
              int dimOX,
              //
              double *masOY,                          //   -  Massive of OY grid nodes. Dimension = dimOY.
              int dimOY,
              //
              double ** mat );




extern double normOfMatrAtL1(
              double *masOX,                          //   -  Massive of OX grid nodes. Dimension = dimOX.
              int dimOX,
              //
              double *masOY,                          //   -  Massive of OY grid nodes. Dimension = dimOY.
              int dimOY,
              //
              double ** mat );





extern double norm_at_L_1(
              double *masOX,                          //   -  Massive of OX grid nodes. Dimension = dimOX.
              int dimOX,
              //
              double *masOY,                          //   -  Massive of OY grid nodes. Dimension = dimOY.
              int dimOY,
              //
              double * mat_asV );




/*
extern double normOfVectAtL2(
              double *masOx,                          //   -  Massive of grid steps. Dimension = numOfGridSteps +1.
              int numOfGridSteps,
              //
              double *vect );
*/
