


extern bool solByEqualVolumes(
              double par_a,                           //   -  Item of left and right setback (parameter "a" in test).
              double par_b,                           //   -  Item of second parameter from "u_funcion".
              //
              double lbDom,                           //   -  Left and right boundaries of rectangular domain.
              double rbDom,
              //
              double bbDom,                           //   -  Botton and upper boundaries of rectangular domain.
              double ubDom,
              //
              double tau,                             //   -  Time step.
              int numOfTSt,                           //   -  A number of time steps.
              //
              double *masOX,                          //   -  Massive of OX points. Dimension = numOfOXSt +1.
              int numOfOXSt,                          //   -  Number of OX steps.
              //
              double *masOY,                          //   -  Massive of OY points. Dimension = numOfOYSt +1.
              int numOfOYSt,                          //   -  Number of OY steps.
              //
              int numOfSolOrd,                        //   -  For print only. Solution order which we want to get.
              //
              double *rhoInCurrTL_asV );              //   -  Rho (solution) in Last Time Level which we will compute.





extern double* solByEqualVolWithVarStepPlusPrint(
              double par_a,                           //   -  Item of left and right setback (parameter "a" in test).
              double par_b,                           //   -  Item of second parameter from "u_funcion".
              //
              double lbDom,                           //   -  Left and right boundaries of rectangular domain.
              double rbDom,
              //
              double bbDom,                           //   -  Botton and upper boundaries of rectangular domain.
              double ubDom,
              //
              double tau,                             //   -  Time step.
              int numOfTSt,                           //   -  A number of time steps.
              //
              int numOfOXSt,                          //   -  Number of OX steps.
              //
              int numOfOYSt,                          //   -  Number of OY steps.
              //
              bool isTimeStShBeChan,                  //   -  Is time step shoulb be change?
              bool isGridStShBeChan,                  //   -  Is grid step shoulb be change?
              //
              int numOfGrStepLayer );                 //   -  How many computations with different grid steps we want to make.

