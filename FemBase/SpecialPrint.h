


extern bool printVectorBy3D(
              char *fileName,
              //
              int n,                                  //   -  характеристика сетки.
              int numOfCurrTimeIter,                  //   -  Ќа какой шаг по времени относ€тс€ эти данные.
              int numOfSolOrd,                        //   -  Solution order which we want to get.
              //
              double numOfTimeSteps,                  //   -  ќбщее число шагов по времени.
              //
              double *masOx,
              double *VectorOfData,
              //
              int dimOfVect );





extern bool printSurface_asV(
              char *fileName,                         //   -  char *fileName,
              //
              int grPar,                              //   -  Grid parameter. Number of abscissa grid steps.
              int iCurrTSt,                           //   -  Index of current time step (layer) IN WHICH we printing solution.
              int numOfSolOrd,                        //   -  Solution order which we want to print.
              //
              double numOfTSt,                        //   -  A number of time steps.
              //
              double *masOX,                          //   -  Massive of abscissa grid points. Dimension = numOfOXSt +1.
              int numOfOXSt,                          //   -  Number of abscissa grid steps.
              //
              double *masOY,                          //   -  Massive of ordinate grid points. Dimension = numOfOYSt +1.
              int numOfOYSt,                          //   -  Number of ordinate grid steps.
              //
              double *matData_asV );                  //   -  Matrix of data.
