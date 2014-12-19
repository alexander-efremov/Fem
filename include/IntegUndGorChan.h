


extern double integOfChan_SLLeftSd(                   //   -  The domain is Channel with Slant Line on the left side.
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
              double *bv,   int wTrPCI,               //   -  Where travel point current (botton vertex) is.
              double *uv,   int wTrPNI,               //   -  Where travel point next (upper vertex) is.
              //
              int * indCurSqOx,                       //   -  Index by OX axis where bv and uv are.
              //
              double rb,  int * indRB,                //   -  Right boundary by Ox. Index by OX axis where rb is.
              //
              int * indCurSqOy,                       //   -  Index of current square by Oy axis.
              //
              double * masOX,                         //   -  Massive of OX steps. Dimension = numOfOXSt +1.
              int numOfOXSt,                          //   -  Number of OX steps.
              //
              double * masOY,                         //   -  Massive of OY steps. Dimension = numOfOYSt +1.
              int numOfOYSt,                          //   -  Number of OY steps.
              //
              double * rhoInPrevTL_asV );





extern double integOfChan_SLRightSd(                  //   -  The domain is Channel with Slant Line on the right side.
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
              double *bv,   int wTrPCI,               //   -  Where travel point current (botton vertex) is.
              double *uv,   int wTrPNI,               //   -  Where travel point next (upper vertex) is.
              //
              int * indCurSqOx,                       //   -  Index by OX axis where bv and uv are.
              //
              double lb,  int * indLB,                //   -  Left boundary by Ox. Index by OX axis where lb is.
              //
              int * indCurSqOy,                       //   -  Index of current square by Oy axis.
              //
              double * masOX,                         //   -  Massive of OX steps. Dimension = numOfOXSt +1.
              int numOfOXSt,                          //   -  Number of OX steps.
              //
              double * masOY,                         //   -  Massive of OY steps. Dimension = numOfOYSt +1.
              int numOfOYSt,                          //   -  Number of OY steps.
              //
              double * rhoInPrevTL_asV );


