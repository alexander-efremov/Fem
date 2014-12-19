


extern double integUnderBottTr(
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
              double * rhoInPrevTL_asV );





extern double integUnderUpperTr(
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
              double * rhoInPrevTL_asV );





extern double integUnderUnunifTr(
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
              double * rhoInPrevTL_asV );





extern int quadrAngleType(
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
              double * thiVsecT );                    //   -  Third vertex of second triangle.





extern double spaceVolumeInPrevTL(
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
              double * rhoInPrevTL_asV );
