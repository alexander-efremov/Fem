/*#include <iostream>

using namespace std;


#include "Const.h"                                    //   -  Helpful used consts.

#include "InitAndBoundaryData.h"                      //   -  Initial and Boundary data.

#include "OriginPrint.h"                              //   -  Printing to file for.

#include "SpecialPrint.h"                             //   -  Printing to file for.

#include "LowOrdOper.h"                               //   -  Operator on low order accurancy.

//   #include "OneCellInteg.h"                               //   -  Richardson extrapolation.





int main()
{



int numOfGrStepLayer = 1;                            //   -  How many computations with different grid steps we want to make.

double buf_D;

bool bul;




//   Computation of necessary initial global variable.

initCompOfGlVar();





buf_D = itemOfInteg_2SpecType(
              1.,   //   Py,
              3.,   //   Qy,
              //
              2.,   //   alpha,
              //
              2.;   //   a,
              9.,   //   b,
              4. ); //   betta



solByEqualVolWithVarStepPlusPrint(
              C_par_a,                                //   -  Item of first parameter.
              C_par_b,                                //   -  Item of second parameter from "u_funcion".
              //
              C_lbDom, C_rbDom,                       //   -  Left and right boundaries of rectangular domain.
              //
              C_bbDom, C_ubDom,                       //   -  Botton and upper boundaries of rectangular domain.
              //
              C_tau,                                  //   -  Time step.
              C_numOfTSt,                             //   -  A number of time steps.
              //
                                             //   -  Massive of OX points. Dimension = C_numOfOXSt +1.
              C_numOfOXSt,                            //   -  Number of OX steps.
              //
                                           //   -  Massive of OY points. Dimension = C_numOfOXSt +1.
              C_numOfOYSt,                            //   -  Number of OY steps.
              //
              true,                                   //   -  Is time step shoulb be change?
              true,                                   //   -  Is grid step shoulb be change?
              //
              numOfGrStepLayer );                     //   -  How many computations with different grid steps we want to make.





bul = RichExtrPlusPrint(
              C_par_a,                                //   -  Item of left and right setback (parameter "a" in test).
              C_par_b,                                //   -  Item of second parameter from "u_funcion".
			     //
              C_numOfGridSteps,                       //   -  Number of segments of domain [0; 1] irregular decomposition.
              //
              C_tau,                                  //   -  Time step.
              C_numOfTimeSteps,                       //   -  A number of time steps.
              //
              2 );                                    //   -  "numOfSolOrd"  -  ������� �������, ������� �� ����� ��������.





//   Memory cleaning.

memClean();

return 0;
}*/
