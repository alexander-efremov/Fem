


#include <iostream>

using namespace std;

#include <math.h>
#include <string.h>
#include <stdio.h>

#include "InitAndBoundaryData.h"                      //   -  Initial and Boundary data.

#include "SpecialPrint.h"

#include "SpaceVolume.h"

//   #include "RichExtrap.h"

bool solByEqualVolumes(
                       double par_a, //   -  Item of first initial or boundary data parameter.
                       double par_b, //   -  Item of second parameter from "u_funcion".
                       //
                       double lbDom, //   -  Left and right boundaries of rectangular domain.
                       double rbDom,
                       //
                       double bbDom, //   -  Botton and upper boundaries of rectangular domain.
                       double ubDom,
                       //
                       double tau, //   -  Time step.
                       int numOfTSt, //   -  A number of time steps.
                       //
                       double *masOX, //   -  Massive of OX points. Dimension = numOfOXSt +1.
                       int numOfOXSt, //   -  Number of OX steps.
                       //
                       double *masOY, //   -  Massive of OY points. Dimension = numOfOYSt +1.
                       int numOfOYSt, //   -  Number of OY steps.
                       //
                       int numOfSolOrd, //   -  For print only. Solution order which we want to get.
                       //
                       double *rhoInCurrTL_asV) //   -  Rho (solution) in Current Time Layer which we will compute.
{



    double *rhoInPrevTL_asV; //   -  Rho (solution) in Previous Time Layer which we have computed.

    double spVolInPrevTL; //   -  Space volume of rho in previous time layer.

    bool isSpVolComp; //   -  Is space volume computed?

    double RPInCurrTL; //   -  Right part in current time layer.



    int iCurrTL; //   -  Index of current time Layer IN WHICH we computing solution.

    int iOfOXN; //   -  Index of current OX node.

    int iOfOYN; //   -  Index of current OY node.

    int iOfThr; //   -  Through OXY plane index.



    double iForPrint = 0.25 * numOfTSt; //   -  multiple of indices for printing to files.

    double anSol; //   -  Analitical solution in one grid node.

    double normOfErr; //   -  Norm of absolute error in one time level.

    double buf_D; //   -  Buffer. Only.

    int j, k; //   -  Indices.

    bool bul;





    //   New memory.

    /*
    rhoInPrevTL = new double * [ numOfOXSt +1 ];

    for( j=0; j< numOfOXSt +1; j++ )
    {
       rhoInPrevTL[j] = new double [ numOfOYSt +1 ];
    }
     */

    rhoInPrevTL_asV = new double [ (numOfOXSt + 1) * (numOfOYSt + 1) ];





    //   Initial data of rho.

    for (k = 0; k < numOfOYSt + 1; k++)
    {
        for (j = 0; j < numOfOXSt + 1; j++)
        {
            rhoInPrevTL_asV[ (numOfOXSt + 1) * k + j ] = initDataOfSol(par_a, lbDom, rbDom, bbDom, ubDom, j, masOX, k, masOY);
        }
    }



    //   Let's do it. It's easy.
    //numOfTSt = 1;
    //printf("tl setted to 1\n");
    for (iCurrTL = 1; iCurrTL < numOfTSt + 1; iCurrTL++)
    {
        //   "iCurrTL" - index of current time layer in which we computing the solution.

        //   I want to see what is computing :)

        //cout<<"\r \t \t \t \t \t \t \t \t \r";

        //cout << "SchOrd = " << numOfSolOrd << ", Nx = " << numOfOXSt;

        //cout << ", indexOfCurrTL = " << iCurrTL << " ( " << numOfTSt << " ) " << flush;



        //   If we know solution on the boundary we can use it.

        for (iOfOXN = 0; iOfOXN < numOfOXSt + 1; iOfOXN++)
        {
            //   Botton boundary.

            //   rhoInCurrTL[ iOfOXN ][ 0 ]  =   bottonBound( par_a,   lbDom, rbDom,   bbDom, ubDom,   tau*iCurrTL, masOX[ iOfOXN ] );

            rhoInCurrTL_asV[ iOfOXN ] = bottonBound(par_a, lbDom, rbDom, bbDom, ubDom, tau*iCurrTL, masOX[ iOfOXN ]);

            //   Upper boundary.

            //   rhoInCurrTL[ iOfOXN ][ numOfOYSt ]  =   upperBound( par_a,   lbDom, rbDom,   bbDom, ubDom,   tau*iCurrTL, masOX[ iOfOXN ] );

            rhoInCurrTL_asV[ (numOfOXSt + 1) * numOfOYSt + iOfOXN ] = upperBound(par_a, lbDom, rbDom, bbDom, ubDom, tau*iCurrTL, masOX[ iOfOXN ]);
        }


        for (iOfOYN = 0; iOfOYN < numOfOYSt + 1; iOfOYN++)
        {
            //   Left boundary.

            //   rhoInCurrTL[ 0 ][ iOfOYN ]  =  leftBound( par_a,   lbDom, rbDom,   bbDom, ubDom,   tau*iCurrTL, masOY[ iOfOYN ] );

            rhoInCurrTL_asV[ (numOfOXSt + 1) * iOfOYN ] = leftBound(par_a, lbDom, rbDom, bbDom, ubDom, tau*iCurrTL, masOY[ iOfOYN ]);

            //   Right boundary.

            //   rhoInCurrTL[ numOfOXSt ][ iOfOYN ]  =  rightBound( par_a,   lbDom, rbDom,   bbDom, ubDom,   tau*iCurrTL, masOY[ iOfOYN ] );

            rhoInCurrTL_asV[ (numOfOXSt + 1) * iOfOYN + numOfOXSt ] = rightBound(par_a, lbDom, rbDom, bbDom, ubDom, tau*iCurrTL, masOY[ iOfOYN ]);
        }



        //   Enumeration from first unknown element to last one.

        for (iOfOYN = 1; iOfOYN < numOfOYSt; iOfOYN++)
        {
            for (iOfOXN = 1; iOfOXN < numOfOXSt; iOfOXN++)
            {
                //   And now the mo-o-o-o-st intresting....!  The volume computation.

                spVolInPrevTL = spaceVolumeInPrevTL(
                                                    par_a, par_b, //   -  Items of parameters.
                                                    //
                                                    lbDom, rbDom, //   -  Left and right boundaries of rectangular domain.
                                                    //
                                                    bbDom, ubDom, //   -  Botton and upper boundaries of rectangular domain.
                                                    //
                                                    tau, iCurrTL, //   -  Time data. Necessary for velocity.
                                                    //
                                                    iOfOXN, masOX, numOfOXSt, //   -  OX data.
                                                    //
                                                    iOfOYN, masOY, numOfOYSt, //   -  OY data.
                                                    //
                                                    isSpVolComp, //   -  Is space volume computed?
                                                    rhoInPrevTL_asV);

                /*  if ((numOfOXSt +1)*iOfOYN + iOfOXN == 12)
                  {
                      printf("\na = %f\n", par_a);
                      printf("b = %f\n", par_b);
                      printf("lbDom = %f\n", lbDom);
                      printf("rbDom = %f\n", rbDom);
                      printf("bbDom = %f\n", bbDom);
                      printf("ubDom = %f\n", ubDom);
                      printf("tau = %f\n", tau);
                      printf("iCurTl = %d\n", iCurrTL);
                       printf("iOfOXN = %d\n", iOfOXN);
                        printf("iOfOYN = %d\n", iOfOYN);
                         printf("numOfOXSt = %d\n", numOfOXSt);
                         printf("numOfOYSt = %d\n", numOfOYSt);
                      printf("spVolInPrevTL = %f\n", spVolInPrevTL);
                  }*/
                //   Item of spVolInPrevTL should be divided by hx and hy;

                //   Let's divide by hx of ununiform grid.

                buf_D = (masOX[iOfOXN + 1] - masOX[iOfOXN - 1]) / 2.;

                spVolInPrevTL = spVolInPrevTL / buf_D;

                //   Let's divide by hy of ununiform grid.

                buf_D = (masOY[iOfOYN + 1] - masOY[iOfOYN - 1]) / 2.;

                spVolInPrevTL = spVolInPrevTL / buf_D;



                //   2. Let's compute new item of rho :)

                RPInCurrTL = f_function(//   -  It's an initial right part of differential equation.
                                        par_a, par_b,
                                        //
                                        lbDom, rbDom,
                                        bbDom, ubDom,
                                        //
                                        tau,
                                        iCurrTL, //   -  Index of current time layer.
                                        //
                                        iOfOXN, //   -  Index of current OX node.
                                        masOX, //   -  Massive of OX steps. Dimension = numOfOXSt +1.
                                        numOfOXSt, //   -  Number of OX steps.
                                        //
                                        iOfOYN, //   -  Index of current OY node.
                                        masOY, //   -  Massive of OY steps. Dimension = numOfOYSt +1.
                                        numOfOYSt); //   -  Number of OY steps.



                rhoInCurrTL_asV[ (numOfOXSt + 1) * iOfOYN + iOfOXN ] = spVolInPrevTL;

                //   rhoInCurrTL[ iOfOXN ][ iOfOYN ]  =  spVolInPrevTL;

                /*
                         if( (iOfOXSt == 0)  ||  (iOfOXSt == numOfOXSt) ){   rhoInCurrTSt[ iOfOXSt ][ iOfOYSt ] *= 2.;   }

                         //   and

                         if( (iOfOYSt == 0)  ||  (iOfOYSt == numOfOYSt) ){   rhoInCurrTSt[ iOfOXSt ][ iOfOYSt ] *= 2.;   }
                 */

                //   rhoInCurrTL[ iOfOXN ][ iOfOYN ] +=  tau * RPInCurrTL;

                rhoInCurrTL_asV[ (numOfOXSt + 1) * iOfOYN + iOfOXN ] += tau * RPInCurrTL;
                //if ((numOfOXSt +1)*iOfOYN + iOfOXN == 12)
                //    {
                //      printf("spVolInPrevTL = %f\n", rhoInCurrTL_asV[ (numOfOXSt +1)*iOfOYN + iOfOXN ]);
                //}
            }
        }



        //   Printing if it's necessary.

        if (((int) ((int) (iCurrTL / iForPrint) * iForPrint) == iCurrTL))
        {

            for (iOfOYN = 0; iOfOYN < numOfOYSt + 1; iOfOYN++)
            {
                for (iOfOXN = 0; iOfOXN < numOfOXSt + 1; iOfOXN++)
                {
                    anSol = analytSolut(par_a, lbDom, rbDom, bbDom, ubDom, iCurrTL *tau, masOX[ iOfOXN ], masOY[ iOfOYN ]);

                    //   Let's use "rhoInPrevTSt" as absolute error.

                    rhoInPrevTL_asV[ (numOfOXSt + 1) * iOfOYN + iOfOXN ] = anSol - rhoInPrevTL_asV[ (numOfOXSt + 1) * iOfOYN + iOfOXN ];
                }
            }



            //   Solution visualization.

//            bul = printSurface_asV(
//                                   "rho", //   -  char *fileName,
//                                   //
//                                   numOfOXSt, //   -  Number of OX steps. Grid parameter.
//                                   iCurrTL, //   -  Index of current time layer IN WHICH we printing solution.
//                                   numOfSolOrd, //   -  Solution order which we want to print.
//                                   //
//                                   numOfTSt, //   -  A number of time steps.
//                                   //
//                                   masOX, //   -  Massive of OX points. Dimension = numOfOXSt +1.
//                                   numOfOXSt, //   -  Number of OX steps.
//                                   //
//                                   masOY, //   -  Massive of OY points. Dimension = numOfOYSt +1.
//                                   numOfOYSt, //   -  Number of OY steps.
//                                   //
//                                   rhoInCurrTL_asV); //   -  Matrix of data.
//
//            //   Error visualization.
//
//            bul = printSurface_asV(
//                                   "rhoErr", //   -  char *fileName,
//                                   //
//                                   numOfOXSt, //   -  Number of OX steps. Grid parameter.
//                                   iCurrTL, //   -  Index of current time layer IN WHICH we printing solution.
//                                   numOfSolOrd, //   -  Solution order which we want to print.
//                                   //
//                                   numOfTSt, //   -  A number of time steps.
//                                   //
//                                   masOX, //   -  Massive of OX points. Dimension = numOfOXSt +1.
//                                   numOfOXSt, //   -  Number of OX steps.
//                                   //
//                                   masOY, //   -  Massive of OY points. Dimension = numOfOYSt +1.
//                                   numOfOYSt, //   -  Number of OY steps.
//                                   //
//                                   rhoInPrevTL_asV); //   -  Matrix of data. We use "rhoInPrevTSt" as absolute error.
//

            //   Error item.

            //   normOfErr  =  MaxModItemOfMatr( rhoInPrevTL, numOfOXSt +1, numOfOYSt +1 );
        }



        //   Data update. Only after printing.

        for (iOfThr = 0; iOfThr < (numOfOXSt + 1) * (numOfOYSt + 1); iOfThr++)
        {
            rhoInPrevTL_asV[ iOfThr ] = rhoInCurrTL_asV[ iOfThr ];
        }

    } //   end of "for( iCurrTL = 1; iCurrTL< numOfTSt +1; iCurrTL++ )".





    //   Memory clearing.
    /*
    for( j=0; j< numOfOXSt +1; j++ )
    {
       delete rhoInPrevTL[j];
    }

    delete rhoInPrevTL;
     */

    delete rhoInPrevTL_asV;



    return true;
}

void print_matrix(int n, int m, double *a, int precision = 8)
{
    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < m; ++j)
        {
            int k = i * n + j;
            switch (precision)
            {
            case 1:
                printf("%.1f ", a[k]);
                break;
            case 2:
                printf("%.2f ", a[k]);
                break;
            case 3:
                printf("%.3f ", a[k]);
                break;
            case 4:
                printf("%.4f ", a[k]);
                break;
            case 5:
                printf("%.5f ", a[k]);
                break;
            case 6:
                printf("%.6f ", a[k]);
                break;
            case 7:
                printf("%.7f ", a[k]);
                break;
            case 8:
                printf("%.8f ", a[k]);
                break;
            }
        }
        printf("\n");
    }
}

double* solByEqualVolWithVarStepPlusPrint1(
                                           double par_a, //   -  Item of first initial or boundary data parameter.
                                           double par_b, //   -  Item of second parameter from "u_funcion".
                                           //
                                           double lbDom, //   -  Left and right boundaries of rectangular domain.
                                           double rbDom,
                                           //
                                           double bbDom, //   -  Botton and upper boundaries of rectangular domain.
                                           double ubDom,
                                           //
                                           double tau, //   -  Time step.
                                           int numOfTSt, //   -  A number of time steps.

                                           int numOfOXSt, //   -  Number of OX steps.

                                           int numOfOYSt,
                                           int numOfGrStepLayer, double &norm) //   -  How many computations with different grid steps we want to make.
{
    int indByNumOfGridSteps; //   -  Index by computation with different grid step.


    double varTau = tau; //   -  Variable time step.

    int varNumOfTSt = numOfTSt; //   -  Variable number of time steps.


    double * varMasOX; //   -  Variable massive of OX nodes.

    int varNumOfOXSt = numOfOXSt; //   -  Variable number of OX steps.


    double * varMasOY; //   -  Variable massive of OY nodes.

    int varNumOfOYSt = numOfOYSt; //   -  Variable number of OY steps.



    double *rhoInCurrTL_asV; //   -  Rho (solution) in current (last one) time level. Matrix as vector.

    double anSol; //   -  Analytical solution i one point.
    double* result = NULL;

    double buf_D;


    int j, k;

    bool bul;


    indByNumOfGridSteps = numOfGrStepLayer;


    //   New time step.


    varTau = tau ;

    varNumOfTSt = numOfTSt ;



    //   New absciss grid.

    varMasOX = new double [ varNumOfOXSt + 1 ];

    buf_D = (rbDom - lbDom) / varNumOfOXSt;

    for (j = 0; j < varNumOfOXSt + 1; j++)
    {
        varMasOX[j] = lbDom + ((double) j) * buf_D;
    }

    varMasOY = new double [ varNumOfOYSt + 1 ];

    buf_D = (ubDom - bbDom) / varNumOfOYSt;

    for (k = 0; k < varNumOfOYSt + 1; k++)
    {
        varMasOY[k] = bbDom + ((double) k) * buf_D;
    }




    rhoInCurrTL_asV = new double [ (varNumOfOXSt + 1) * (varNumOfOYSt + 1) ];


    //   Computation of solution.

    bul = solByEqualVolumes(
                            par_a, par_b,
                            //
                            lbDom, rbDom,
                            //
                            bbDom, ubDom,
                            //
                            varTau, //   -  Time step.
                            varNumOfTSt, //   -  A number of time steps.
                            //
                            varMasOX, //   -  Massive of abscissa grid points. Dimension = varNumOfOxGrSt +1.
                            varNumOfOXSt, //   -  Variable number of abscissa grid steps.
                            //
                            varMasOY, //   -  Massive of ordinate grid points. Dimension = varNumOfOyGrSt +1.
                            varNumOfOYSt, //   -  Variable number of ordinate grid steps.
                            //
                            0, //   -  For print only. Solution order which we want to get.
                            //
                            rhoInCurrTL_asV); //   -  Rho (solution) in Current (Last) Time Level.

    //  print_matrix(varNumOfOXSt + 1, varNumOfOYSt + 1, rhoInCurrTL_asV);

    // copy result to tmp
    result = new double[(varNumOfOXSt + 1)*(varNumOfOYSt + 1)];
    memcpy(result, rhoInCurrTL_asV, (varNumOfOXSt + 1)*(varNumOfOYSt + 1) * sizeof (double));
    //   :) Now we know solution! It's great. What we need more?

    //   May be error distribution and it's maximum? It's easy :)

    buf_D = varNumOfTSt * varTau;

    for (k = 0; k < varNumOfOYSt + 1; k++)
    {
        for (j = 0; j < varNumOfOXSt + 1; j++)
        {
            anSol = analytSolut(par_a, lbDom, rbDom, bbDom, ubDom, buf_D, varMasOX[j], varMasOY[k]);
            rhoInCurrTL_asV[ (varNumOfOXSt + 1) * k + j ] = anSol - rhoInCurrTL_asV[ (varNumOfOXSt + 1) * k + j ];
        }
    }
    norm = norm_at_L_1(
                        varMasOX, //   -  Massive of OX grid nodes. Dimension = dimOX.
                        varNumOfOXSt + 1,
                        varMasOY, //   -  Massive of OY grid nodes. Dimension = dimOY.
                        varNumOfOYSt + 1,
                        rhoInCurrTL_asV);
    delete varMasOX;
    delete varMasOY;
    delete rhoInCurrTL_asV;


    return result;
}

double* solByEqualVolWithVarStepPlusPrint(
                                          double par_a, //   -  Item of first initial or boundary data parameter.
                                          double par_b, //   -  Item of second parameter from "u_funcion".
                                          //
                                          double lbDom, //   -  Left and right boundaries of rectangular domain.
                                          double rbDom,
                                          //
                                          double bbDom, //   -  Botton and upper boundaries of rectangular domain.
                                          double ubDom,
                                          //
                                          double tau, //   -  Time step.
                                          int numOfTSt, //   -  A number of time steps.

                                          int numOfOXSt, //   -  Number of OX steps.

                                          int numOfOYSt, //   -  Number of OY steps.
                                          //
                                          bool isTimeStShBeChan, //   -  Is time step shoulb be change?
                                          bool isGridStShBeChan, //   -  Is grid step shoulb be change?
                                          //
                                          int numOfGrStepLayer, double *norm) //   -  How many computations with different grid steps we want to make.
{
    int indByNumOfGridSteps; //   -  Index by computation with different grid step.

    double *maxModAbsErr; //   -  Maximum of absolute error on corresponding grid.

    double *masOfGrStepItem; //   -  Massive of grid steps (parameters).

    double *ordOfErr; //   -  Order of error computed on two different grids.



    double varTau = tau; //   -  Variable time step.

    int varNumOfTSt = numOfTSt; //   -  Variable number of time steps.


    double * varMasOX; //   -  Variable massive of OX nodes.

    int varNumOfOXSt = numOfOXSt; //   -  Variable number of OX steps.


    double * varMasOY; //   -  Variable massive of OY nodes.

    int varNumOfOYSt = numOfOYSt; //   -  Variable number of OY steps.



    double *rhoInCurrTL_asV; //   -  Rho (solution) in current (last one) time level. Matrix as vector.

    double anSol; //   -  Analytical solution i one point.

    //   double **absErr; - To save memory we can use "rhoInCurrTL".  //   -  Absolute error in time for error solution estamation.

    double* result = NULL;

    double buf_D;

    //   int buf_I;

    int j, k;

    bool bul;







    //   Memory declaraion.

    maxModAbsErr = new double [ numOfGrStepLayer ];

    masOfGrStepItem = new double [ numOfGrStepLayer - 1 ];

    ordOfErr = new double [ numOfGrStepLayer - 1 ];





    for (indByNumOfGridSteps = 0; indByNumOfGridSteps < numOfGrStepLayer; indByNumOfGridSteps++)
    {
        //   New time step.

        if ((isTimeStShBeChan == true) || (indByNumOfGridSteps == 0))
        {
            varTau = tau / pow(2, indByNumOfGridSteps);

            varNumOfTSt = numOfTSt * pow(2, indByNumOfGridSteps);
        }



        if ((isGridStShBeChan == true) || (indByNumOfGridSteps == 0))
        {
            //   New absciss grid steps.

            varNumOfOXSt = numOfOXSt * pow(2, indByNumOfGridSteps);

            //   New absciss grid.

            varMasOX = new double [ varNumOfOXSt + 1 ];

            buf_D = (rbDom - lbDom) / varNumOfOXSt;

            for (j = 0; j < varNumOfOXSt + 1; j++)
            {
                varMasOX[j] = lbDom + ((double) j) * buf_D;
            }



            //   New ordinate grid steps.

            varNumOfOYSt = numOfOYSt * pow(2, indByNumOfGridSteps);

            //   New absciss grid.

            varMasOY = new double [ varNumOfOYSt + 1 ];

            buf_D = (ubDom - bbDom) / varNumOfOYSt;

            for (k = 0; k < varNumOfOYSt + 1; k++)
            {
                varMasOY[k] = bbDom + ((double) k) * buf_D;
            }



            //   New matrices for solution and for error.
            /*
            rhoInCurrTL = new double * [varNumOfOXSt +1];

            for( j=0; j< varNumOfOXSt +1; j++ )
            {
               rhoInCurrTL[j] = new double [varNumOfOYSt +1];
            }
             */

            rhoInCurrTL_asV = new double [ (varNumOfOXSt + 1) * (varNumOfOYSt + 1) ];
        }

        //   Computation of solution.

        bul = solByEqualVolumes(
                                par_a, par_b,
                                //
                                lbDom, rbDom,
                                //
                                bbDom, ubDom,
                                //
                                varTau, //   -  Time step.
                                varNumOfTSt, //   -  A number of time steps.
                                //
                                varMasOX, //   -  Massive of abscissa grid points. Dimension = varNumOfOxGrSt +1.
                                varNumOfOXSt, //   -  Variable number of abscissa grid steps.
                                //
                                varMasOY, //   -  Massive of ordinate grid points. Dimension = varNumOfOyGrSt +1.
                                varNumOfOYSt, //   -  Variable number of ordinate grid steps.
                                //
                                0, //   -  For print only. Solution order which we want to get.
                                //
                                rhoInCurrTL_asV); //   -  Rho (solution) in Current (Last) Time Level.

        //  print_matrix(varNumOfOXSt + 1, varNumOfOYSt + 1, rhoInCurrTL_asV);

        // copy result to tmp
        result = new double[(varNumOfOXSt + 1)*(varNumOfOYSt + 1)];
        memcpy(result, rhoInCurrTL_asV, (varNumOfOXSt + 1)*(varNumOfOYSt + 1) * sizeof (double));
        //   :) Now we know solution! It's great. What we need more?

        //   May be error distribution and it's maximum? It's easy :)

        buf_D = varNumOfTSt * varTau;

        for (k = 0; k < varNumOfOYSt + 1; k++)
        {
            for (j = 0; j < varNumOfOXSt + 1; j++)
            {
                anSol = analytSolut(par_a, lbDom, rbDom, bbDom, ubDom, buf_D, varMasOX[j], varMasOY[k]);
                rhoInCurrTL_asV[ (varNumOfOXSt + 1) * k + j ] = anSol - rhoInCurrTL_asV[ (varNumOfOXSt + 1) * k + j ];
            }
        }
        // printf("\nError\n");
        // print_matrix(varNumOfOXSt + 1, varNumOfOYSt + 1, rhoInCurrTL_asV);
        //   maxModAbsErr[ indByNumOfGridSteps ]  =  MaxModItemOfMatr( rhoInCurrTL, varNumOfOXSt +1, varNumOfOYSt +1);



        maxModAbsErr[ indByNumOfGridSteps ] = norm_at_L_1(
                                                          varMasOX, varNumOfOXSt + 1,
                                                          varMasOY, varNumOfOYSt + 1, rhoInCurrTL_asV);
        *norm = maxModAbsErr[indByNumOfGridSteps];
        //  printf("\nindByNumOfGridSteps = %d\n Norm %f\n", indByNumOfGridSteps, maxModAbsErr[ indByNumOfGridSteps ]);

        //   Let's compute order of numerical solution.

        if (indByNumOfGridSteps > 0)
        {
            masOfGrStepItem[ indByNumOfGridSteps - 1] = indByNumOfGridSteps - 0.5;


            buf_D = maxModAbsErr[ indByNumOfGridSteps - 1] / maxModAbsErr[ indByNumOfGridSteps ];

            ordOfErr[ indByNumOfGridSteps - 1] = log(buf_D) / log(2.);
        }

        //   Memory clearing.

        if (isGridStShBeChan == true)
        {
            delete varMasOX;
            delete varMasOY;
            delete rhoInCurrTL_asV;
        }
    }

    if (numOfGrStepLayer > 1)
    {
        //   Let's write data to file.

//        bul = printVectorBy3D(
//                              "ordOfRhoErr", //   -  char *fileName,
//                              //
//                              numOfOXSt, //   -  �������������� �����.
//                              numOfGrStepLayer, //   -  � �������� �����.
//                              1. / 3., //   -  ��������� ���� �� ������� � ���� �� ������������.
//                              numOfGrStepLayer, //   -  ����� ����� ����� �� ����.
//                              //
//                              masOfGrStepItem, //   -  double *masOx,
//                              ordOfErr, //   -  double *VectorOfData,
//                              //
//                              numOfGrStepLayer - 1); //   -  dimOfVect.
    }

    if (isGridStShBeChan == false)
    {
        delete varMasOX;
        delete varMasOY;
        delete rhoInCurrTL_asV;
    }
    delete maxModAbsErr;
    delete masOfGrStepItem;
    delete ordOfErr;
    return result;
}

