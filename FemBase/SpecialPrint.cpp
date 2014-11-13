


#include <iostream>

using namespace std;

#include<stdio.h>

#include<stdlib.h>

#include<string.h>

#include "OriginPrint.h"                              //   -  Data printing to file in special form for software "TecPlot".



bool printVectorBy3D(
              char *fileName,
              //
              int n,                                  //   -  характеристика сетки.
              int numOfCurrTimeIter,                  //   -  На какой шаг по времени относятся эти данные.
              int numOfSolOrd,                        //   -  Solution order which we want to get.
              //
              double numOfTimeSteps,                  //   -  Общее число шагов по времени.
              //
              double *masOx,
              double *VectorOfData,
              //
              int dimOfVect )
{

char strOfName[50];

char str1[3] = "N=", str2[7];

char str3[3] = "k=", str4[7];

char str5[5] = ".dat", str6[2];




//   Создание строки c названием файла куда будем записывать информацию.

strcpy( strOfName, fileName );



//itoa( numOfSolOrd, str6, 10 );
sprintf(str6,"%d",numOfSolOrd);


strcat( strOfName, str6 );


strcat( strOfName, str1 );

//itoa(   n, str2, 10 );
sprintf(str2,"%d",n);

strcat( strOfName, str2 );


strcat( strOfName, str3 );

//itoa( numOfCurrTimeIter, str4, 10 );
sprintf(str4,"%d",numOfCurrTimeIter);

strcat( strOfName, str4 );


strcat( strOfName, str5 );




double *mas_Empty = new double [ 1 ];

mas_Empty[0] = static_cast<double>(static_cast<int>(numOfCurrTimeIter *1000. /numOfTimeSteps)) /10;


//	Вывод данных.

print_TecPlot_3D(	strOfName,
					//
					masOx,	dimOfVect,
					//
					mas_Empty, 1,
					VectorOfData
					);


delete[] mas_Empty;


return true;
}







bool printSurface_asV(
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
              double *matData_asV )                   //   -  Matrix of data.
{



char strOfName[50];

char str1[4] = "Nx=", str2[7];

char str3[3] = "k=", str4[7];

char str5[5] = ".dat", str6[2];





//   Let's make a string name of the file.

strcpy( strOfName, fileName );


//itoa( numOfSolOrd, str6, 10 );
sprintf(str6,"%d",numOfSolOrd);
strcat( strOfName, str6 );


strcat( strOfName, str1 );

//itoa(   numOfOXSt, str2, 10 );
sprintf(str2,"%d",numOfOXSt);

strcat( strOfName, str2 );


strcat( strOfName, str3 );

//itoa(    iCurrTSt, str4, 10 );
sprintf(str4,"%d",iCurrTSt);

strcat( strOfName, str4 );


strcat( strOfName, str5 );



//	  Data printing.

printSurfaceByMatrix_asV(
              strOfName,
              //
              masOX,                                  //   -  Massive of abscissa grid points. Dimension = numOfOXSt +1.
              numOfOXSt,                              //   -  Number of abscissa grid steps.
              //
              masOY,                                  //   -  Massive of ordinate grid points. Dimension = numOfOYSt +1.
              numOfOYSt,                              //   -  Number of ordinate grid steps.
              //
              matData_asV );                          //   -  Matrix of data.



return true;
}
