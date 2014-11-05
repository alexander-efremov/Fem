


#include <iostream>

using namespace std;

#include<stdio.h>

#include<stdlib.h>



void print_TecPlot_3D (
            char *filename,
            double *mas_0x,
            int len_x,
            double *mas_0y,
            int len_y,
            double *mas_0z )
{
/**********************************************************************************

filename  - ��� ����� ���� ����� ������ ����������.
*mas_0x   - ������ �������� �� ��� �������.
len_x   - (length) ������ ������� �� ��� �������.
*mas_0y   - ������ �������� �� ��� �������.
len_y   - (length) ������ ������� �� ��� �������.
*mas_0z   - ������ �������� �� ��� z, �������� � ���� ������,
        ����������� (�� ���������) len_x*len_y.

������ mas_0z ������ ���� ���������� ��������� �������:


    y[0]    y[1]    y[2]    y[3]
    |     |     |     |
    |     |     |     |
    |     |     |     |
x[0]-------------------------------------------------
    |z[0] |z[5] |z[10]|
    |     |     |     |
    |     |     |     |
x[1]-------------------------------------------------
    |z[1] |z[6] |z[11]|
    |     |     |     |
    |     |     |     |
x[2]-------------------------------------------------
    |z[2] |z[7] |...  |
    |     |     |     |
    |     |     |     |
x[3]-------------------------------------------------
    |z[3] |z[8] |     |
    |     |     |     |
    |     |     |     |
x[4]-------------------------------------------------
    |z[4] |z[9] |     |
    |     |     |     |
    |     |     |     |



************************************************************************************/

FILE *pfile;
char answ;

// Cheking of existance such file.
/*
pfile = fopen( filename, "r");

if( pfile != NULL )
  {
  cout<<"\nBe carefully, the file ";
  cout<<filename;
  cout<<" is already exist and all datas at one will be lose";

question:

  cout<<"\nDo you want to continue, 'y' or 'n'? "<<endl;
  cin>>answ;
  if ( answ == 'n')
    {
    cout<<"Nothig not write to the file ";
    cout<<filename;
    fclose(pfile);
    return ;
    }
  if ( answ == 'y') cout<<"Ok, all datas at one will be lose! ";

  if ( answ != 'y') goto question;

  fclose(pfile);
  }
*/
// ������ ����� � ����

pfile = fopen( filename, "w");
if( pfile == NULL ) { cout<<"\nError can not creat-open file ", filename; exit(1); }

// ��������� � ����� ������� ���������� ��� ��������� TecPlot, ������� � ����� ������������ ���������!

/*
TITLE     = 'DEM DATA | DEM DATA | DEM DATA | DEM DATA'
VARIABLES = 'X' 'Y' 'E'
ZONE T='SubZone'
I=51, J=5, K=1, ZONETYPE=Ordered
DATAPACKING=POINT
DT=(SINGLE SINGLE SINGLE )
*/

fprintf(pfile,  "TITLE = 'DEM DATA | DEM DATA | DEM DATA | DEM DATA'" );
fprintf(pfile,  "\nVARIABLES = 'X' 'Y' 'E' " );
fprintf(pfile,  "\nZONE T='SubZone'" );
fprintf(pfile,  "\nI=%d J=%d K=%d ZONETYPE=Ordered", len_y, len_x, 1 );
fprintf(pfile,  "\nDATAPACKING=POINT" );
fprintf(pfile,  "\nDT=(SINGLE SINGLE SINGLE )" );

// ����� ��������, ������ ���������� ���� ������, ��� ������� ��������� �������:

/*

x[0]    y[0]    z[0]
x[0]    y[1]    z[5]
x[0]    y[2]    z[6]
x[0]    y[3]    z[7]
x[1]    y[0]    z[1]
x[1]    y[1]    z[6]
x[1]    y[2]    z[11]
x[1]    y[3]    z[16]
...
x[i]    y[j]    z[ 5*j + i ]
...


    y[0]    y[1]    y[2]    y[3]
    |     |     |     |
    |     |     |     |
    |     |     |     |
x[0]-------------------------------------------------
    |z[0]   |z[5]   |z[10]    |
    |     |     |     |
    |     |     |     |
x[1]-------------------------------------------------
    |z[1]   |z[6]   |z[11]    |
    |     |     |     |
    |     |     |     |
x[2]-------------------------------------------------
    |z[2]   |z[7]   |...    |
    |     |     |     |
    |     |     |     |
x[3]-------------------------------------------------
    |z[3]   |z[8]   |     |
    |     |     |     |
    |     |     |     |
x[4]-------------------------------------------------
    |z[4]   |z[9]   |     |
    |     |     |     |
    |     |     |     |

*/


for (int i=0; i<len_x; i++ )
  for (int j=0; j<len_y; j++ )
    fprintf(pfile, "\n%-30.20g  %-30.20g %-30.20g", mas_0x[i], mas_0y[j], mas_0z[ j*len_x + i ] );

fclose(pfile);
// cout<<"\nAll dates are written to the file: "; cout<<filename;
return ;
}







void printSurfaceByMatrix_asV(
              char *filename,
              //
              double *masOX,                          //   -  Massive of abscissa grid points. Dimension = numOfOXSt +1.
              int numOfOXSt,                          //   -  Number of abscissa grid steps.
              //
              double *masOY,                          //   -  Massive of ordinate grid points. Dimension = numOfOYSt +1.
              int numOfOYSt,                          //   -  Number of ordinate grid steps.
              //
              double *matData_asV )                   //   -  Matrix of data.
{
/*


    y[0]    y[1]    y[2]    y[3]
    |     |     |     |
    |     |     |     |
    |     |     |     |
x[0]-------------------------------------------------
    |z[0] |z[5] |z[10]|
    |     |     |     |
    |     |     |     |
x[1]-------------------------------------------------
    |z[1] |z[6] |z[11]|
    |     |     |     |
    |     |     |     |
x[2]-------------------------------------------------
    |z[2] |z[7] |...  |
    |     |     |     |
    |     |     |     |
x[3]-------------------------------------------------
    |z[3] |z[8] |     |
    |     |     |     |
    |     |     |     |
x[4]-------------------------------------------------
    |z[4] |z[9] |     |
    |     |     |     |
    |     |     |     |



************************************************************************************/

FILE *pfile;

char answ;



pfile = fopen( filename, "w");

if( pfile == NULL ) { cout<<"\nError can not create-open file ", filename; exit(1); }

// ��������� � ����� ������� ���������� ��� ��������� TecPlot, ������� � ����� ������������ ���������!

/*
TITLE     = 'DEM DATA | DEM DATA | DEM DATA | DEM DATA'
VARIABLES = 'X' 'Y' 'E'
ZONE T='SubZone'
I=51, J=5, K=1, ZONETYPE=Ordered
DATAPACKING=POINT
DT=(SINGLE SINGLE SINGLE )
*/

fprintf(pfile,  "TITLE = 'DEM DATA | DEM DATA | DEM DATA | DEM DATA'" );
fprintf(pfile,  "\nVARIABLES = 'X' 'Y' 'E' " );
fprintf(pfile,  "\nZONE T='SubZone'" );
fprintf(pfile,  "\nI=%d J=%d K=%d ZONETYPE=Ordered", numOfOYSt, numOfOXSt, 1 );
fprintf(pfile,  "\nDATAPACKING=POINT" );
fprintf(pfile,  "\nDT=(SINGLE SINGLE SINGLE )" );

// ����� ��������, ������ ���������� ���� ������, ��� ������� ��������� �������:

/*

x[0]    y[0]    z[0]
x[0]    y[1]    z[5]
x[0]    y[2]    z[6]
x[0]    y[3]    z[7]
x[1]    y[0]    z[1]
x[1]    y[1]    z[6]
x[1]    y[2]    z[11]
x[1]    y[3]    z[16]
...
x[i]    y[j]    z[ 5*j + i ]
...


    y[0]    y[1]    y[2]    y[3]
    |     |     |     |
    |     |     |     |
    |     |     |     |
x[0]-------------------------------------------------
    |z[0]   |z[5]   |z[10]    |
    |     |     |     |
    |     |     |     |
x[1]-------------------------------------------------
    |z[1]   |z[6]   |z[11]    |
    |     |     |     |
    |     |     |     |
x[2]-------------------------------------------------
    |z[2]   |z[7]   |...    |
    |     |     |     |
    |     |     |     |
x[3]-------------------------------------------------
    |z[3]   |z[8]   |     |
    |     |     |     |
    |     |     |     |
x[4]-------------------------------------------------
    |z[4]   |z[9]   |     |
    |     |     |     |
    |     |     |     |

*/


for (int i=0; i< numOfOXSt; i++ )
{
  for (int j=0; j< numOfOYSt; j++ )
   {
    //   fprintf(pfile, "\n%-30.20g  %-30.20g %-30.20g", masOX[i], masOY[j], matData[ i ][ j ] );

      fprintf(pfile, "\n%-30.20g  %-30.20g %-30.20g", masOX[i], masOY[j], matData_asV[ (numOfOXSt +1)*j + i ] );
   }
}


fclose(pfile);



return ;
}
