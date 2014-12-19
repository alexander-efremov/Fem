


extern void print_TecPlot_3D (
            char *filename,
            double *mas_0x,
            int len_x,
            double *mas_0y,
            int len_y,
            double *mas_0z );





void printSurfaceByMatrix_asV(
              char *filename,
              //
              double *masOX,                          //   -  Massive of abscissa grid points. Dimension = numOfOXSt +1.
              int numOfOXSt,                          //   -  Number of abscissa grid steps.
              //
              double *masOY,                          //   -  Massive of ordinate grid points. Dimension = numOfOYSt +1.
              int numOfOYSt,                          //   -  Number of ordinate grid steps.
              //
              double *matData_asV );                  //   -  Matrix of data.
