#include <math.h>

//   #include "InitAndBoundaryData.h"                      //   -  Initial and Boundary data.

#include "OneCellInteg.h"


double integOfChan_SLLeftSd(//   -  The domain is Channel with Slant Line on the left side.
	double par_a, //   -  Solution parameter.
	//
	double lbDom, //   -  Left and right boundaries of rectangular domain.
	double rbDom,
	//
	double bbDom, //   -  Botton and upper boundaries of rectangular domain.
	double ubDom,
	//
	double tau,
	int iCurrTL, //   -  Index of current time layer.
	//
	double* bv, int wTrPCI, //   -  Where travel point current (botton vertex) is.
	double* uv, int wTrPNI, //   -  Where travel point next (upper vertex) is.
	//
	int* indCurSqOx, //   -  Index by OX axis where bv and uv are.
	//
	double rb, int* indRB, //   -  Right boundary by Ox. Index by OX axis where rb is.
	//
	int* indCurSqOy, //   -  Index of current square by Oy axis.
	//
	double* masOX, //   -  Massive of OX steps. Dimension = numOfOXSt +1.
	int numOfOXSt, //   -  Number of OX steps.
	//
	double* masOY, //   -  Massive of OY steps. Dimension = numOfOYSt +1.
	int numOfOYSt, //   -  Number of OY steps.
	//
	double* rhoInPrevTL_asV)
{
	double lv[2], mv[2]; //   -  Left and middle vertices.

	int wMvI; //   -  Where middle vertex is.

	int indCurSqOxToCh[2]; //   -  Indices of current square by Ox axis to be changed. Under which we want to integrate.


	double h = masOX[1] - masOX[0];

	double a_SL, b_SL; //   -  Coefficients of slant line: x = a_SL *y  +  b_SL.


	double Gx, Hx; //   -  Left and right boundary for each integration.


	double integ = 0.;

	double buf_D;

	int j;


	//   Let's compute helpful values.

	if (uv[0] <= bv[0])
	{
		lv[0] = uv[0];
		lv[1] = uv[1];

		mv[0] = bv[0];
		mv[1] = bv[1];
		wMvI = wTrPCI;
	}

	if (uv[0] > bv[0])
	{
		lv[0] = bv[0];
		lv[1] = bv[1];

		mv[0] = uv[0];
		mv[1] = uv[1];
		wMvI = wTrPNI;
	}


	if ((fabs(uv[1] - bv[1])) <= 1.e-12)
	{
		//   Computation is impossible. Too smale values. Let's return some approximate value.

		//   buf_D  =  (uv[1] - bv[1])  *  (rb  - (uv[0] + bv[0]) /2.) * rhoInPrevTL[ indCurSqOx[0] ][ indCurSqOy[0] ];

		return fabs(uv[1] - bv[1]); //   fabs(uv[1] - bv[1]);
	}


	//   Integration. First step: under [ indCurSqOx[0]; indCurSqOx[1] ] square.

	//   A. Under triangle.

	if ((fabs(uv[1] - bv[1])) > 1.e-12)
	{
		//   Coefficients of slant line: x = a_SL *y  +  b_SL.

		a_SL = (uv[0] - bv[0]) / (uv[1] - bv[1]);

		b_SL = bv[0] - a_SL * bv[1];


		//   Integration under one cell triangle.
		if (fabs(a_SL) > 1.e-12)
		{
			buf_D = integUnderLeftTr_OneCell(
				par_a, //   -  Solution parameter.
				//
				lbDom, rbDom, //   -  Left and right boundaries of rectangular domain.
				//
				bbDom, ubDom, //   -  Botton and upper boundaries of rectangular domain.
				//
				bv[1], //   -  double Py,
				uv[1], //   -  double Qy,
				//
				a_SL,
				b_SL,
				mv[0], //   -  double Hx,
				//
				tau, iCurrTL, //   -  Index of current time layer.
				//
				indCurSqOx, //   -  Index of current square by Ox axis.
				indCurSqOy, //   -  Index of current square by Oy axis.
				//
				masOX, //   -  Massive of OX steps. Dimension = numOfOXSt +1.
				numOfOXSt, //   -  Number of OX steps.
				//
				masOY, //   -  Massive of OY steps. Dimension = numOfOYSt +1.
				numOfOYSt, //   -  Number of OY steps.
				//
				rhoInPrevTL_asV);


			/*
   double buf1_D_Ch, buf2_D_Ch;


   indCurSqOx[0] = 0;   indCurSqOx[1] = 1;

   indCurSqOy[0] = 0;   indCurSqOy[1] = 1;

   buf_D = integUnderRectAng_OneCell(
              par_a,                                  //   -  Solution parameter.
              //
              lbDom, rbDom,                           //   -  Left and right boundaries of rectangular domain.
              //
              bbDom, ubDom,                           //   -  Botton and upper boundaries of rectangular domain.
              //
              ( masOY[1] - masOY[0] ) /3.,            //   -  bv[1]
              2. * ( masOY[1] - masOY[0] ) /3.,       //   -  uv[1],
              //
              ( masOX[1] - masOX[0] ) /3.,            //   -  lv[0],                                  //   double Gx,
              2. * ( masOX[1] - masOX[0] ) /3.,       //   -  mv[0],                                  //   double Hx,
              //
              tau, iCurrTL,                            //   -  Index of current time layer.
              //
              indCurSqOx,                             //   -  Index of current square by Ox axis.
              indCurSqOy,                             //   -  Index of current square by Oy axis.
              //
              masOX,                                  //   -  Massive of OX steps. Dimension = numOfOXSt +1.
              numOfOXSt,                              //   -  Number of OX steps.
              //
              masOY,                                  //   -  Massive of OY steps. Dimension = numOfOYSt +1.
              numOfOYSt,                              //   -  Number of OY steps.
              //
              rhoInPrevTL );



   indCurSqOx[0] =-1;   indCurSqOx[1] = 0;

   indCurSqOy[0] = 0;   indCurSqOy[1] = 1;

   a_SL  =  1.;   b_SL  = -0.1;

   buf1_D_Ch = integUnderLeftTr_OneCell(
              par_a,                                  //   -  Solution parameter.
              //
              lbDom, rbDom,                           //   -  Left and right boundaries of rectangular domain.
              //
              bbDom, ubDom,                           //   -  Botton and upper boundaries of rectangular domain.
              //
              ( masOY[1] - masOY[0] ) /3.,            //   -  bv[1]
              2. * ( masOY[1] - masOY[0] ) /3.,       //   -  uv[1],
              //
              a_SL,
              b_SL,
              -( masOX[1] - masOX[0] ) /3.,           //   -  double Hx,
              //
              tau, iCurrTL,                            //   -  Index of current time layer.
              //
              indCurSqOx,                             //   -  Index of current square by Ox axis.
              indCurSqOy,                             //   -  Index of current square by Oy axis.
              //
              masOX,                                  //   -  Massive of OX steps. Dimension = numOfOXSt +1.
              numOfOXSt,                              //   -  Number of OX steps.
              //
              masOY,                                  //   -  Massive of OY steps. Dimension = numOfOYSt +1.
              numOfOYSt,                              //   -  Number of OY steps.
              //
              rhoInPrevTL );




   buf2_D_Ch = integUnderRightTr_OneCell(
              par_a,                                  //   -  Solution parameter.
              //
              lbDom, rbDom,                           //   -  Left and right boundaries of rectangular domain.
              //
              bbDom, ubDom,                           //   -  Botton and upper boundaries of rectangular domain.
              //
              ( masOY[1] - masOY[0] ) /3.,            //   -  bv[1]
              2. * ( masOY[1] - masOY[0] ) /3.,       //   -  uv[1],
              //
              a_SL,
              b_SL,
              -2. * ( masOX[1] - masOX[0] ) /3.,      //   -  double Gx,
              //
              tau, iCurrTL,                            //   -  Index of current time layer.
              //
              indCurSqOx,                             //   -  Index of current square by Ox axis.
              indCurSqOy,                             //   -  Index of current square by Oy axis.
              //
              masOX,                                  //   -  Massive of OX steps. Dimension = numOfOXSt +1.
              numOfOXSt,                              //   -  Number of OX steps.
              //
              masOY,                                  //   -  Massive of OY steps. Dimension = numOfOYSt +1.
              numOfOYSt,                              //   -  Number of OY steps.
              //
              rhoInPrevTL );
*/


			integ += buf_D;
		}
	}


	//   B. Under rectangle. Need to be cheking.

	if (wMvI == 1)
	{
		if (indCurSqOx[0] == indRB[0])
		{
			Hx = rb;
		}

		if (indCurSqOx[0] < indRB[0])
		{
			if (indCurSqOx[1] >= 0)
			{
				Hx = masOX[indCurSqOx[1]];
			}

			if (indCurSqOx[1] < 0)
			{
				Hx = h * indCurSqOx[1];
			}
		}


		buf_D = integUnderRectAng_OneCell(
			par_a, //   -  Solution parameter.
			//
			lbDom, rbDom, //   -  Left and right boundaries of rectangular domain.
			//
			bbDom, ubDom, //   -  Botton and upper boundaries of rectangular domain.
			//
			bv[1], //   -  double Py,
			uv[1], //   -  double Qy,
			//
			mv[0], //   -  double Gx,
			Hx, //   -  double Hx,
			//
			tau, iCurrTL, //   -  Index of current time layer.
			//
			indCurSqOx, //   -  Index of current square by Ox axis.
			indCurSqOy, //   -  Index of current square by Oy axis.
			//
			masOX, //   -  Massive of OX steps. Dimension = numOfOXSt +1.
			numOfOXSt, //   -  Number of OX steps.
			//
			masOY, //   -  Massive of OY steps. Dimension = numOfOYSt +1.
			numOfOYSt, //   -  Number of OY steps.
			//
			rhoInPrevTL_asV);

		integ += buf_D;
	}


	//   Second step: from "masOX[ indCurSqOx[1] ]" to "rb" by iteration.


	indCurSqOxToCh[0] = indCurSqOx[0] + 1;

	indCurSqOxToCh[1] = indCurSqOxToCh[0] + 1;


	for (j = indCurSqOx[0] + 1; j < indRB[0] + 1; j++)
	{
		//   If this is first cell we should integrate under triangle only.


		if (indCurSqOxToCh[1] > 0)
		{
			Gx = masOX[indCurSqOxToCh[0]];

			Hx = masOX[indCurSqOxToCh[1]];
		}


		if (indCurSqOxToCh[1] <= 0)
		{
			Gx = h * indCurSqOxToCh[0];

			Hx = h * indCurSqOxToCh[1];
		}


		if (j == indRB[0])
		{
			Hx = rb;
		}


		buf_D = integUnderRectAng_OneCell(
			par_a, //   -  Solution parameter.
			//
			lbDom, rbDom, //   -  Left and right boundaries of rectangular domain.
			//
			bbDom, ubDom, //   -  Botton and upper boundaries of rectangular domain.
			//
			bv[1], //   -  double Py,
			uv[1], //   -  double Qy,
			//
			Gx, //   -  double Gx,
			Hx, //   -  double Hx,
			//
			tau, iCurrTL, //   -  Index of current time layer.
			//
			indCurSqOxToCh, //   -  Index of current square by Ox axis.
			indCurSqOy, //   -  Index of current square by Oy axis.
			//
			masOX, //   -  Massive of OX steps. Dimension = numOfOXSt +1.
			numOfOXSt, //   -  Number of OX steps.
			//
			masOY, //   -  Massive of OY steps. Dimension = numOfOYSt +1.
			numOfOYSt, //   -  Number of OY steps.
			//
			rhoInPrevTL_asV);

		integ += buf_D;


		indCurSqOxToCh[0] += 1;

		indCurSqOxToCh[1] = indCurSqOxToCh[0] + 1;
	}


	return integ;
}


double integOfChan_SLRightSd(//   -  The domain is Channel with Slant Line on the right side.
	double par_a, //   -  Solution parameter.
	//
	double lbDom, //   -  Left and right boundaries of rectangular domain.
	double rbDom,
	//
	double bbDom, //   -  Botton and upper boundaries of rectangular domain.
	double ubDom,
	//
	double tau,
	int iCurrTL, //   -  Index of current time layer.
	//
	double* bv, int wTrPCI, //   -  Where travel point current (botton vertex) is.
	double* uv, int wTrPNI, //   -  Where travel point next (upper vertex) is.
	//
	int* indCurSqOx, //   -  Index by OX axis where bv and uv are.
	//
	double lb, int* indLB, //   -  Left boundary by Ox. Index by OX axis where lb is.
	//
	int* indCurSqOy, //   -  Index of current square by Oy axis.
	//
	double* masOX, //   -  Massive of OX steps. Dimension = numOfOXSt +1.
	int numOfOXSt, //   -  Number of OX steps.
	//
	double* masOY, //   -  Massive of OY steps. Dimension = numOfOYSt +1.
	int numOfOYSt, //   -  Number of OY steps.
	//
	double* rhoInPrevTL_asV)
{
	double mv[2], rv[2]; //   -  Middle and right vertices.

	int wMvI; //   -  Where middle vertex is.

	int indCurSqOxToCh[2]; //   -  Indices of current square by Ox axis to be changed. Under which we want to integrate.


	double h = masOX[1] - masOX[0];

	double a_SL, b_SL; //   -  Coefficients of slant line: x = a_SL *y  +  b_SL.


	double Gx, Hx; //   -  Left boundary for each integration.


	double integ = 0.;

	double buf_D;

	int j;


	//   Let's compute helpful values.

	if (uv[0] <= bv[0])
	{
		mv[0] = uv[0];
		mv[1] = uv[1];
		wMvI = wTrPNI;

		rv[0] = bv[0];
		rv[1] = bv[1];
	}

	if (uv[0] > bv[0])
	{
		mv[0] = bv[0];
		mv[1] = bv[1];
		wMvI = wTrPCI;

		rv[0] = uv[0];
		rv[1] = uv[1];
	}


	//   buf_D = fabs(mv[0] - lb) * fabs(uv[1] - bv[1])  +   fabs(uv[1] - bv[1]) * fabs(rv[0] - mv[0]) / 2.;

	//   return  buf_D;


	if ((fabs(uv[1] - bv[1])) <= 1.e-12)
	{
		//   Computation is impossible. Too smale values. Let's return some approximate value.

		//   buf_D  =  (uv[1] - bv[1])  *  ((uv[0] + bv[0]) /2.  -  lb) * rhoInPrevTL[ indCurSqOx[0] ][ indCurSqOy[0] ];

		return fabs(uv[1] - bv[1]); //   fabs(uv[1] - bv[1]);
	}


	//   First step: from "lb" to "masOX[ indCurSqOx[0] ]" by iteration.


	//   integ  += fabs( mv[0] - lb) * fabs(uv[1] - bv[1]);


	indCurSqOxToCh[0] = indLB[0];

	indCurSqOxToCh[1] = indCurSqOxToCh[0] + 1;

	for (j = indLB[0]; j < indCurSqOx[0]; j++)
	{
		//   If this is first cell we should integrate under rectangle only.


		if (indCurSqOxToCh[0] >= 0)
		{
			Gx = masOX[indCurSqOxToCh[0]];

			Hx = masOX[indCurSqOxToCh[1]];
		}


		if (indCurSqOxToCh[0] < 0)
		{
			Gx = h * indCurSqOxToCh[0];

			Hx = h * indCurSqOxToCh[1];
		}


		if (j == indLB[0])
		{
			Gx = lb;
		}


		buf_D = integUnderRectAng_OneCell(
			par_a, //   -  Solution parameter.
			//
			lbDom, rbDom, //   -  Left and right boundaries of rectangular domain.
			//
			bbDom, ubDom, //   -  Botton and upper boundaries of rectangular domain.
			//
			bv[1], //   -  double Py,
			uv[1], //   -  double Qy,
			//
			Gx, //   -  double Gx,
			Hx, //   -  double Hx,
			//
			tau, iCurrTL, //   -  Index of current time layer.
			//
			indCurSqOxToCh, //   -  Index of current square by Ox axis.
			indCurSqOy, //   -  Index of current square by Oy axis.
			//
			masOX, //   -  Massive of OX steps. Dimension = numOfOXSt +1.
			numOfOXSt, //   -  Number of OX steps.
			//
			masOY, //   -  Massive of OY steps. Dimension = numOfOYSt +1.
			numOfOYSt, //   -  Number of OY steps.
			//
			rhoInPrevTL_asV);

		integ += buf_D;


		indCurSqOxToCh[0] += 1;

		indCurSqOxToCh[1] = indCurSqOxToCh[0] + 1;
	}


	//   Integration. Second step: under [ indCurSqOx[0]; indCurSqOx[1] ] square.

	//   A. Under rectangle.

	if (wMvI == 1)
	{
		if (indCurSqOx[0] == indLB[0])
		{
			Gx = lb;
		}

		if (indCurSqOx[0] > indLB[0])
		{
			if (indCurSqOx[0] >= 0)
			{
				Gx = masOX[indCurSqOx[0]];
			}

			if (indCurSqOx[0] < 0)
			{
				Gx = h * indCurSqOx[0];
			}
		}


		buf_D = integUnderRectAng_OneCell(
			par_a, //   -  Solution parameter.
			//
			lbDom, rbDom, //   -  Left and right boundaries of rectangular domain.
			//
			bbDom, ubDom, //   -  Botton and upper boundaries of rectangular domain.
			//
			bv[1], //   -  double Py,
			uv[1], //   -  double Qy,
			//
			Gx, //   -  double Gx,
			mv[0], //   -  double Hx,
			//
			tau, iCurrTL, //   -  Index of current time layer.
			//
			indCurSqOx, //   -  Index of current square by Ox axis.
			indCurSqOy, //   -  Index of current square by Oy axis.
			//
			masOX, //   -  Massive of OX steps. Dimension = numOfOXSt +1.
			numOfOXSt, //   -  Number of OX steps.
			//
			masOY, //   -  Massive of OY steps. Dimension = numOfOYSt +1.
			numOfOYSt, //   -  Number of OY steps.
			//
			rhoInPrevTL_asV);

		integ += buf_D;
	}


	//   B. Under triangle.

	if ((fabs(uv[1] - bv[1])) > 1.e-12)
	{
		//   integ += fabs(uv[1] - bv[1]) * (rv[0] - mv[0]) /2.;

		//   Coefficients of slant line: x = a_SL *y  +  b_SL.

		a_SL = (uv[0] - bv[0]) / (uv[1] - bv[1]);

		b_SL = bv[0] - a_SL * bv[1];


		//   Integration under one cell triangle.

		if (fabs(a_SL) > 1.e-12)
		{
			buf_D = integUnderRightTr_OneCell(
				par_a, //   -  Solution parameter.
				//
				lbDom, rbDom, //   -  Left and right boundaries of rectangular domain.
				//
				bbDom, ubDom, //   -  Botton and upper boundaries of rectangular domain.
				//
				bv[1], //   -  double Py,
				uv[1], //   -  double Qy,
				//
				a_SL,
				b_SL,
				mv[0], //   -  double Gx,
				//
				tau, iCurrTL, //   -  Index of current time layer.
				//
				indCurSqOx, //   -  Index of current square by Ox axis.
				indCurSqOy, //   -  Index of current square by Oy axis.
				//
				masOX, //   -  Massive of OX steps. Dimension = numOfOXSt +1.
				numOfOXSt, //   -  Number of OX steps.
				//
				masOY, //   -  Massive of OY steps. Dimension = numOfOYSt +1.
				numOfOYSt, //   -  Number of OY steps.
				//
				rhoInPrevTL_asV);

			integ += buf_D;
		}
	}


	return integ;
}