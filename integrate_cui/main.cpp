#include <stdio.h>
#include <math.h>

double func1(double x)
{
	return exp(x) + 1;
}

double func2(double x)
{
	return pow(2.,x) + 1/log(2.);
}

double func3(double x)
{
	return x*exp(x);
}

int get_n_rect_trapez(double a, double b, double eps)
{
	return ((int) ((b - a) / sqrt(eps))) + 1;
}

int get_n_simpson(double a, double b, double eps)
{
	double h = 2.*pow(eps, 1.0 / 4.);
	double g = ((int) (b - a) / h);
	return g + 1;
}

double solve_rectangle(double a, double b, double eps, double (*f)(double))
{
	double n = get_n_rect_trapez(a, b, eps);
	double h = (b - a) / n;

	double x = 0;
	double prev = 0;
	for (int k = 1; k <= n; k++)
	{
		x = a + (k - 0.5)*h;
		prev += h * (*f)(x);
	}

	double curr = 0;
	while (fabs(prev - curr) / 3. > eps)
	{
		for (int k = 1; k <= n; k++)
		{
			x = a + (k - 0.5)*h;
			curr += h * (*f)(x);
		}
		n *= 2;
		prev = curr;
	};

	return prev;
}

double solve_trapezium(double a, double b, double eps, double(*f)(double))
{
	double n = get_n_rect_trapez(a, b, eps);
	double h = (b - a) / n;

	double x = 0;
	double prev = 0;
	for (int i = 1; i < n; i++)
	{
		x = a + i*h;
		prev += (*f)(x);
	}
	prev += (*f)(a) / 2;
	prev += (*f)(b) / 2;
	prev = prev*h;

	double curr = 0;
	while (fabs(prev - curr) / 3. > eps)
	{
		for (int i = 1; i < n; i++)
		{
			x = a + i*h;
			curr += (*f)(x);
		}
		curr += (*f)(a) / 2;
		curr += (*f)(b) / 2;
		curr *= h;
		n *= 2;
		prev = curr;
	};

	return prev;
}

double solve_simpson(double a, double b, double eps, double(*f)(double))
{
	double n = 2*get_n_simpson(a, b, eps);
	double h = (b - a) / n;

	double x = 0;
	double prev = 0;
	for (int i = 0; i <= n; i++)
	{
		x = a + i*h;
		prev += (*f)(x) * ((i == 0 || i == n) ? 1 : ((i & 1) == 0) ? 2 : 4);
	}
	
	prev *= h/3;

	double curr = 0;
	while (fabs(prev - curr) / 15. > eps)
	{
		for (int i = 0; i <= n; i++)
		{
			x = a + i*h;
			curr += (*f)(x) * ((i == 0 || i == n) ? 1 : ((i & 1) == 0) ? 2 : 4);
		}
		curr *= h/3;
		n *= 2;
		prev = curr;
	};
	return prev;
}

int main()
{
	const double eps = 1e-4;
	const double A = 0;
	const double B = 1;
	double r = solve_rectangle(A, B, eps, func1);
	double t = solve_trapezium(A, B, eps, func1);
	double s = solve_simpson(A, B, eps, func1);
	printf("Func1 = e^x + 1\n");
	printf("Integral by using the rectangle method = %f\n", r);
	printf("Integral by using the trapezium method = %f\n", t);
	printf("Integral by using the Simpson method = %f\n", s);
	r = solve_rectangle(A, B, eps, func2);
	t = solve_trapezium(A, B, eps, func2);
	s = solve_simpson(A, B, eps, func2);
	printf("\nFunc2 = 2^x + 1/ln2\n");
	printf("Integral by using the rectangle method = %f\n", r);
	printf("Integral by using the trapezium method = %f\n", t);
	printf("Integral by using the Simpson method = %f\n", s);

	r = solve_rectangle(A, B, eps, func3);
	t = solve_trapezium(A, B, eps, func3);
	s = solve_simpson(A, B, eps, func3);
	printf("\nFunc3 = x*e^x\n");
	printf("Integral by using the rectangle method = %f\n", r);
	printf("Integral by using the trapezium method = %f\n", t);
	printf("Integral by using the Simpson method = %f\n", s);
	return 0;
}