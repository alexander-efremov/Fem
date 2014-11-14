#include <stdio.h>
#include <math.h>

double func(double x)
{
	return exp(x) + 1;
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

double solve_rectangle(double a, double b, double eps)
{
	double n = get_n_rect_trapez(a, b, eps);
	double h = (b - a) / n;

	double x = 0;
	double prev = 0;
	for (int k = 1; k <= n; k++)
	{
		x = a + (k - 0.5)*h;
		prev += h * func(x);
	}

	double curr = 0;
	while (fabs(prev - curr) / 3. > eps)
	{
		for (int k = 1; k <= n; k++)
		{
			x = a + (k - 0.5)*h;
			curr += h * func(x);
		}
		n *= 2;
		prev = curr;
	};

	return prev;
}

double solve_trapezium(double a, double b, double eps)
{
	double n = get_n_rect_trapez(a, b, eps);
	double h = (b - a) / n;

	double x = 0;
	double prev = 0;
	for (int i = 1; i < n; i++)
	{
		x = a + i*h;
		prev += func(x);
	}
	prev += func(a) / 2;
	prev += func(b) / 2;
	prev = prev*h;

	double curr = 0;
	while (fabs(prev - curr) / 3. > eps)
	{
		for (int i = 1; i < n; i++)
		{
			x = a + i*h;
			curr += func(x);
		}
		curr += func(a) / 2;
		curr += func(b) / 2;
		curr *= h;
		n *= 2;
		prev = curr;
	};

	return prev;
}

double solve_simpson(double a, double b, double eps)
{
	double n = 2*get_n_simpson(a, b, eps);
	double h = (b - a) / n;

	double x = 0;
	double prev = 0;
	for (int i = 0; i <= n; i++)
	{
		x = a + i*h;
		prev += func(x) * ((i == 0 || i == n) ? 1 : ((i & 1) == 0) ? 2 : 4);
	}
	
	prev *= h/3;

	double curr = 0;
	while (fabs(prev - curr) / 15. > eps)
	{
		for (int i = 0; i <= n; i++)
		{
			x = a + i*h;
			curr += func(x) * ((i == 0 || i == n) ? 1 : ((i & 1) == 0) ? 2 : 4);
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
	double r = solve_rectangle(A, B, eps);
	double t = solve_trapezium(A, B, eps);
	double s = solve_simpson(A, B, eps);
	printf("Rectangle = %f\n", r);
	printf("Trapezium = %f\n", t);
	printf("Simpson = %f\n", s);
	return 0;
}