#ifndef TEST_UTILS_H
#define	TEST_UTILS_H

struct dp_t1 {
	double x;
	double y;

	dp_t1(double x, double y)
		: x(x), y(y) {
	}

	dp_t1()
		: x(0), y(0) {
	}

	inline double operator[](int i) {
		return (i == 0) ? x : y;
	}

	inline dp_t1& operator=(const dp_t1& p) {
		if (this != &p) {
			x = p.x;
			y = p.y;
		}

		return *this;
	}

	inline dp_t1 operator+(const dp_t1& p) const {
		return dp_t1(p.x + x, p.y + y);
	}

	inline bool operator==(const dp_t1& p) const {
		return (x == p.x && y == p.y);
	}

	inline int operator<(const dp_t1& p) {
		return ((x < p.x) || ((x == p.x) && (y < p.y)));
	}

	inline int operator>(const dp_t1& p) {
		return ((x > p.x) || ((x == p.x) && (y > p.y)));
	}
};

inline void _print_matrix_to_file(int n, int m, double* data, char* filename)
{
	FILE* f = fopen(filename, "w");
	for (int i = 0; i < n; ++i)
	{
		for (int j = 0; j < m; ++j)
		{
			int k = i * n + j;
			fprintf(f, "%20.14le ", data[k]);
		}
		fprintf(f, "\n ");
	}
	fclose(f);
}

inline void _print_matrix(double* a, int n, int m, int precision = 8)
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



#endif /* TEST_UTILS_H */