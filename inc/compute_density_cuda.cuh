#ifndef COMPUTE_DENSITY_CUDA_CUH
#define	COMPUTE_DENSITY_CUDA_CUH

#ifdef __CUDACC__
#define HOST __host__
#define DEVICE __device__
#else
#define HOST
#define DEVICE
#endif

extern double* compute_density_cuda_internal(double b, double lb, double rb, 						 double bb, double ub,
                        double tau, int time_step_count, int ox_length, int oy_length, double& norm, float& time);

extern double* compute_density_quad_cuda_internal(double b, double lb, double rb, double bb, double ub,
                        double tau, int time_step_count, int ox_length, int oy_length, double& norm, float& time);

struct c_ip_t
{
	int x;
	int y;


	DEVICE c_ip_t(int x, int y)
		: x(x), y(y)
	{
	}

	DEVICE c_ip_t()
		: x(0), y(0)
	{
	}


	DEVICE c_ip_t& operator=(const c_ip_t& p)
	{
		if (this != &p)
		{
			x = p.x;
			y = p.y;
		}

		return *this;
	}

	inline DEVICE  c_ip_t operator+(const c_ip_t& p) const
	{
		return c_ip_t(p.x + x, p.y + y);
	}

	inline DEVICE  c_ip_t& operator+=(int a)
	{
		x += a;
		y += a;
		return *this;
	}

	inline DEVICE c_ip_t& operator-=(int a)
	{
		x -= a;
		y -= a;
		return *this;
	}

	inline DEVICE c_ip_t& operator-(int a)
	{
		x -= a;
		y -= a;
		return *this;
	}

	inline DEVICE c_ip_t& operator+(int a)
	{
		x += a;
		y += a;
		return *this;
	}

	inline DEVICE bool operator==(const c_ip_t& p) const
	{
		return (x == p.x && y == p.y);
	}

	inline DEVICE int operator<(const c_ip_t& p)
	{
		return ((x < p.x) || ((x == p.x) && (y < p.y)));
	}

	inline DEVICE int operator>(const c_ip_t& p)
	{
		return ((x > p.x) || ((x == p.x) && (y > p.y)));
	}
};

struct c_dp4_t
{
	double x;
	double y;
	double x_initial;
	double y_initial;

	DEVICE c_dp4_t(double x, double y, double x_init, double y_init)
		: x(x), y(y), x_initial(x_init), y_initial(y_init)
	{
	}

	DEVICE c_dp4_t(double x, double y)
		: x(x), y(y), x_initial(0), y_initial(0)
	{
	}

	DEVICE c_dp4_t()
		: x(0), y(0), x_initial(0), y_initial(0)
	{
	}

	inline DEVICE c_dp4_t& operator=(const c_dp4_t& p)
	{
		if (this != &p)
		{
			x = p.x;
			y = p.y;
			x_initial = p.x_initial;
			y_initial = p.y_initial;
		}

		return *this;
	}
};

struct c_dp_t
{
	double x;
	double y;

	DEVICE c_dp_t(double x, double y)
		: x(x), y(y)
	{
	}

	DEVICE c_dp_t()
		: x(0), y(0)
	{
	}

	inline DEVICE c_dp_t& operator=(const c_dp_t& p)
	{
		if (this != &p)
		{
			x = p.x;
			y = p.y;
		}

		return *this;
	}

	inline DEVICE c_dp_t& operator=(const c_dp4_t& p)
	{
		x = p.x;
		y = p.y;


		return *this;
	}

	inline DEVICE bool operator==(const c_dp_t& p) const
	{
		return (x == p.x && y == p.y);
	}

	inline DEVICE int operator<(const c_dp_t& p)
	{
		return ((x < p.x) || ((x == p.x) && (y < p.y)));
	}

	inline DEVICE int operator>(const c_dp_t& p)
	{
		return ((x > p.x) || ((x == p.x) && (y > p.y)));
	}
};

#endif /* COMPUTE_DENSITY_CUDA_CUH */