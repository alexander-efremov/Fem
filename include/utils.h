#ifndef UTILS_H
#define	UTILS_H

void print_matrix_to_file(int n, int m, double *data, std::string filename) {
    FILE *f = fopen(filename.c_str(), "w");
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < m; ++j) {
            int k = i * n + j;
            fprintf(f, "%20.14le ", data[k]);
        }
        fprintf(f, "\n ");
    }
    fclose(f);
}

/*void compute_diff_write_to_file(double *result, int tl, int n, int m, double tau)
{
    double c_h = 1. / n;
    double *diff = new double[(n + 1) * (m + 1)]();

    std::ostringstream oss;
    oss << "diff_cpu_" << tl << ".bin";

    std::string name(oss.str());

    for (int j = 0; j <= m; j++)
    {
        for (int i = 0; i <= n; i++)
        {
            int opt = (n + 1) * j + i;
            double f = analytical_solution(tl*tau, i*c_h, j * c_h);
            diff [opt] = fabs(result[opt] - f);
        }
    }
    print_matrix_to_file(n + 1, m + 1, diff, name);
    delete[] diff;
}*/

void print_matrix(int n, int m, double *a, int precision = 8) {
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < m; ++j) {
            int k = i * n + j;
            switch (precision) {
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

#endif	/* UTILS_H */

