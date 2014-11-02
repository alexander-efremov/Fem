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

void print_matrix(double *a, int n, int m, int precision = 8) {
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

