#ifndef UTILS_H
#define	UTILS_H


inline void print_params(double b,
        double lb,
        double rb,
        double bb,
        double ub,
        double tau,
        int tl_count,
        int ox_length,
        int oy_length) {
    printf("b = %f\n", b);
    printf("lbDom = %f\n", lb);
    printf("rbDom = %f\n", rb);
    printf("bbDom = %f\n", bb);
    printf("ubDom = %f\n", ub);
    printf("tau = %f\n", tau);
    printf("Time level count = %d\n", tl_count);
    printf("ox length = %d\n", ox_length + 1);
    printf("oy length = %d\n", oy_length + 1);
}

inline void print_params(int index, int needed_index,
        double b,
        double lb,
        double rb,
        double bb,
        double ub,
        double tau,
        int tl,
        int tl_count,
        int ox_length,
        int oy_length,
        int cur_x,
        int cur_y,
        double value) {
    if (index == needed_index) {
        printf("index = %d\n", index);
        printf("b = %f\n", b);
        printf("lbDom = %f\n", lb);
        printf("rbDom = %f\n", rb);
        printf("bbDom = %f\n", bb);
        printf("ubDom = %f\n", ub);
        printf("tau = %f\n", tau);
        printf("Time level count = %d\n", tl_count);
        printf("current time level = %d\n", tl);
        printf("current x = %d\n", cur_x);
        printf("current y = %d\n", cur_y);
        printf("ox length = %d\n", ox_length + 1);
        printf("oy length = %d\n", oy_length + 1);
        printf("value = %f\n", value);
    }
}

inline void print_matrix_to_file(int n, int m, double *data, char* filename) {
    FILE *f = fopen(filename, "w");
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < m; ++j) {
            int k = i * n + j;
            fprintf(f, "%20.14le ", data[k]);
        }
        fprintf(f, "\n ");
    }
    fclose(f);
}

inline void print_matrix(double *a, int n, int m, int precision = 8) {
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
