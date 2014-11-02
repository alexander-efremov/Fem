#ifndef UTILS_H
#define	UTILS_H

void print_matrix_to_file(int n, int m, double *data, std::string filename)
{
    FILE *f = fopen(filename.c_str(), "w");
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

#endif	/* UTILS_H */

