#include <iostream>
#include <fstream>
#include <gtest/gtest.h>
#include <gtest/gtest-spi.h>
#include "timer.h"
#include "common.h"

class TestBase : public testing::Test {
protected:
    double _accuracy;

    TestBase() {
        _accuracy = 1.0e-8;
        initCompOfGlVar();
    }

    virtual ~TestBase() {

        memClean();
    }

    virtual void SetUp() {

    }

    virtual void TearDown() {

    }

    int GetSize() {
        return C_numOfOXSt * C_numOfOYSt;
    }

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

    void print_matrix_to_file(int n, int m, double *data, std::string file_name, int precision = 8) {

        FILE * pFile = fopen(file_name.c_str(), "w");
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < m; ++j) {
                int k = i * n + j;
                switch (precision) {

                    case 8:
                        fprintf(pFile, "%20.14le ", data[k]);
                        break;
                }
            }
            fprintf(pFile, "\n ");
        }
        fclose(pFile);
    }
};

class cpu : public TestBase {
protected:

    cpu() {
    }

    virtual ~cpu() {
    }

    double *GetCpuToLevel(int level) {
        return solve_cpu_test(C_par_a, C_par_b, C_lbDom, C_rbDom, C_bbDom,
                C_ubDom, C_tau, C_numOfTSt, masOX, C_numOfOXSt, masOY,
                C_numOfOYSt, level, false);
    }

};

TEST_F(cpu, main_test) {
    const int finishLevel = 1;
    const int startLevel = 0;

    for (int level = startLevel; level < finishLevel; ++level) {
        std::cout << "level = " << level << std::endl;
        double *data = GetCpuToLevel(level);

        delete[] data;
    }
}