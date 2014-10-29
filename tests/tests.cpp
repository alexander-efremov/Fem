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
        // You can do clean-up work that doesn't throw exceptions here.
        memClean();
    }

    // If the constructor and destructor are not enough for setting up
    // and cleaning up each test, you can define the following methods:

    virtual void SetUp() {
        // Code here will be called immediately after the constructor (right
        // before each test).
    }

    virtual void TearDown() {
        // Code here will be called immediately after each test (right
        // before the destructor).
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

class cputest : public TestBase {
protected:

    cputest() {
    }

    virtual ~cputest() {
    }

    double *GetCpuToLevel(int level) {
        return solve_cpu_test(C_par_a, C_par_b, C_lbDom, C_rbDom, C_bbDom,
                C_ubDom, C_tau, C_numOfTSt, masOX, C_numOfOXSt, masOY,
                C_numOfOYSt, level, false);
    }

};

//TEST_F(gputest, main_test)
//{

//	const int finishLevel = 7;
//	const int startLevel = 0;
//	const bool isComputeDiff = true;
//	const bool isOneTl = false;
//
//	for (int level = startLevel; level < finishLevel; ++level)
//	{
//                std::cout << "level = " << level << std::endl;
//		double *data = NULL;
//
//		ComputeParameters *p = new ComputeParameters(level, true, isComputeDiff);
//		ASSERT_TRUE(p->result != NULL);
//		std::cout << *p << std::endl;
//
//		float gpu_time = solve_at_gpu(p, isOneTl, isComputeDiff);
//		ASSERT_TRUE(gpu_time != -1);
//		if (isOneTl)
//		{		data = GetCpuToLevel1TL(level, isComputeDiff); }
//		else
//		{  data = GetCpuToLevel(level, isComputeDiff); }
//
//		// data = _modelDataProvider.GetModelData(level);
//
//		printf("%s\n", "Checking is started");
//		for (int i = 0; i < p->get_real_matrix_size(); i++)
//		{
//			ASSERT_NEAR(data[i], p->result[i], __FLT_EPSILON__) << "i = " <<  i << std::endl;
//		}
//                printf("%s\n\n", "Checking is done");
//
//		std::ostringstream oss; 
//
//		char *s = "diff_gpu_"; 
//		oss << s << p->t_count << ".bin"; 
//
//		std::string name(oss.str());
//
//		if (isComputeDiff)
//		{
//			print_matrix_to_file(p->get_real_x_size(), p->get_real_y_size(), p->diff, name); 
//		}
//
//		std::ostringstream oss1; 
//		char *s1 = "gpu_result_"; 
//		oss1 << s1 << p->t_count << ".bin"; 
//		std::string name1(oss1.str());
//
//		print_matrix_to_file(p->get_real_x_size(), p->get_real_y_size(), p->result, name1); 
//	
//
//		std::ostringstream oss2; 
//		char *s2 = "cpu_result_"; 
//		oss2 << s2 << p->t_count << ".bin"; 
//		std::string name2(oss2.str());
//
//		print_matrix_to_file(p->get_real_x_size(), p->get_real_y_size(), data, name2); 
//
//		delete p;
//		delete[] data;
//	}
//}