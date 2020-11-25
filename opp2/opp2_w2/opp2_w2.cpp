#include <iostream>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <omp.h>

void matrixAndVectorParallelMultiplication(int size, double* matrix, double* vector, double* result, int first_process_line, int lines_per_process) {
    // parallel
    // result = matrix * vector
    for (int i = 0; i < lines_per_process; i++) {
        double current_value = 0;
        for (int j = 0; j < size; j++) {
            current_value += matrix[size * (first_process_line + i) + j] * vector[j];
        }
        // every process fills its own part of result vector
        result[first_process_line + i] = current_value;
    }
}

void parallelScalarProduct(int size, double* vector1, double* vector2, double* result_parts, int first_process_element, int elements_per_process) {
    // parallel
    // sum of result_parts elements = (vector1, vector2)
    double sum = 0;
    int rank = first_process_element / elements_per_process;
    for (int i = first_process_element; i < first_process_element + elements_per_process; i++) {
        sum += vector1[i] * vector2[i];
    }
    result_parts[rank] = sum;
}

void vectorAndVectorDifference(int size, double* vector1, double* vector2, double* result) {
    // returns result = vector1 - vector2
    for (int i = 0; i < size; i++) {
        result[i] = vector1[i] - vector2[i];
    }
}

void scalarAndVectorMultiplication(int size, double scalar, double* vector, double* result) {
    // returns result = scalar * vector
    for (int i = 0; i < size; i++) {
        result[i] = scalar * vector[i];
    }
}

double vectorLength(int size, double* vector) {
    // returns ||vector||
    double length = 0;
    for (int i = 0; i < size; i++) {
        length += vector[i] * vector[i];
    }
    return std::sqrt(length);
}

void fillMatrixAndVector(int test_type, int size, double* matrix, double* vector) {
    switch (test_type) {
        // filling matrix and vector according to test type
        case 1:
            // filling matrix A
            // result:
            // [2][1][1]...
            // [1][2][1]...
            // [1][1][2]...
            // ...
            for (int i = 0; i < size; i++) {
                for (int j = 0; j < size; j++) {
                    if (i == j) {
                        matrix[i * size + j] = 2;
                    }
                    else {
                        matrix[i * size + j] = 1;
                    }
                }
            }
            // filling vector b
            for (int i = 0; i < size; i++) {
                vector[i] = size + 1;
            }
            break;
        case 2:
            // filling matrix A
            // result:
            // [1][0][0]...
            // [1][1][0]...
            // [1][1][1]...
            // ...
            for (int i = 0; i < size; i++) {
                for (int j = 0; j < size; j++) {
                    if (j < i + 1) {
                        matrix[i * size + j] = 1;
                    }
                    else {
                        matrix[i * size + j] = 0;
                    }
                }
            }
            // filling vector b
            for (int i = 0; i < size; i++) {
                vector[i] = i + 1;
            }
            break;
        default:
            break;
    }
}

void printVector(int size, double* vector) {
    for (int i = 0; i < size; i++) {
        std::cout << vector[i] << " ";
    }
}

double SoLESolutionParallelTest(int num_threads, int size, double epsilon, int test_type) {
    // finding solution of A*x - b = 0

    double* A = (double*) malloc(size * size * sizeof(double));
    double* x = (double*) calloc(size, sizeof(double)); // filling with 0s
    double* b = (double*) malloc(size * sizeof(double));
    double* y = (double*) malloc(size * sizeof(double));
    double* mult = (double*) malloc(size * sizeof(double));

    // parts of scalar products
    double* sp1_parts = (double*) malloc(num_threads * sizeof(double));
    double* sp2_parts = (double*) malloc(num_threads * sizeof(double));

    fillMatrixAndVector(test_type, size, A, b);

    double tau;

    double b_length = vectorLength(size, b);
    double y_length;

    double start_time = omp_get_wtime();

#pragma omp parallel
    {
        // parallel section BEGINS
        int num_threads = omp_get_num_threads();    // number of processes
        int rank = omp_get_thread_num();    // id of current process
        int lines_per_process = size / num_threads;    // number of matrix lines per single process
        int first_process_line = rank * lines_per_process;  // first matrix line current process will work with

        while (true) {
            // STEP 1/3: y^n = Ax^n - b
            // mult = A * x
            matrixAndVectorParallelMultiplication(size, A, x, mult, first_process_line, lines_per_process);
#pragma omp barrier
            {
                // all the processes reached this point => continuing
            };
#pragma omp single
            {
                // this block works on the only one process

                // y = mult - b
                vectorAndVectorDifference(size, mult, b, y);

                // finish criteria: ||A * x - b|| / ||b|| < epsilon
                y_length = vectorLength(size, y);
            };

            if (y_length / b_length < epsilon) { break; }

            // STEP 2/3: tau^n = (y^n, Ay^n) / (Ay^n, Ay^n)
            // mult = A * y
            matrixAndVectorParallelMultiplication(size, A, y, mult, first_process_line, lines_per_process);
            // counting first scalar product's parts
            parallelScalarProduct(size, y, mult, sp1_parts, first_process_line, lines_per_process);
            // counting second scalar product's parts
            parallelScalarProduct(size, mult, mult, sp2_parts, first_process_line, lines_per_process);
#pragma omp barrier
            {
                // all the processes reached this point => continuing
            };
#pragma omp single
            {
                // this block works on the only one process

                double sp1 = 0;
                double sp2 = 0;
                // uniting parts of scalar products together
                for (int i = 0; i < num_threads; i++) {
                    sp1 += sp1_parts[i];
                    sp2 += sp2_parts[i];
                }
                tau = sp1 / sp2;    // counting tau

                // STEP 3/3: x^(n+1) = x^n - tau^n * y^n
                // mult = tau * y
                scalarAndVectorMultiplication(size, tau, y, mult);
                // x = x - mult
                vectorAndVectorDifference(size, x, mult, x);
            };
        };
        // parallel section ENDS
    }

    double time = omp_get_wtime() - start_time;

//    std::cout << "Result: ";
//    printVector(size, x);
//    std::cout << std::endl;

    printf("Time: %lf sec\n", time);

    return time;
}

void SoLESolutionPrintStatistic(int num_threads, int size, double epsilon, int repeats, int test_type) {
    double current_time, best_time = size;

    for (int i = 1; i <= repeats; i++) {
        std::cout << "Try " << i << "/" << repeats << std::endl;

        // parallel test
        current_time = SoLESolutionParallelTest(num_threads, size, epsilon, test_type);

        best_time = (current_time < best_time) ? current_time : best_time;
    }

    std::cout << "STATISTIC" << std::endl;
    std::cout << "> number of processes = " << num_threads << std::endl;
    std::cout << "> size = " << size << std::endl;
    std::cout << "> epsilon = " << epsilon << std::endl;
    std::cout << "> repeats = " << repeats << std::endl;
    printf("> BEST TIME: %lf sec\n", best_time);
}

int main(int argc, char* argv[]) {
    // ARGUMENTS:
    // argv[1] == number of processes
    // argv[2] == size
    // argv[3] == number of tests
    // argv[4] == test ID

    omp_set_num_threads(atoi(argv[1])); // setting the number of processes

    SoLESolutionPrintStatistic(atoi(argv[1]), atoi(argv[2]), 0.0001, atoi(argv[3]), atoi(argv[4]));

    return 0;
}