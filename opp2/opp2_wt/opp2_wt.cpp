#include <iostream>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <time.h>

void matrixAndVectorMultiplication(int size, double* matrix, double* vector, double* result) {
    // non-parallel
    // result = matrix * vector
    for (int i = 0; i < size; i++) {
        double current_value = 0;
        for (int j = 0; j < size; j++) {
            current_value += matrix[size * i + j] * vector[j];
        }
        result[i] = current_value;
    }
}

double scalarProduct(int size, double* vector1, double* vector2) {
    // non-parallel
    // returns (vector1, vector2)
    double sum = 0;
    for (int i = 0; i < size; i++) {
        sum += vector1[i] * vector2[i];
    }
    return sum;
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

double SoLESolutionNonParallelTest(int size, double epsilon, int test_type) {
    // finding solution of A*x - b = 0

    double* A = (double*) malloc(size * size * sizeof(double));
    double* x = (double*) calloc(sizeof(double), size); // filling with 0s
    double* b = (double*) malloc(size * sizeof(double));
    double* y = (double*) malloc(size * sizeof(double));

    // for multiplications
    double* mult = (double*) malloc(size * sizeof(double));

    fillMatrixAndVector(test_type, size, A, b);

    double tau;
    double b_length = vectorLength(size, b);

    struct timespec start_time, end_time;

    clock_gettime(CLOCK_MONOTONIC_RAW, &start_time);

    while (true) {
        // STEP 1/3: y^n = Ax^n - b
        // mult = A * x
        matrixAndVectorMultiplication(size, A, x, mult);
        // y = mult - b
        vectorAndVectorDifference(size, mult, b, y);

        // finish criteria: ||A * x - b|| / ||b|| < epsilon
        double y_length = vectorLength(size, y);
        if (y_length / b_length < epsilon) { break; }

        // STEP 2/3: tau^n = (y^n, Ay^n) / (Ay^n, Ay^n)
        // mult = A * y
        matrixAndVectorMultiplication(size, A, y, mult);
        // counting tau
        tau = scalarProduct(size, y, mult) / scalarProduct(size, mult, mult);

        // STEP 3/3: x^(n+1) = x^n - tau^n * y^n
        // mult = tau * y
        scalarAndVectorMultiplication(size, tau, y, mult);
        // x = x - mult
        vectorAndVectorDifference(size, x, mult, x);
    }

    clock_gettime(CLOCK_MONOTONIC_RAW, &end_time);

    double time = end_time.tv_sec - start_time.tv_sec + 0.000000001 * (end_time.tv_nsec - start_time.tv_nsec);

//    std::cout << "Result: ";
//    printVector(size, x);
//    std::cout << std::endl;

    printf("Time: %lf sec\n", time);

    return time;
}

void SoLESolutionPrintStatistic(int size, double epsilon, int repeats, int test_type) {
    double current_time, best_time = size;

    for (int i = 1; i <= repeats; i++) {
        std::cout << "Try " << i << "/" << repeats << std::endl;

        // non-parallel test
        current_time = SoLESolutionNonParallelTest(size, epsilon, test_type);

        best_time = (current_time < best_time) ? current_time : best_time;
    }

    std::cout << "STATISTIC" << std::endl;
    std::cout << "> size = " << size << std::endl;
    std::cout << "> epsilon = " << epsilon << std::endl;
    std::cout << "> repeats = " << repeats << std::endl;
    printf("> BEST TIME: %lf sec\n", best_time);
}

int main(int argc, char* argv[]) {
    // ARGUMENTS:
    // argv[1] == size
    // argv[2] == number of repeats
    // argv[3] == test ID

    SoLESolutionPrintStatistic(atoi(argv[1]), 0.0001, atoi(argv[2]), atoi(argv[3]));

    return 0;
}