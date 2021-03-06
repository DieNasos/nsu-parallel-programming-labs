#include <time.h>
#include "Geometry.h"

double SoLESolutionTest2(int size, double epsilon) {
    Matrix A(size, size);
    A.fill(2, 0);   // MATRIX_TYPE::SOLE_A_2

    Vector x(size);
    x.fill(1);  // VECTOR_TYPE::SOLE_X

    Vector b(size);
    b.fill(3);  // VECTOR_TYPE::SOLE_B_2

    Vector y(size);
    Vector multiplication(size);

    double tau;
    double b_length = b.length();

    struct timespec start_time, end_time;

    clock_gettime(CLOCK_MONOTONIC_RAW, &start_time);

    while (true) {
        // STEP 1/3: y^n = Ax^n - b
        y = A * x - b;

        // finish criteria: ||A * x - b|| / ||b|| < epsilon
        if (y.length() / b_length < epsilon) {
            break;
        }

        // STEP 2/3: tau^n = (y^n, Ay^n) / (Ay^n, Ay^n)
        multiplication = A * y;
        tau = y.scalarProduct(multiplication) / multiplication.scalarProduct(multiplication);

        // STEP 3/3: x^(n+1) = x^n - tau^n * y^n
        x = x - y * tau;
    }

    clock_gettime(CLOCK_MONOTONIC_RAW, &end_time);

    double time = end_time.tv_sec - start_time.tv_sec + 0.000000001 * (end_time.tv_nsec - start_time.tv_nsec);

    // std::cout << "Result: ";
    // x.println();
    printf("Time: %lf sec\n", time);

    return time;
}

double SoLESolutionTest1(int size, double epsilon) {
    Matrix A(size, size);
    A.fill(1, 0);   // MATRIX_TYPE::SOLE_A_1

    Vector x(size);
    x.fill(1);  // VECTOR_TYPE::SOLE_X

    Vector b(size);
    b.fill(2);  // VECTOR_TYPE::SOLE_B1

    Vector y(size);
    Vector multiplication(size);

    double tau;
    double b_length = b.length();

    struct timespec start_time, end_time;

    clock_gettime(CLOCK_MONOTONIC_RAW, &start_time);

    while (true) {
        // STEP 1/3: y^n = Ax^n - b
        y = A * x - b;

        // finish criteria: ||A * x - b|| / ||b|| < epsilon
        if (y.length() / b_length < epsilon) {
            break;
        }

        // STEP 2/3: tau^n = (y^n, Ay^n) / (Ay^n, Ay^n)
        multiplication = A * y;
        tau = y.scalarProduct(multiplication) / multiplication.scalarProduct(multiplication);

        // STEP 3/3: x^(n+1) = x^n - tau^n * y^n
        x = x - y * tau;
    }

    clock_gettime(CLOCK_MONOTONIC_RAW, &end_time);

    double time = end_time.tv_sec - start_time.tv_sec + 0.000000001 * (end_time.tv_nsec - start_time.tv_nsec);

    // std::cout << "Result: ";
    // x.println();
    printf("Time: %lf sec\n", time);

    return time;
}

void SoLESolutionPrintStatistic(int size, double epsilon, int repeats, int test_type) {
    double current_time;
    double best_time = size;

    for (int i = 1; i <= repeats; i++) {
        std::cout << "Try " << i << "/" << repeats << std::endl;

        switch (test_type) {
            case 1:
                current_time = SoLESolutionTest1(size, epsilon);
                break;
            case 2:
                current_time = SoLESolutionTest2(size, epsilon);
                break;
            default:
                break;
        }
        best_time = (current_time < best_time) ? current_time : best_time;
    }

    std::cout << "STATISTIC" << std::endl;
    std::cout << "> size = " << size << std::endl;
    std::cout << "> epsilon = " << epsilon << std::endl;
    std::cout << "> repeats = " << repeats << std::endl;
    printf("> Best time: %lf sec\n", best_time);
}

int main(int argc, char* argv[]) {
    SoLESolutionPrintStatistic(atoi(argv[1]), 0.0001, atoi(argv[2]), atoi(argv[3]));
    return 0;
}