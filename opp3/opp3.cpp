#include "Matrix.h"

double MatrixAndMatrixMultiplicationTest(int size1, int size2, int size3, int gridSizeX, int gridSizeY) {
    int rank;
    // every process gets its rank
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // creating matrices
    Matrix A(size1, size2, gridSizeX, gridSizeY);
    Matrix B(size2, size3, gridSizeX, gridSizeY);
    Matrix C(size1, size3, gridSizeX, gridSizeY);

    if (rank == 0) {
        // main process fills matrices A & B
        A.fill(1);  // MATRIX_TYPE::A
        B.fill(2);  // MATRIX_TYPE::B
    }

    // every process fills its copy of C with 0s
    C.fill(0);  // MATRIX_TYPE::ZERO

    struct timespec start_time, end_time;

    clock_gettime(CLOCK_MONOTONIC_RAW, &start_time);

    // counting multiplication of matrices
    C = A * B;

    clock_gettime(CLOCK_MONOTONIC_RAW, &end_time);

    double time = end_time.tv_sec - start_time.tv_sec + 0.000000001 * (end_time.tv_nsec - start_time.tv_nsec);

    if (rank == 0) {
        // printing result matrix and time elapsed
        //C.print();
        printf("Time: %lf sec\n", time);
    }

    return time;
}

void MultiplicationPrintStatistic(int size1, int size2, int size3, int gridSizeX, int gridSizeY, int repeats) {
    int current_process;
    MPI_Comm_rank(MPI_COMM_WORLD, &current_process);

    double current_time;
    double best_time = std::numeric_limits<double>::max();

    for (int i = 1; i <= repeats; i++) {
        if (current_process == 0) {
            std::cout << "Try " << i << "/" << repeats << std::endl;
        }

        current_time = MatrixAndMatrixMultiplicationTest(size1, size2, size3, gridSizeX, gridSizeY);

        if (current_process == 0) {
            best_time = (current_time < best_time) ? current_time : best_time;
        }
    }

    if (current_process == 0) {
        std::cout << "STATISTIC" << std::endl;
        std::cout << "> size1 = " << size1 << std::endl;
        std::cout << "> size2 = " << size1 << std::endl;
        std::cout << "> size3 = " << size1 << std::endl;
        std::cout << "> gridSizeX = " << gridSizeX << std::endl;
        std::cout << "> gridSizeY = " << gridSizeY << std::endl;
        std::cout << "> repeats = " << repeats << std::endl;
        printf("> Best time: %lf sec\n", best_time);
    }
}

int main(int argc, char **argv) {
    // argv[1], argv[2], argv[3] == matrices sizes
    // argv[4], argv[5] == 2D-grid sizes
    // argv[6] == number of test repeats

    MPI_Init(&argc, &argv); // initializing MPI

    MultiplicationPrintStatistic(atoi(argv[1]), atoi(argv[2]), atoi(argv[3]), atoi(argv[4]), atoi(argv[5]), atoi(argv[6]));

    MPI_Finalize();

    return 0;
}