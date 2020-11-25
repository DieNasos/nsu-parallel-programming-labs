#include <cstring>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <mpi.h>

void printOmega(double* omega, const int N, const double min_x, const double min_y, const double min_z) {
    // printing values in region's points
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            for (int k = 0; k < N; k++) {
                std::cout << min_x + i << ", " << min_y + j << ", " << min_z + k << ": " << omega[i*N*N + j*N + k] << std::endl;
            }
        }
    }
}

double phi(double x, double y, double z) {
    return x*x + y*y + z*z;
}

double ro(double x, double y, double z, const double a) {
    return 6 - a * phi(x,y,z);
}

double recalculate_layer(
        int base_z, int height, double* omega_part, double* tmp_omega_part,
        double delta_x, double delta_y, double delta_z,
        const int N, const double min_x, const double min_y, const double min_z, const double a) {
    // absolute layer number
    int abs_z = base_z + height;

    // if it's boundary layer
    if (abs_z == 0 || abs_z == N - 1) {
        // copying in new array without re-calculation
        memcpy(tmp_omega_part + height*N*N, omega_part + height*N*N, N * N * sizeof(double));
        return 0;
    }

    // else re-calculating every layer-element using the Jacobi iterative formula
    double max_delta = 0;
    double z = min_z + abs_z * delta_z;

    for (int i = 0; i < N; i++) {
        double x = min_x + i * delta_x;

        for (int j = 0; j < N; j++) {
            double y = min_y + j * delta_y;

            int cell = height*N*N + i*N + j;

            // boundary element => no re-calculation
            if (i == 0 || i == N-1 || j == 0 || j == N-1) {
                tmp_omega_part[cell] = omega_part[cell];
                continue;
            }

            // Jacobi formula:
            // {[phi^m(i+1, j, k) - phi^m(i-1, j, k)] / delta_x^2
            // + [phi^m(i, j+1, k) - phi^m(i, j-1, k)] / delta_y^2
            // + [phi^m(i, j, k+1) - phi^m(i, j, k-1)] / delta_z^2
            // - ro(x, y, z)}
            // / (2/delta_x^2 + 2/delta_y^2 + 2/delta_z^2)

            tmp_omega_part[cell] = ((omega_part[height*N*N + (i+1)*N + j] + omega_part[height*N*N + (i-1)*N + j]) / (delta_x * delta_x)
                                    + (omega_part[height*N*N + i*N + (j+1)] + omega_part[height*N*N + i*N + (j-1)]) / (delta_y * delta_y)
                                    + (omega_part[(height+1)*N*N + i*N + j]+omega_part[(height-1)*N*N + i*N + j]) / (delta_z * delta_z)
                                    - ro(x, y, z, a)) / (2/(delta_x*delta_x) + 2/(delta_y*delta_y) + 2/(delta_z*delta_z) + a);;
            max_delta = std::max(max_delta, std::abs(tmp_omega_part[cell] - omega_part[cell]));
        }
    }

    return max_delta;
}

double JacobiMethod(
        const double epsilon, const double a, const int N,
        const double min_x, const double min_y, const double min_z,
        const double max_x, const double max_y, const double max_z) {
    /*
     * finding the function phi satisfying:
     * d^2(phi)/d^2(x) + d^2(phi)/d^2(y) + d^2(phi)/d^2(z) - a*phi = ro
     */

    int number_of_processes, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &number_of_processes);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // if the matrix size isn't entirely divided by the number of processes, then we finish
    if (N % number_of_processes && rank == 0) {
        std::cout << "Invalid number of processes" << std::endl;
        return -1;
    }

    double start_time = MPI_Wtime();

    // distances between adjacent points
    double delta_x = (max_x - min_x) / (N - 1);
    double delta_y = (max_y - min_y) / (N - 1);
    double delta_z = (max_z - min_z) / (N - 1);

    int part_height = N / number_of_processes;  // height of area parts
    int part_base_z = rank * part_height - 1; // z-coordinate of part's base for current process

    // 3-dimensional array for storing function values at points of the part
    // +2 - adding an image of the upper and bottom outer boundary layers
    double* omega = new double[(part_height + 2) * N * N];
    double* tmp_omega = new double[(part_height + 2) * N * N];

    int iterationsCounter = 0;

    // STEP 1
    // setting the values of the desired function on the boundary of the area omega:
    // phi(x, y, z) = F(x, y, z) for i = 0, i = Nx, j = 0, j = Ny, k = 0 and k = Nz
    // + setting the initial approximation in the inner part of the area omega:
    // phi^0 = 0 for i = 1, ..., Nx-1, j = 0, ..., Ny-1, k = 0, ..., Nz-1.

    for (int i = 0; i < part_height + 2; i++) {
        int omega_z = i + part_base_z;    // Oz-coordinate in area omega
        double real_z = min_z + delta_z * omega_z;    // real Oz-coordinate

        for (int j = 0; j < N; j++) {
            double x = min_x + delta_x * j; // Ox-coordinate

            for (int k = 0; k < N; k++) {
                double y = min_y + delta_y * k;   // Oy-coordinate

                if (omega_z == 0 || omega_z == N-1 || j == 0 || j == N-1 || k == 0 || k == N-1) {
                    omega[i*N*N + j*N + k] = phi(x, y, real_z);    // calculating phi value on the border
                } else {
                    omega[i*N*N + j*N + k] = 0;   // setting the initial approximation in the inner part
                }
            }
        }
    }

    double max_delta_shared;

    // STEP 2: in cycle calculating the next approximation of the desired function until the finish criteria
    do {
        double max_delta = 0;

        // re-calculating boundary top layer
        double tmp_delta = recalculate_layer(part_base_z, 1, omega, tmp_omega, delta_x, delta_y, delta_z, N, min_x, min_y, min_z, a);
        max_delta = std::max(max_delta, tmp_delta);

        // re-calculating boundary bottom layer
        tmp_delta = recalculate_layer(part_base_z, part_height, omega, tmp_omega, delta_x, delta_y, delta_z, N, min_x, min_y, min_z, a);
        max_delta = std::max(max_delta, tmp_delta);

        // array for info about requests
        MPI_Request requests[4];

        if (rank != 0) {
            // sending top boundary layer to the previous process
            MPI_Isend(tmp_omega + N*N, N*N, MPI_DOUBLE, rank-1, 0, MPI_COMM_WORLD, &requests[0]);
            // receiving top boundary outer layer from the previous process
            MPI_Irecv(tmp_omega, N*N, MPI_DOUBLE, rank-1, 0, MPI_COMM_WORLD, &requests[2]);
            // waiting for the finish of operations
            MPI_Wait(&requests[0], MPI_STATUS_IGNORE);
            MPI_Wait(&requests[2], MPI_STATUS_IGNORE);
        }

        if (rank != number_of_processes - 1) {
            // sending bottom boundary layer to the next process
            MPI_Isend(tmp_omega + part_height*N*N, N*N, MPI_DOUBLE, rank+1, 0, MPI_COMM_WORLD, &requests[1]);
            // receiving bottom outer boundary layer from the next process
            MPI_Irecv(tmp_omega + (part_height+1)*N*N, N*N, MPI_DOUBLE, rank+1, 0, MPI_COMM_WORLD, &requests[3]);
            // waiting for the finish of operations
            MPI_Wait(&requests[1], MPI_STATUS_IGNORE);
            MPI_Wait(&requests[3], MPI_STATUS_IGNORE);
        }

        // re-calculating all elements of non-boundary layers
        for (int i = 2; i < part_height; i++) {
            double tmp_delta = recalculate_layer(part_base_z, i, omega, tmp_omega, delta_x, delta_y, delta_z, N, min_x, min_y, min_z, a);
            max_delta = std::max(max_delta, tmp_delta);
        }

        // “archiving” the fully recalculated area for this process
        memcpy(omega, tmp_omega, (part_height+2)*N*N * sizeof(double));

        // finding the maximum delta of all processes
        MPI_Reduce(&max_delta, &max_delta_shared, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
        // sending maximum delta to all other processes
        MPI_Bcast(&max_delta_shared, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

        iterationsCounter++;  // + ineration
    } while (max_delta_shared >= epsilon);

    delete[] tmp_omega;

    double* fullResult;
    if (rank == 0) {
        fullResult = new double[N*N*N];
    }

    // gathering counted area-parts
    MPI_Gather(omega + N*N, part_height*N*N, MPI_DOUBLE, fullResult, part_height*N*N, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    double end_time = MPI_Wtime();

    if (rank == 0) {
        // counting maximum deviation from the standard response
        double max_delta = 0;

        for (int layer = 0; layer < N; layer++){
            double z = min_z + layer * delta_z;

            for (int j = 0; j < N; j++) {
                double x = min_x + j * delta_x;

                for (int k = 0; k < N; k++) {
                    double y = min_y + k * delta_y;

                    max_delta = std::max(max_delta, std::abs(fullResult[layer*N*N + j*N + k] - phi(x, y, z)));
                }
            }
        }

        std::cout << "Answer: delta = " << max_delta << std::endl;
        printf("Time: %lf\n", end_time - start_time);
        std::cout << iterationsCounter << " cycle iterations" << std::endl;
        // printOmega(fullResult, N, min_x, min_y, min_z);

        delete[] fullResult;
    }

    delete[] omega;

    return end_time - start_time;
}

void JacobiMethodTest(
        const int repeats,
        const double epsilon, const double a, const int N,
        const double min_x, const double min_y, const double min_z,
        const double max_x, const double max_y, const double max_z) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    double current_time, best_time = std::numeric_limits<double>::max();

    for (int i = 1; i <= repeats; i++) {
        if (rank == 0) {
            std::cout << "Try " << i << "/" << repeats << std::endl;
        }

        current_time = JacobiMethod(epsilon, 1e5, N, min_x, min_y, min_z, max_x, max_y, max_z);

        if (rank == 0) {
            best_time = (current_time < best_time) ? current_time : best_time;
        }
    }

    if (rank == 0) {
        std::cout << "N = " << N << std::endl;
        std::cout << "epsilon = " << epsilon << std::endl;
        printf("Best time: %lf sec\n", best_time);
    }
}

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);

    /*
     * simulation area: [-1; 1] × [-1; 1] × [-1; 1]
     * desired function: phi(x, y, z) = x^2 + y^2 + z^2
     * right side of the equation: ro(x, y, z) = 6 - a * phi(x, y, z)
     * equation parameter: a = 10^5 (1e5)
     * convergence threshold: epsilon = 10^-8 (?)
     * initial approximation: phi^0 = 0
     */
    JacobiMethodTest(std::atoi(argv[1]), std::atof(argv[3]), 1e5, std::atoi(argv[2]), -1, -1, -1, 1, 1, 1);

    MPI_Finalize();

    return 0;
}