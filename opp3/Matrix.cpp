#include "Matrix.h"

Matrix::Matrix(int rows, int columns) {
    this->rows = rows;
    this->columns = columns;
    // matrix itself is empty for now
}

Matrix::Matrix(int rows, int columns, int gridSizeX, int gridSizeY) {
    this->rows = rows;
    this->columns = columns;
    this->gridSizeX = gridSizeX;
    this->gridSizeY = gridSizeY;
    // matrix itself is empty for now
}

Matrix::Matrix(int rows, int columns, double *source) {
    this->rows = rows;
    this->columns = columns;
    // filling matrix with elements of source-array
    if (source != NULL) {
        for (int i = 0; i < this->rows; i++) {
            std::vector<double> line;
            this->matrix.push_back(line);

            for (int j = 0; j < this->columns; j++) {
                this->matrix[i].push_back(source[i * rows + j]);
            }
        }
    } else {
        //std::cout << "WARNING :: MATRIX :: CONSTRUCTOR :: SOURCE ARRAY IS EMPTY" << std::endl;
        this->fill(0);  // MATRIX_TYPE::ZERO
    }
}

Matrix& Matrix::operator=(Matrix other_matrix) {
    // step 0: try of self equating => no need to continue
    if (this == &other_matrix) {
        return *this;
    }

    // step 1: current matrix is bigger than other => removing cells

    // number of current matrix's rows > other's
    while (this->matrix.size() > other_matrix.matrix.size()) {
        // removing elements from last row
        while (!this->matrix[this->matrix.size() - 1].empty()) {
            this->matrix[this->matrix.size() - 1].pop_back();
        }

        // removing last row
        this->matrix.pop_back();
    }

    // checking every row
    for (int i = 0; i < this->rows; i++) {
        // current matrix's row is bigger than other's
        if (this->matrix[i].size() > other_matrix.matrix[i].size()) {
            // removing extra elements
            while (this->matrix[i].size() != other_matrix.matrix[i].size()) {
                this->matrix[i].pop_back();
            }
        }
    }

    // step 2: other matrix is bigger than current => adding cells (with 0s for now)
    for (int i = this->rows - 1; i < other_matrix.rows; i++) {

        // adding row
        if (i > this->rows - 1) {
            std::vector<double> line;
            this->matrix.push_back(line);
        }

        // adding elements to row
        for (int j = this->matrix[i].size(); j < other_matrix.matrix[i].size(); j++) {
            this->matrix[i].push_back(0);
        }
    }

    // step 3: copying values of sizes
    this->rows = other_matrix.rows;
    this->columns = other_matrix.columns;

    // step 4: filling current matrix with values from other
    for (int i = 0; i < this->rows; i++) {
        for (int j = 0; j < this->matrix[i].size(); j++) {
            this->matrix[i][j] = other_matrix.matrix[i][j];
        }
    }

    return *this;
}

Matrix Matrix::operator*(Matrix second_matrix) {
    // creating new communicator for grid
    int p[2];
    p[0] = this->gridSizeX;
    p[1] = this->gridSizeY;
    MPI_Comm grid;
    MPI_Cart_create(MPI_COMM_WORLD, 2, p, new int[2]{0}, 0, &grid);

    // getting rank & coordinates of current process
    int grid_rank, grid_coords[2];
    MPI_Comm_rank(grid, &grid_rank);
    MPI_Cart_coords(grid, grid_rank, 2, grid_coords);

    // creating communicators for sub-grids
    MPI_Comm sub_grids[2];
    MPI_Cart_sub(grid, new int[2]{1, 0}, &sub_grids[1]);
    MPI_Cart_sub(grid, new int[2]{0, 1}, &sub_grids[0]);

    // getting sizes of stripes in mat1 & mat2 + sub-matrices of result matrix
    int sub_sizes[2];
    sub_sizes[0] = this->rows / p[0];
    sub_sizes[1] = second_matrix.columns / p[1];

    // arrays for content of matrices
    double *mat1, *mat2, *mult = NULL;

    double *sub_mat1 = new double[sub_sizes[0] * this->columns];
    double *sub_mat2 = new double[this->columns * sub_sizes[1]];
    double *sub_mult = new double[sub_sizes[0] * sub_sizes[1]];
    // sub_mat1 == local sub-matrix of mat1, including gridSizeX horizontal stripes
    // sub_mat2 == local sub-matrix of mat2, including gridSizeY vertical stripes
    // sub_mult == local sub-matrix of multiplication matrix

    // counts & displacements
    int *count_2, *displs_2, *count_mult, *displs_mult;

    // array of data-types
    MPI_Datatype types[3] {MPI_DOUBLE, MPI_UB};

    if (grid_rank == 0) {
        // displacements
        MPI_Aint displs[2]{0, (MPI_Aint) (sizeof(double) * sub_sizes[1])};

        // init
        displs_2 = new int[p[1]];
        count_2 = new int[p[1]];
        displs_mult = new int[p[0]*p[1]];
        count_mult = new int[p[0]*p[1]];
        for (int i = 0; i < p[0]; i++) {
            for (int j = 0; j < p[1]; j++) {
                if (i == 0) {
                    // executes once
                    displs_2[j] = j;
                    count_2[j] = 1;
                }
                displs_mult[i * p[1] + j] = i * p[1] * sub_sizes[0] + j;
                count_mult[i * p[1] + j] = 1;
            }
        }

        // creating new data-type (types[2])
        MPI_Type_vector(this->columns, sub_sizes[1], second_matrix.columns, MPI_DOUBLE, &types[0]);
        MPI_Type_struct(2, new int[2]{1, 1}, displs, types, &types[2]);
        MPI_Type_commit(&types[2]);

        // getting content of matrices & allocating memory for result
        mat1 = this->getAsArray();
        mat2 = second_matrix.getAsArray();
        mult = new double[this->rows * second_matrix.columns];
    }

    if (!grid_coords[0]) {
        // main process sends horizontal mat2-stripes in y-coordinate
        MPI_Scatterv(mat2, count_2, displs_2, types[2], sub_mat2, this->columns * sub_sizes[1], MPI_DOUBLE, 0, sub_grids[0]);
    }

    if (!grid_coords[1]) {
        // main process sends horizontal mat1-stripes in x-coordinate
        MPI_Scatter(mat1, sub_sizes[0] * this->columns, MPI_DOUBLE, sub_mat1, sub_sizes[0] * this->columns, MPI_DOUBLE, 0, sub_grids[1]);
    }

    // sending sub-matrices of mat1 in y dimension
    MPI_Bcast(sub_mat1, sub_sizes[0] *  this->columns, MPI_DOUBLE, 0, sub_grids[0]);
    // sending sub-matrices of mat2 in x dimension
    MPI_Bcast(sub_mat2, this->columns * sub_sizes[1], MPI_DOUBLE, 0, sub_grids[1]);

    // counting sub-matrices of result matrix in every process
    for (int i = 0; i < sub_sizes[0]; i++) {
        for (int j = 0; j < sub_sizes[1]; j++) {
            sub_mult[sub_sizes[1] * i + j] = 0;
            for (int k = 0; k < this->columns; k++) {
                sub_mult[sub_sizes[1] * i + j] += sub_mat1[this->columns * i + k] * sub_mat2[sub_sizes[1] * k + j];
            }
        }
    }

    // gathering sub-matrices of multiplication in main process
    MPI_Gatherv(sub_mult, sub_sizes[0]*sub_sizes[1], MPI_DOUBLE, mult, count_mult, displs_mult, types[2], 0, grid);

    delete[] sub_mat1;
    delete[] sub_mat2;
    delete[] sub_mult;
    MPI_Comm_free(&grid);
    MPI_Comm_free(&sub_grids[0]);
    MPI_Comm_free(&sub_grids[1]);

    // creating multiplication matrix object
    // will be EMPTY on every process except of MAIN (0)!
    Matrix multMat(this->rows, second_matrix.columns, mult);

    if (grid_rank == 0) {
        delete[] mat1;
        delete[] mat2;
        delete[] mult;
        delete[] count_2;
        delete[] count_mult;
        delete[] displs_2;
        delete[] displs_mult;
        MPI_Type_free(&types[0]);
        MPI_Type_free(&types[2]);
    }

    return multMat;
}

void Matrix::setContent(double *source) {
    // filling matrix with elements of source-array
    if (source != NULL) {
        for (int i = 0; i < this->rows; i++) {
            std::vector<double> line;
            this->matrix.push_back(line);

            for (int j = 0; j < this->columns; j++) {
                this->matrix[i].push_back(source[i * rows + j]);
            }
        }
    } else {
        //std::cout << "WARNING :: MATRIX :: CONSTRUCTOR :: SOURCE ARRAY IS EMPTY" << std::endl;
        this->fill(0);  // MATRIX_TYPE::ZERO
    }
}

void Matrix::setGridSizes(int x, int y) {
    this->gridSizeX = x;
    this->gridSizeY = y;
}

double* Matrix::getAsArray() {
    // returns matrix in array form

    double* array = new double[this->rows * this->columns];

    for (int i = 0; i < this->rows; i++) {
        for (int j = 0; j < this->columns; j++) {
            // filling array with values of matrix cells
            array[i * this->rows + j] = this->matrix[i][j];
        }
    }

    return array;
}

double Matrix::get(int row, int column) {
    return this->matrix[row][column];
}

void Matrix::set(int row, int column, double value) {
    this->matrix[row][column] = value;
}

void Matrix::fill(int type) {
    // matrix is empty => filling it with 0s
    if (this->matrix.empty()) {
        for (int i = 0; i < this->rows; i++) {
            std::vector<double> line;
            this->matrix.push_back(line);

            for (int j = 0; j < this->columns; j++) {
                this->matrix[i].push_back(0);
            }
        }
    }

    switch (type) {
        case 0: // MATRIX_TYPE::ZERO
            for (int i = 0; i < this->rows; i++) {
                for (int j = 0; j < this->columns; j++) {
                    this->matrix[i][j] = 0;
                }
            }
            break;
        case 1: // MATRIX_TYPE::A
            for (int i = 0; i < this->rows; i++) {
                for (int j = 0; j < this->columns; j++) {
                    if (i == j) {
                        this->matrix[i][j] = 1;
                    } else {
                        this->matrix[i][j] = 2;
                    }
                }
            }
            break;
        case 2: // MATRIX_TYPE::B
            for (int i = 0; i < this->rows; i++) {
                for (int j = 0; j < this->columns; j++) {
                    if (i == j) {
                        this->matrix[i][j] = 2;
                    } else {
                        this->matrix[i][j] = 3;
                    }
                }
            }
            break;
        default:
            break;
    }
}

void Matrix::print() {
    for (int i = 0; i < this->matrix.size(); i++) {
        for (int j = 0; j < this->matrix[i].size(); j++) {
            std::cout << this->matrix[i][j] << " ";
        }
        std::cout << std::endl;
    }
}