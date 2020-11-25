#include "Geometry.h"

Matrix::Matrix(int rows, int columns) {
    // initializing sizes variables
    this->rows = rows;
    this->columns = columns;

    // filling all rows*columns-rectangle with 0s
    for (int i = 0; i < this->rows; i++) {
        std::vector<double> vector;
        matrix.push_back(vector);

        for (int j = 0; j < columns; j++) {
            this->matrix[i].push_back(0);
        }
    }
}

Matrix::~Matrix() {
    for (int i = 0; i < this->rows; i++) {
        this->matrix[i].clear();    // setting every row free
    }

    this->matrix.clear();
}

Matrix Matrix::operator=(Matrix other_matrix) {
    // step 0: try of self equating => no need to continue
    if (this == &other_matrix) {
        return *this;
    }

    // step 1: current matrix is bigger than other => removing cells

    // step 1.1: number of current matrix's rows > other's
    while (this->matrix.size() > other_matrix.matrix.size()) {
        // removing elements from last row
        while (this->matrix[this->matrix.size() - 1].size() > 0) {
            this->matrix[this->matrix.size() - 1].pop_back();
        }

        // removing last row
        this->matrix.pop_back();
    }

    // step 1.2: checking every row
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
        for (int j = 0; j < this->columns; j++) {
            this->matrix[i][j] = other_matrix.matrix[i][j];
        }
    }

    return *this;
}

Vector Matrix::operator*(Vector vector) {
    // returns multiplication = matrix * vector

    int number_of_processes, current_process;
    MPI_Comm_size(MPI_COMM_WORLD, &number_of_processes);
    MPI_Comm_rank(MPI_COMM_WORLD, &current_process);

    int lines_per_process = this->columns / number_of_processes;

    // building full vector if it is divided

    Vector vector_full(this->columns);
    vector_full = vector.unite();

    // now we have full vector => starting to count multiplication

    Vector multiplication(this->columns);

    double current_value;

    for (int i = 0; i < lines_per_process; i++) {
        current_value = 0;

        for (int j = 0; j < this->columns; j++) {
            current_value += this->matrix[i][j] * vector_full[j];
        }

        // counting multiplication part
        multiplication.vector[i] = current_value;
    }

    if (current_process != 0) {
        // sending part to 0-process
        MPI_Send(&multiplication.vector[0], lines_per_process, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
        // receiving full from 0-process
        MPI_Recv(&multiplication.vector[0], this->columns, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
    else {
        int beginning = lines_per_process;  // index in full array from which we'll add new received part

        for (int i = 1; i < number_of_processes; i++) {
            // receiving parts from all other processes
            MPI_Recv(&multiplication.vector[beginning], lines_per_process, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            beginning += lines_per_process; // moving index
        }

        // sending full to all processes
        for (int i = 1; i < number_of_processes; i++) {
            MPI_Send(&multiplication.vector[0], this->columns, MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
        }
    }

    return multiplication;
}

double Matrix::get(int row, int column) {
    if (row >= this->rows) {
        std::cout << "ERROR :: MATRIX :: GET :: ROW (" << row << ") >= ROWS OF MATRIX (" << this->rows << ")" << std::endl;
    }
    if (columns >= this->columns) {
        std::cout << "ERROR :: MATRIX :: GET :: COLUMN (" << column << ") >= COLUMNS OF MATRIX (" << this->columns << ")" << std::endl;
    }

    return this->matrix[row][column];
}

void Matrix::set(int row, int column, double value) {
    if (row >= this->rows) {
        std::cout << "ERROR :: MATRIX :: SET :: ROW (" << row << ") >= ROWS OF MATRIX (" << this->rows << ")" << std::endl;
    }
    if (column >= this->columns) {
        std::cout << "ERROR :: MATRIX :: SET :: COLUMN (" << column << ") >= COLUMNS OF MATRIX (" << this->columns << ")" << std::endl;
    }

    this->matrix[row][column] = value;
}

void Matrix::fill(int type, int first_process_line) {
    switch (type) {
        case 0: // MATRIX_TYPE::ZERO
            // result:
            // [0][0][0]...
            // [0][0][0]...
            // [0][0][0]...
            // ...
            for (int i = 0; i < this->rows; i++) {
                for (int j = 0; j < this->columns; j++) {
                    this->matrix[i][j] = 0;
                }
            }
            break;
        case 1: // MATRIX_TYPE::SOLE_A_1
            // result:
            // [2][1][1]...
            // [1][2][1]...
            // [1][1][2]...
            // ...
            for (int i = 0; i < this->rows; i++) {
                for (int j = 0; j < this->columns; j++) {
                    if (i == j - first_process_line) {
                        this->matrix[i][j] = 2;
                    }
                    else {
                        this->matrix[i][j] = 1;
                    }
                }
            }
            break;
        case 2: // MATRIX_TYPE::SOLE_A_2
            // result:
            // [1][0][0]...
            // [1][1][0]...
            // [1][1][1]...
            // ...
            for (int i = 0; i < this->rows; i++) {
                for (int j = 0; j < this->columns; j++) {
                    if (j < i + first_process_line + 1) {
                        this->matrix[i][j] = 1;
                    }
                    else {
                        this->matrix[i][j] = 0;
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