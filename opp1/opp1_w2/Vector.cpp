#include "Geometry.h"

Vector::Vector(int size) {
    // filling with 0s
    for (int i = 0; i < size; i++) {
        this->vector.push_back(0);
    }
}

Vector::~Vector() {
    this->vector.clear();
}

Vector Vector::operator=(Vector other_vector) {
    // step 0: try of self equating => no need to continue
    if (this == &other_vector) {
        return *this;
    }

    // step 1: current vector is bigger than other => removing elements
    while (this->vector.size() > other_vector.vector.size()) {
        this->vector.pop_back();
    }

    // step 2: other vector is bigger than current => adding elements (0s for now)
    while (this->vector.size() < other_vector.vector.size()) {
        this->vector.push_back(0);
    }

    // step 3: filling current vector with values from other
    for (int i = 0; i < this->vector.size(); i++) {
        this->vector[i] = other_vector.vector[i];
    }

    return *this;
}

Vector Vector::operator+(Vector second_vector) {
    // returns *this + second_vector

    Vector sum(0);  // empty for now

    int minSize = (this->vector.size() < vector.size()) ? this->vector.size() : vector.size();

    // first minSize elements == sum of vectors
    for (int i = 0; i < minSize; i++) {
        sum.vector.push_back(this->vector[i] + second_vector.vector[i]);
    }

    // pushing tale of biggest vector

    // 1) if current > second
    for (int i = second_vector.vector.size(); i < this->vector.size(); i++) {
        sum.vector.push_back(this->vector[i]);
    }

    // 2) if current < second
    for (int i = this->vector.size(); i < second_vector.vector.size(); i++) {
        sum.vector.push_back(second_vector.vector[i]);
    }

    return sum;
}

Vector Vector::operator-(Vector second_vector) {
    // returns *this - second_vector

    int number_of_processes, current_process, elements_per_process, first_process_element;
    MPI_Comm_size(MPI_COMM_WORLD, &number_of_processes);
    MPI_Comm_rank(MPI_COMM_WORLD, &current_process);

    Vector difference_full(0);   // full difference

    // current vector is bigger than second =>
    // => current is full, second is divided by processes
    if (this->vector.size() > second_vector.vector.size()) {
        elements_per_process = this->vector.size() / number_of_processes;
        first_process_element = current_process * elements_per_process;

        // now we know size of full difference => filling with 0s for now
        for (int i = 0; i < this->vector.size(); i++) {
            difference_full.vector.push_back(0);
        }

        Vector difference_part(0);   // empty for now

        for (int i = 0; i < second_vector.vector.size(); i++) {
            // counting difference part
            difference_part.vector.push_back(this->vector[first_process_element + i] - second_vector.vector[i]);
        }

        if (current_process != 0) {
            // sending part to 0-process
            MPI_Send(&difference_part.vector[0], second_vector.vector.size(), MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
            // receiving full from 0-process
            MPI_Recv(&difference_full.vector[0], this->vector.size(), MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
        else {
            // first elements_per_process cells of full = part from 0-process
            for (int i = 0; i < elements_per_process; i++) {
                difference_full.vector[i] = difference_part[i];
            }

            int beginning = elements_per_process;  // index in full array from which we'll add new received part

            for (int i = 1; i < number_of_processes; i++) {
                // receiving parts from all other processes
                MPI_Recv(&difference_full.vector[beginning], elements_per_process, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

                beginning += elements_per_process; // moving index
            }

            // sending full to all processes
            for (int i = 1; i < number_of_processes; i++) {
                MPI_Send(&difference_full.vector[0], this->vector.size(), MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
            }
        }
    }

    // second vector is bigger than current =>
    // => second is full, current is divided by processes
    if (second_vector.vector.size() > this->vector.size()) {
        elements_per_process = this->vector.size() / number_of_processes;
        first_process_element = current_process * elements_per_process;

        for (int i = 0; i < this->vector.size(); i++) {
            difference_full.vector.push_back(this->vector[i] - second_vector.vector[first_process_element + i]);
        }
    }

    // vectors' sizes are equal
    if (this->vector.size() == second_vector.vector.size()) {
        for (int i = 0; i < this->vector.size(); i++) {
            // just simple difference
            difference_full.vector.push_back(this->vector[i] - second_vector.vector[i]);
        }
    }

    return difference_full;
}

Vector Vector::operator*(const double scalar) {
    // returns *this * scalar

    if (scalar == 1) {
        // no need to continue => the end
        return *this;
    }

    Vector multiplication(0);   // empty for now

    for (int i = 0; i < this->vector.size(); i++) {
        multiplication.vector.push_back(this->vector[i] * scalar);
    }

    return multiplication;
}

double Vector::operator[](const int index) {
    // element getter

    if (index >= this->vector.size()) {
        std::cout << "ERROR :: VECTOR :: OPERATOR[] :: INDEX (" << index << ") >= SIZE OF VECTOR (" << this->vector.size() << ")" << std::endl;
    }

    return this->vector[index];
}

double Vector::scalarProduct(Vector second_vector) {
    // returns (*this, second_vector)

    int number_of_processes, current_process;
    MPI_Comm_size(MPI_COMM_WORLD, &number_of_processes);
    MPI_Comm_rank(MPI_COMM_WORLD, &current_process);

    int elements_per_process = this->vector.size() / number_of_processes;
    int first_process_element = current_process * elements_per_process;

    double sum_full = 0;
    double sum_part = 0;

    for (int i = first_process_element; i < first_process_element + elements_per_process; i++) {
        sum_part += this->vector[i] * second_vector.vector[i]; // counting sum part
    }

    if (current_process != 0) {
        // sending sum part to 0-process
        MPI_Send(&sum_part, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
        // receiving full sum from 0-process
        MPI_Recv(&sum_full, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
    else {
        sum_full += sum_part;   // adding to full sum part that was counted on 0-process

        for (int i = 1; i < number_of_processes; i++) {
            // receiving sum part from other processes
            MPI_Recv(&sum_part, 1, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

            sum_full += sum_part;
        }

        // sending full sum to all processes
        for (int i = 1; i < number_of_processes; i++) {
            MPI_Send(&sum_full, 1, MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
        }
    }

    return sum_full;
}

double* Vector::get_ptr(int index) {
    return &this->vector[index];
}

double Vector::get(int index) {
    if (index >= this->vector.size()) {
        std::cout << "ERROR :: VECTOR :: GET :: INDEX (" << index << ") >= SIZE OF VECTOR (" << this->vector.size() << ")" << std::endl;
    }
    return this->vector[index];
}

void Vector::set(int index, double value) {
    if (index >= this->vector.size()) {
        std::cout << "ERROR :: VECTOR :: SET :: INDEX (" << index << ") >= SIZE OF VECTOR (" << this->vector.size() << ")" << std::endl;
    }
    this->vector[index] = value;
}

void Vector::push(double value) {
    this->vector.push_back(value);
}

void Vector::pop() {
    this->vector.pop_back();
}

int Vector::size() {
    return this->vector.size();
}

double Vector::length() {
    // ||vector|| = sqrt(x_1^2 + x_2^2 + ... + x_size^2)

    int number_of_processes, current_process;
    MPI_Comm_size(MPI_COMM_WORLD, &number_of_processes);
    MPI_Comm_rank(MPI_COMM_WORLD, &current_process);

    double sum_full = 0;
    double sum_part = 0;

    for (int i = 0; i < this->vector.size(); i++) {
        sum_part += this->vector[i] * this->vector[i];
    }

    if (current_process != 0) {
        // sending sum part to 0-process
        MPI_Send(&sum_part, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
        // receiving full sum from 0-process
        MPI_Recv(&sum_full, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
    else {
        sum_full += sum_part;   // adding to full sum part that was counted on 0-process

        for (int i = 1; i < number_of_processes; i++) {
            // receiving sum part from other processes
            MPI_Recv(&sum_part, 1, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            sum_full += sum_part;
        }

        // sending full sum to all processes
        for (int i = 1; i < number_of_processes; i++) {
            MPI_Send(&sum_full, 1, MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
        }
    }

    return std::sqrt(sum_full);
}

void Vector::fill(int type) {
    int number_of_processes, current_process, el;
    MPI_Comm_size(MPI_COMM_WORLD, &number_of_processes);
    MPI_Comm_rank(MPI_COMM_WORLD, &current_process);

    switch (type) {
        case 0: // VECTOR_TYPE::ZERO
            for (int i = 0; i < this->vector.size(); i++) {
                this->vector[i] = 0;
            }
            break;
        case 1: // VECTOR_TYPE::SOLE_X
            for (int i = 0; i < this->vector.size(); i++) {
                this->vector[i] = 0;
            }
            break;
        case 2: // VECTOR_TYPE::SOLE_B_1
            for (int i = 0; i < this->vector.size(); i++) {
                this->vector[i] = number_of_processes * this->vector.size() + 1;
            }
            break;
        case 3: // VECTOR_TYPE::SOLE_B_2
            for (int i = 0; i < this->vector.size(); i++) {
                this->vector[i] = this->vector.size() * current_process + i + 1;
            }
            break;
        case 4: // VECTOR_TYPE::SOLE_U
            for (int i = 0; i < this->vector.size(); i++) {
                this->vector[i] = sin(2 * M_PI * i / (number_of_processes * this->vector.size()));
            }
            break;
        default:
            break;
    }
}

Vector Vector::unite() {
    // unites together parts of vector from different processes

    int number_of_processes, current_process;
    MPI_Comm_size(MPI_COMM_WORLD, &number_of_processes);
    MPI_Comm_rank(MPI_COMM_WORLD, &current_process);

    Vector full(number_of_processes * this->vector.size());

    if (current_process != 0) {
        // sending part to 0-process
        MPI_Send(&this->vector[0], this->vector.size(), MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
        // receiving full from 0-process
        MPI_Recv(&full.vector[0], number_of_processes * this->vector.size(), MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
    else {
        // first this->vector.size() elements of full = part from 0-process
        for (int i = 0; i < this->vector.size(); i++) {
            full.vector[i] = this->vector[i];
        }

        int beginning = this->vector.size();  // index in full array from which we'll add new received part

        for (int i = 1; i < number_of_processes; i++) {
            // receiving parts from all other processes
            MPI_Recv(&full.vector[beginning], this->vector.size(), MPI_DOUBLE, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            beginning += this->vector.size(); // moving index
        }

        // sending full to all processes
        for (int i = 1; i < number_of_processes; i++) {
            MPI_Send(&full.vector[0], number_of_processes * this->vector.size(), MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
        }
    }

    return full;
}

void Vector::print() {
    for (int i = 0; i < this->vector.size(); i++) {
        std::cout << this->vector[i] << " ";
    }
}

void Vector::println() {
    for (int i = 0; i < this->vector.size(); i++) {
        std::cout << this->vector[i] << " ";
    }

    std::cout << std::endl;
}