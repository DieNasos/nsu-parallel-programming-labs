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

    Vector difference(0);   // empty for now

    int minSize = (this->vector.size() < second_vector.vector.size()) ? this->vector.size() : second_vector.vector.size();

    // first minSize elements == difference of vectors
    for (int i = 0; i < minSize; i++) {
        difference.vector.push_back(this->vector[i] - second_vector.vector[i]);
    }

    // pushing tale of current vector if it is bigger than second
    for (int i = second_vector.vector.size(); i < this->vector.size(); i++) {
        difference.vector.push_back(this->vector[i]);
    }

    return difference;
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

    double sum = 0;

    for (int i = 0; i < this->vector.size(); i++) {
        sum += this->vector[i] * second_vector[i];
    }

    return sum;
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
    double sum = 0;

    for (int i = 0; i < this->vector.size(); i++) {
        sum += this->vector[i] * this->vector[i];
    }

    return std::sqrt(sum);
}

void Vector::fill(int type) {
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
                this->vector[i] = this->vector.size() + 1;
            }
            break;
        case 3: // VECTOR_TYPE::SOLE_B_2
            for (int i = 0; i < this->vector.size(); i++) {
                this->vector[i] = i + 1;
            }
            break;
        case 4: // VECTOR_TYPE::SOLE_U
            for (int i = 0; i < this->vector.size(); i++) {
                this->vector[i] = sin(2 * M_PI * i / this->vector.size());
            }
            break;
        default:
            break;
    }
}

Vector Vector::unite() {
    return *this;   // no vector dividing => nothing to unite
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