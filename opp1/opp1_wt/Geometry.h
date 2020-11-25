#ifndef OPP1_W1_GEOMETRY_H
#define OPP1_W1_GEOMETRY_H

#include <iostream>
#include <vector>
#include <cmath>
#include <cstdio>
#include <cstdlib>

// enum class VECTOR_TYPE {ZERO = 0, SOLE_X, SOLE_B_1, SOLE_B_2, SOLE_U};
// enum class MATRIX_TYPE {ZERO = 0, SOLE_A_1, SOLE_A_2};

class Vector {
private:
    std::vector<double> vector;

public:
    Vector(int size); // fills with 0s
    ~Vector();

    Vector operator= (Vector other_vector);
    Vector operator+ (Vector second_vector);
    Vector operator- (Vector second_vector);
    Vector operator* (const double scalar);
    double operator[] (const int index);    // only getter!

    double scalarProduct(Vector second_vector); // returns (*this, second_vector)

    double* get_ptr(int index); // element's pointer getter
    double get(int index);  // element getter
    void set(int index, double value);  // element setter

    void push(double value);    // push_back
    void pop(); // pop_back

    int size(); // returns amount of elements
    double length();    // returns ||*this||

    void fill(int type);  // fills vector according to type
    Vector unite(); // unites parts of vector from different processes and returns it full

    void print();
    void println(); // print + newline

    friend class Matrix;
};

class Matrix {
private:
    int rows;
    int columns;
    std::vector<std::vector<double> > matrix;

public:
    Matrix(int rows, int columns);  // fills matrix with 0s
    ~Matrix();

    Matrix operator=(Matrix other_matrix);
    Vector operator* (Vector vector);   // returns matrix * vector

    double get(int row, int column);    // element getter
    void set(int row, int column, double value);    // element setter
    void fill(int type, int first_process_line);    // fills according to type

    void print();

    friend class Vector;
};

#endif //OPP1_W1_GEOMETRY_H