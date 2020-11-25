#pragma once

#include <iostream>
#include <vector>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <limits>
#include <mpi.h>

class Matrix {
private:
    int rows;
    int columns;
    std::vector< std::vector<double> > matrix;

    // sizes of 2D-grid
    int gridSizeX;
    int gridSizeY;

public:
    Matrix(int rows, int columns);
    Matrix(int rows, int columns, int gridSizeX, int gridSizeY);
    Matrix(int rows, int columns, double *source);

    Matrix& operator=(Matrix other_matrix);
    Matrix operator* (Matrix second_matrix);

    void setGridSizes(int x, int y);    // grid sizes setter
    void setContent(double *source);    // content setter
    double* getAsArray();   // returns matrix as array
    double get(int row, int column);    // element getter
    void set(int row, int column, double value);    // element setter
    void fill(int type);    // fills matrix according to type
    void print();   // prints matrix in console
};