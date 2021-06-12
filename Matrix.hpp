//
//  Matrix.hpp
//  getting+back+in+the+groove
//
//  Created by Maria  on 4/1/21.
//

#ifndef Matrix_hpp
#define Matrix_hpp

#include <stdio.h>

#include <iostream>
#include <vector>
using namespace std;

class Matrix {
    
private:
    
    int num_rows;
    int num_cols;
    vector<vector<float>> table;
    

public:
    Matrix();
    Matrix(int num_rows, int num_cols, const float& initial);
    
    
    void set (int row, int col, float val);
    float get (int row, int col) const;
    
    int get_num_rows() const;
    int get_num_cols() const; 

    
    void  delete_column(int col);
    void  delete_row(int row);
    void  print_matrix();
    Matrix  solve_NSLP(Matrix& A, Matrix& L, Matrix& b);

    bool upper_triangular();
   
    Matrix gauss_inverse();
    void Construct_Gauss_Transforms(int col, Matrix& A);
    Matrix identity(int n);
    
    Matrix operator*(Matrix & B);
    Matrix operator+(Matrix & B);
    Matrix operator*(float x);
    
    Matrix& operator=(const Matrix& right);
    Matrix& operator *= ( Matrix& right);

    Matrix forward_sub(Matrix& b);
    Matrix backward_sub(Matrix& y);
    Matrix deflate(Matrix& D, Matrix& A);
    Matrix transpose(Matrix& A);
    
    void swap_rows(Matrix& A, int row1, int row2);
    void div_row(Matrix& A, int x, int row);
    void add_mul_row(Matrix&A, int x, int row1, int row2);
    void sub_rows(Matrix& A, int row1, int row2);
    Matrix augmented_matrix(Matrix& A, Matrix & B);
    Matrix RREF(Matrix& A);
    Matrix particular_solution(Matrix&RREF);
    Matrix solve_GLSP(Matrix&A, Matrix&b);
    
    
    
  
   
};


#endif /* Matrix_hpp */
