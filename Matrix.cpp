//
//  Matrix.cpp
//  getting+back+in+the+groove
//
//  Created by Maria  on 4/1/21.
//

#include "Matrix.hpp"
#include <fstream>
#include <vector>


using namespace std;



// Constructors, Getters and Setters and some other housekeeping --------------------

Matrix::Matrix(int m_num_rows, int m_num_cols, const float& initial){
    table.resize(m_num_rows);
    for(int i = 0; i < table.size(); i++){
        table[i].resize(m_num_cols, initial);
    }
    num_rows = m_num_rows;
    num_cols = m_num_cols;
}


void Matrix::set(int row, int col, float val){
    if(row > num_rows || col > num_cols) cout << "out of bounds" << endl;
    table[row][col] = val;
}


float Matrix::get(int row, int col) const {
    return table[row][col];
}


int Matrix::get_num_rows() const {
    return num_rows;
}


int Matrix::get_num_cols() const {
    return num_cols;
}

void Matrix::print_matrix(){
    
    for(int i = 0; i < num_rows; i++){
        cout << endl;
        for(int j = 0; j < num_cols; j++){
            cout << table[i][j] << " ";
        }
    }
    
    
    cout << endl;
    
}

void Matrix::delete_column(int col){
    
    for (int i = 0; i < num_rows; ++i)
    {
        if (num_cols > col)
        {
            table[i].erase(table[i].begin() + col);
        }
    }
    num_cols--;
}


void Matrix::delete_row(int row){
    
    if (num_cols > row)
    {
        table.erase( table.begin() + row);
    }
    
    num_rows--;
}




// Operator Overloading -----------------------------------------------


// Matrix - Matrix multiplication

Matrix Matrix::operator*(Matrix & B){
    assert (num_cols == B.num_rows);
    Matrix product(num_rows, B.get_num_cols(), 0.0);
    if(num_cols == B.get_num_rows()){
        for(int i = 0; i < num_rows; i++){
            for(int j = 0; j < B.get_num_cols(); j++){
                float temp = 0.0;
                for(int k = 0; k < num_cols; k++){
                    temp += table[i][k] * B.get(k,j);
                }
                product.set(i,j,temp);
            }
        }
    }
   
    return product;
}


// Scalar Matrix Multiplication

Matrix Matrix::operator*(float x){
    Matrix scalar_mult(num_rows, num_cols, 0);
    for(int i = 0; i < num_rows; i++){
        for(int j = 0; j < num_cols; j++){
            float temp = 0.0;
            for(int k = 0; k < num_cols; k++){
                temp += x * table[i][k];
            }
            scalar_mult.set(i,j, temp);
        }
    }
    
    return scalar_mult;
}

Matrix& Matrix::operator=(const Matrix &right){
    num_rows = right.num_rows;
    num_cols = right.num_cols;
    table.resize(num_rows);
    for(int i = 0; i < table.size(); i++){
        table[i].resize(num_cols, 0);
    }
    
    for (int i = 0; i < num_rows; i++) {
        for (int j = 0; j < num_cols; j++) {
            table[i][j] = right.table[i][j];
        }
    }
    return *this;
}

Matrix Matrix::operator+(Matrix & B){
    Matrix add(num_rows, num_cols, 0);
    for(int i = 0; i < num_rows; i++){
        for(int j = 0; j < num_cols; j++){
            float temp = table[i][j] + B.table[i][j];
            add.set(i, j, temp);
            
            
        }
    }
    return add;
    
}





Matrix Matrix::gauss_inverse(){
    Matrix inverse = identity(num_cols);
    for(int i = 0; i < num_rows; i++){
        for(int j = 0; j < num_cols; j++){
            if ( i != j && table[i][j] != 0){
                inverse.set(i,j, -1*table[i][j]);
            }
        }
    }
    return inverse;
}


// boolean function to determine whether a matrix is upper triangular
// returns 1 if upper, 0 otherwise

bool Matrix::upper_triangular(){
    for (int i = 0; i < num_rows; i++)
    for(int j = 0; j < num_cols; j++){
        {
            if(i > j && table[i][j] != 0) return 0;
            
        }
        
    }
    return 1;
}

// constructs identity matrix

Matrix Matrix::identity(int n) {
    Matrix Id(n, n, 0);
    
    for(int i = 0; i < n; i++){
        Id.set(i,i,1);
    }
    return Id;
}



// step 1: Find lower triangular matrix L
/*
 - calculate gauss transforms L1, L2 ... Ln and mutliply matrix A on the right until A is an upper triangular matrix
 - Each gauss transform is constructed as an identity matrix with the entries under the pivot equal to - (entry of A
 you're trying to zero out)/pivot.
 */

void Matrix::Construct_Gauss_Transforms(int pivot_index, Matrix& A){
    float val = 0.000;
    cout << endl <<  "HEY I'm ACTUALLY GOING THROUGH THIS FUNCTION -----------------------------------------" << endl;
    cout <<  endl << "Pivot_index:   " << pivot_index << endl;
    
    for(int i = pivot_index + 1; i < num_rows; i++){
        cout << endl << "value I'm trying to pivot:   " << A.get(i, pivot_index) << endl;
        cout << endl << "Pivot: " << A.get(pivot_index,pivot_index) << endl;
        
        val = -1 * A.get(i, pivot_index)/A.get(pivot_index,pivot_index);
        table[i][pivot_index] = val;
        cout <<  endl << "val:   " << val << endl;
        
    }
    
}


Matrix Matrix::forward_sub(Matrix& b){
    Matrix x(num_rows, 1, 0);
    float temp = b.get(0, 0)/table[0][0];
    x.set(0,0, temp);
    for(int i = 1; i <= num_rows - 1; i++){
        temp = b.get(i, 0);
        
        for(int j = 0; j < i; j++){
            temp -= table[i][j] * x.get(j, 0);
        }
        temp /= table[i][i];
        x.set(i, 0, temp);
    }
    return x;
}


Matrix Matrix::backward_sub(Matrix& y){
    Matrix x(num_rows, 1, 0);
    
    float temp = y.get(num_rows-1, 0)/table[num_rows - 1][num_rows - 1];
    x.set(num_rows - 1, 0, temp);
    
    for(int i = num_rows -2 ; i >= 0; i--){
        temp = y.get(i, 0);
        
        for(int j = num_cols -1 ; j > 0 ; j--){
            temp -= table[i][j] * x.get(j, 0);
        }
        temp /= table[i][i];
        x.set(i, 0, temp);
    }
    return x;
}


Matrix Matrix::solve_NSLP(Matrix& U, Matrix& L, Matrix& b){
    Matrix gauss_transform(num_rows, num_cols, 0.0);
    gauss_transform = identity(num_rows);
    
    Matrix inverse(num_rows, num_cols, 0.0);
    Matrix temp(num_rows, num_cols, 0.0);
    
    int col = 0;
    while(!U.upper_triangular() && col < U.get_num_cols()){
        
        gauss_transform.Construct_Gauss_Transforms(col,U);
        U = gauss_transform*U;
        
        cout << endl << "Gauss transform" << endl;
        gauss_transform.print_matrix();
        inverse = gauss_transform.gauss_inverse();
        
        
        cout << endl << "Inverse" << endl;
        inverse.print_matrix();
        temp = gauss_transform * inverse;
        
        cout << endl << "Verify Inverse:" << endl;
        temp.print_matrix();
        
        cout << endl << "U :" << endl;
        U.print_matrix();
        
        temp = L * inverse;
        L = temp;
        cout << endl << "L" << endl;
        L.print_matrix();
        
        gauss_transform = identity(num_rows);
        
        col++;
        
    }
    
    cout <<  endl << "U: " << endl;
    U.print_matrix();
    
    Matrix x(num_rows, 1, 0);
    Matrix y(num_rows, 1, 0);
    
    
    y =  L.forward_sub(b);
    x = U.backward_sub(y);
    
    cout << endl << "y: " << endl;
    y.print_matrix();
    
    cout << endl << "x: " << endl;
    x.print_matrix();
    
    return x; 
}


Matrix Matrix::transpose(Matrix& A){
    
    Matrix A_T (A.num_cols, A.num_rows, 0);
    for(int i = 0; i < A.num_rows; i++){
        for(int j = 0; j < A.num_cols; j++){
            A_T.set(j, i, A.get(i,j));
        }
    }
    
    return A_T;
    
}

Matrix Matrix::deflate(Matrix& D, Matrix& A){
    
    Matrix A_deflated(A.num_rows, A.num_cols -1,0);
    
    if (D.num_cols != D.num_rows){
        A.delete_row(A.num_rows - 1);
        
    } else{
        A.delete_column(A.num_cols -1);
        
    }
    
    
    
    
    A_deflated = A;
    return A_deflated;
}




