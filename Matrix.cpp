//
//  Matrix.cpp
//  getting+back+in+the+groove
//
//  Created by Maria  on 4/1/21.
//

#include "Matrix.hpp"
#include <fstream>
#include <vector>
#include <algorithm>


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


//  ------------------- ###### OVERLOADED OPERATORS ######  -------------------
   
  



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


//  Matrix Scalar Multiplication

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


//  Matrix - Matrix Addition

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

//  Matrix Assignment Operator

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



//  ------------------- ###### SOLVING NONSINGULAR LINEAR SYSTEMS PROBLEMS ######  -------------------


//                       -------------------  Helper Functions -------------------



//  Construct a Gauss transform

/*  *************************************************************
 *  Given : an integer pivot index, and a Matrix A
 *  Returns: a gauss transform (elementary matrix)
 *  *************************************************************
*/
void Matrix::Construct_Gauss_Transforms(int pivot_index, Matrix& A){
    float val = 0.000;
    for(int i = pivot_index + 1; i < num_rows; i++){
        val = -1 * A.get(i, pivot_index)/A.get(pivot_index,pivot_index);
        table[i][pivot_index] = val;
        
    }
    
}

   
//  Finds the inverse of a Gauss transform

/*  *************************************************************
 *  Given : A Gauss transform matrix
 *  Returns: the inverse of the guass transform
 *  *************************************************************
*/


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


//  Determines whether a matrix is upper triangular

/*  *************************************************************
 *  Given : A matrix
 *  Returns: 1 if uppertriangular, 0 otherwise
 *  *************************************************************
*/

bool Matrix::upper_triangular(){
    for (int i = 0; i < num_rows; i++)
    for(int j = 0; j < num_cols; j++){
        {
            if(i > j && table[i][j] != 0) return 0;
            
        }
        
    }
    return 1;
}


//  Construct an identity matrix

/*  *************************************************************
 *  Given : an integer n, corresponding to the desired dimension
 *  Returns: an n x n identity matrix
 *  *************************************************************
*/

Matrix Matrix::identity(int n) {
    Matrix Id(n, n, 0);
    
    for(int i = 0; i < n; i++){
        Id.set(i,i,1);
    }
    return Id;
}

// Forward substitution algorithm for solving linear systems
/*  *************************************************************
 *  Given : a Matrix (L) in lower triangular form, a vector b
 *  Returns: solution  vector x to Lx = b
 *  *************************************************************
*/


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

// Backwards substitution algorithm for solving linear systems
/*  *************************************************************
 *  Given : a Matrix (U) in upper triangular form, a vector b
 *  Returns: solution  vector x to Ux = b
 *  *************************************************************
*/

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

//        -------------------  Nonsingular Linear Systems Solver -------------------


Matrix Matrix::solve_NSLP(Matrix& U, Matrix& L, Matrix& b){
    Matrix gauss_transform(num_rows, num_cols, 0.0);
    gauss_transform = identity(num_rows);
    Matrix inverse(num_rows, num_cols, 0.0);
    Matrix temp(num_rows, num_cols, 0.0);
    
    

    int col = 0;
    while(!U.upper_triangular() && col < U.get_num_cols()){
        
        gauss_transform.Construct_Gauss_Transforms(col,U);
        U = gauss_transform*U;

        inverse = gauss_transform.gauss_inverse();
        temp = gauss_transform * inverse;
        temp = L * inverse;
        L = temp;
        gauss_transform = identity(num_rows);
        
        col++;
        
    }
    
   
    Matrix x(num_rows, 1, 0);
    Matrix y(num_rows, 1, 0);
    
    
    y =  L.forward_sub(b);
    x = U.backward_sub(y);

    
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


//  ------------------- ###### SOLVING GENERAL LINEAR SYSTEMS PROBLEMS ######  -------------------

//                       -------------------  Helper Functions -------------------


void Matrix::swap_rows(Matrix& A, int row1, int row2){
    
    for(int j = 0; j < num_cols; j++){
        int temp = A.get(row2,j);
        A.set(row2, j, A.get(row1,j));
        A.set(row1, j, temp);
    }
}

void Matrix::div_row(Matrix& A, int x, int row){
    for(int j = 0; j < num_cols; j++)
    A.set(row, j,  x *A.get(row, j));
}

void Matrix::sub_rows(Matrix& A, int row1, int row2){
    for(int j = 0; j < num_cols; j++)
    A.set(row1, j, A.get(row1,j) - A.get(row2,j));
}


void Matrix::add_mul_row(Matrix&A,int x, int row1, int row2){
    for(int j = 0; j < num_cols; j++){
        A.set(row1, j, A.get(row1,j) + x*A.get(row2,j));
    }
}


Matrix Matrix::augmented_matrix(Matrix& A, Matrix & B){
    assert (A.num_rows == B.num_rows);
    int num_cols_to_add = B.num_cols;
    Matrix A_B(A.num_rows, A.num_cols + B.num_cols, 0);
    
    for(int i = 0; i < A.num_rows; i++){
        for(int j = 0; j < A.num_cols + num_cols_to_add; j++){
           if (j < A.num_cols) A_B.set(i,j,A.get(i,j));
           else {
               A_B.set(i,j, B.get(i, j - A.num_cols));
           }
        }
    }
    
    A = A_B;
    return A_B;
    
}

// transforms A into reduced row echelon form
Matrix Matrix::RREF(Matrix& A){
    Matrix U(A.num_rows, A.num_cols, 0);
    int lead = 0;
    int m = 0;
    for ( m = 0; m < num_rows; m++){
        if (num_cols <= lead){
            break;
        }
        int i = m;
        while(A.get(i, lead) == 0){
            i++;
            if (num_rows == i){
                i = m;
                lead++;
                if (num_cols  == lead){
                    break;
                }
            }
        }
        swap(i,m);
        A.div_row(A, A.get(m, lead), m);
        
        
        for(int i = 0; i < num_rows; i++){
            if(i != m){
                A.add_mul_row(A, -1 * A.get(i, lead), i, m);
                
            }
        }
        lead++;
    }
    
    U = A;
    return U;
}

// returns the number of linearlly independent columns of A
int Matrix::rank(Matrix& A){
    int r = 0;
    A = A.RREF(A);

    for(int i = 0; i < num_rows; i++){
        for(int j = 0; j < num_cols; j++){
            if(A.get(i,j) == 1){
                r++;
                break;
            }

        }
 
    }
    
    return r;
}

Matrix Matrix::particular_solution(Matrix &RREF){
    Matrix x (num_cols - 1, 1, 0);
    
    // isolate y : the last column of augmented matrix in RREF form
    
    Matrix y (num_rows, 1, 0);
    for(int i = 0; i < num_rows; i++){
        y.set(i,0, RREF.get(i, num_cols- 1));
    }
    
    // express vector b in terms of the pivot columns of RREF
    for (int i = 0; i  < num_rows; i++){
        for(int j = 0; j < num_cols - 1; j++){
            if(RREF.get(i,j) == 1) x.set(j, 0, y.get(i,0));
        }
    }
    
    
    
    
    return x;
}

Matrix Matrix::solve_GLSP(Matrix&A, Matrix&b){
   
    // Step (1) - create the augmented matrix [A | b ]
    
    Matrix Ay (A.num_rows, A.num_cols + 1, 0);
    Ay = Ay.augmented_matrix(A, b);
    
    // Step (2) - calculate U, the reduced row echelon form
    // of augmented matrix Ay
    
    Matrix U(A.num_rows, A.num_cols + 1, 0);
    U = U.RREF(Ay);
    
    // Step (3) - solve for the particular solution
    
    // solve for particular solution :
    Matrix x( Ay.get_num_cols(), 1, 0);
    x = U.particular_solution(U);
    return x;
  
}
    
    






