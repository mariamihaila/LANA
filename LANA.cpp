
//
//  LANA.cpp
//  getting+back+in+the+groove
//
//  Created by Maria  on 3/31/21.
//

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>


#include "LANA.hpp"
#include  "Matrix.hpp"


LANA::LANA(int ng, int mr, int mv, int mi, int ground_node, const Matrix& V_source,
           const Matrix& Resistor, const Matrix& I_source,const Matrix& A_rg,
           const Matrix& A_vg, const Matrix& A_ig, const Matrix& A_r0,
           const Matrix& A_v0, const Matrix& A_i0, Matrix& D, Matrix&pg,
           const Matrix& u_0, const Matrix& i_r, const Matrix& i_v,
           const Matrix& v_i, const Matrix& v_r) : ng(ng), mr(mr), mv(mv), mi(mi),
                                                   ground_node(ground_node), V_source(V_source),
                                                   Resistor(Resistor), I_source(I_source),
                                                   A_rg(A_rg), A_vg(A_vg), A_ig(A_ig), A_r0(A_r0),
                                                   A_v0(A_v0), A_i0(A_i0), D(D), pg(pg),
                                                   u_0(u_0), i_r(i_r), i_v(i_v), v_i(v_i), v_r(v_r){}



void delete_col(vector<vector<int>> &vect, int col_to_delete){
    for (int i = 0; i < vect.size(); ++i)
    {
        if (vect[i].size() > col_to_delete)
        {
            vect[i].erase(vect[i].begin() + col_to_delete);
        }
    }
}





void LANA ::voltage_source_deflation_matrix(){
    
        // D: the grounded voltage-source deflation matrix
        /*
        -   The columns of D form the basis of the nullspace of Avg
            (the voltage source incidence matrix).
         
        -   In the last step, the column corresponding to the ground node
            is deleted from the D matrix. */
    
    
    // Set U equal to the reduced row echelon form of Avg
    
    Matrix U(0,0,0);
    U = A_vg;
    U = U.RREF(U);
    
    // Take the transpose of U
    U = U.transpose(U);
    
   
    // equivalent to nonpivot col of U
    // equivalent to lin dependent col of Avg)
    // true only if all entries in the row are 0
    bool zero_row = true;
    
    // equivalent to pivot col of U
    // equivalent to lin independent col of Avg)
    // true only ig row contains entry = to 1 && col index is > than in the previous row
    bool lead_found = false;
    
   
    
    // The number of pivot columns of U
    // equivalent to number of independent supernodes
    //int r = A_vg.rank(A_vg);
    
    // nullity: The dimension of the nullspace
    // equals to the number of lin dependent columns of matrix (# of non pivot)
    // equivalent to index of dependent supernode
    int nullity = 0;
    Matrix null_basis(ng, ng, 0);
    
    
    
    // used to check if an entry = to 1 is a lead entry or not
    // set first index to -1 for easy comparison
    int prev_l_index  = -1;
    
    // array to keep track of lead entry row indexes
    vector<int> row_l_index;
    
    int rank = 0;
    for(int i = 0; i < ng; i++){
        lead_found = false;
        zero_row = true; // reset each row
        nullity = i - rank;
        for(int j = 0; j < mv; j++){
            
            if (U.get(i,j) == 1 && !lead_found && j > prev_l_index){
                lead_found = true;
                prev_l_index = j;
                row_l_index.push_back(i);
                zero_row = false;
                rank++;
            } else if (U.get(i,j) != 0){
                zero_row = false;
                D.set(i, nullity, 1);
                D.set(row_l_index[j], nullity, 1);
                }
            }
        
        if (zero_row){
            D.set(i,nullity,1);
        }
    }
    
    D.delete_column(ground_node);
    
    for (int n = 0; n < ng; n++) {
        if ( n > nullity){
            D.delete_column(n);
        }
    }
    D.print_matrix();
}

// a particular solution to the set of voltage-source KVLs, where
// Av_g* pg = V_source.
void LANA::particular_solution(){
    pg = pg.solve_GLSP(A_vg, V_source);
}


void LANA::node_voltage_potentials(int u, Matrix&K, Matrix& D, Matrix& RHS, Matrix& p_0){
    
    Matrix u_solution(u, 1, 0);
    Matrix L(u, u, 0);
    L = L.identity(u);
    u_solution = u_solution.solve_NSLP(K,L,RHS);
    
    u_0 = D * u_solution + p_0;
    cout << endl << "The voltage drop across all nodes: -------------" << endl;
    cout << endl << "u_0 = " << endl;
    u_0.print_matrix();
    
}

void LANA::voltage_drop_across_resistors(){
    v_r = A_r0 * u_0; // Kirchoff's voltage law in node potential form
    cout << endl << "The voltage drop across the resistors: -------------" << endl;
    cout << endl << "v_r = " << endl;
    v_r.print_matrix();
}

void LANA::voltage_drop_across_current_sources(){
    v_i = A_i0*u_0;
    cout << endl << "The voltage drop across the current sources:  -------------" << endl;
    cout << endl << "v_i = " << endl;
    v_i.print_matrix();
}


void LANA::current_through_resistors(Matrix& G){
    i_r = G * v_r; // Ohm's Law in conductance form
    cout << endl << "The current running through resistors:  -------------" << endl;
    cout << endl << "i_r = " << endl;
    i_r.print_matrix();
    
}

void LANA::current_through_voltage_sources(){
    
    Matrix y(ng-1, 1, 0);
    y =  ((A_r0.transpose(A_r0)* i_r));
    y = (A_i0.transpose(A_i0) * I_source) + y;
    y = y * -1;
    
    Matrix i_v(mv, 1, 0);
    Matrix Av0_T(mv, ng - 1, 0);
    Av0_T = A_v0;
    Av0_T = Av0_T.transpose(Av0_T);
    i_v = i_v.solve_GLSP(Av0_T, y);
    cout << endl << "The current running through the voltage sources :  -------------" << endl;
    cout << endl << "i_v = " << endl;
    i_v.print_matrix();
    
}

void LANA::solve() {
    
    voltage_source_deflation_matrix();
    particular_solution();
    
    int u = D.get_num_cols();
    // set D to deflated voltage dependency matrix D by deleting the ground node row
    Matrix D(ng -1, u ,0);
    D = D;
    D.delete_row(ground_node);
    

    // Maximally deflate the circuit equations
    
    Matrix Ar(mr, u, 0);
    Ar = A_r0 * D;
    
    Matrix Av(mv, u, 0);
    Av = A_v0 * D;
    
    Matrix Ai(mi, u, 0);
    Ai = A_i0 * D;
    
    
    // construct the linear-algebraic nodal analysis coefficient matrix: K = Ar^T*G*Ar
    
    Matrix G(mr, mr, 0);
    G = G.identity(mr);
    
    Matrix Ar_T(mr, u, 0);
    Ar_T = Ar;
    Ar_T = Ar_T.transpose(Ar_T);
    Matrix K(u, u, 0);
    
    K = (Ar_T*G) * Ar;

    // construct voltage source forcing terms at the nodes of the circuit: b = -A_r0 * p_0
    
    Matrix b(ng - 1, 1, 0);
    Matrix p_0(ng, 1, 0);
    p_0 = pg;
    p_0.delete_row(ground_node);
    
    b = (A_r0  * p_0);
    b = b * -1;
    
    
    // construct current source forcing terms at the nodes of the circuit: Ai^T * I_source
    Matrix f(u, 1, 0);
    
    Matrix Ai_T(mi, u, 0);
    Ai_T = Ai;
    Ai_T = Ai_T.transpose(Ai_T);
    
    f = (Ai_T * I_source) * -1;
  
 
    // The righthand side to the linear algebraic nodal analysis equations : RHS = Ar^T * G * b - f
    // (The linear algebraic nodal analysis equation in full : (Ar^T * G * Ar) * u = Ar^T * G * b - f
    
    Matrix RHS(u, 1, 0);
    RHS = (Ar_T * G *b) + f;
  
    // Solve for all the circuit variables
    node_voltage_potentials(u, K, D, RHS, p_0);
    voltage_drop_across_resistors();
    voltage_drop_across_current_sources();
    current_through_resistors(G);
    current_through_voltage_sources();
      

    
}


