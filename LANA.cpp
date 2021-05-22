
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


LANA::LANA(int ng, int mr, int mv, int mi, int ground_node, const Matrix& V_source, const Matrix& Resistor, const Matrix& I_source,
           const Matrix& A_rg, const Matrix& A_vg, const Matrix& A_ig, const Matrix& A_r0, const Matrix& A_v0, const Matrix& A_i0,      Matrix& D_vg, Matrix&pg) : ng(ng), mr(mr), mv(mv), mi(mi), ground_node(ground_node), V_source(V_source),
Resistor(Resistor), I_source(I_source),
A_rg(A_rg), A_vg(A_vg), A_ig(A_ig), A_r0(A_r0), A_v0(A_v0), A_i0(A_i0), D_vg(D_vg), pg(pg) {}



int _in_ (int key, vector <int> vect){
    int index = -1;
    for(int i = 0; i < vect.size(); i++){
        if(vect[i] == key) index = i;
    }
    return index;
}


void delete_col(vector<vector<int>> &vect, int col_to_delete){
    for (int i = 0; i < vect.size(); ++i)
    {
        if (vect[i].size() > col_to_delete)
        {
            vect[i].erase(vect[i].begin() + col_to_delete);
        }
    }
}



void LANA ::construct_Dvg(vector<int> ind_super_nodes, vector<int> dep_super_nodes, vector<int> ordinary_nodes, vector<int> non_essential_nodes){
    //int ground_index = ng - 1;

    for (int m = 0; m < ng; m++){
        for (int n = 0; n < ng; n++){
            
            
            // if node n is the independent node from the supernode set or if it is ordinary (these are the independent nodes
            // to be calculated in the end)
            
            if ( _in_(n, ordinary_nodes) != -1 ){
                
                D_vg.set(n, n, 1);
            }
            
            if (_in_(n, ind_super_nodes) != -1){
                D_vg.set(n, n ,1);
                int index = _in_ (n, ind_super_nodes);
                int x = dep_super_nodes[index];
                D_vg.set(x, n, 1);
                
                
            }
            // if node n is an dependent from the supernode set
            
            if (_in_(n, dep_super_nodes) != -1){
                int index = _in_(n, dep_super_nodes);
                int x = ind_super_nodes[index];
                D_vg.set(n, x, 1);
            }
        }
        
    }
    
    int num_cols_deleted = 0;
    for (int n = 0; n < ng; n++){
        if (_in_(n, dep_super_nodes) != -1 || n == ground_node  || _in_(n, non_essential_nodes) != -1){
            D_vg.delete_column(n - num_cols_deleted++);
        }
    }
    
    cout << endl << "Node dependency matrix D_Vg:" << endl;
    D_vg.print_matrix();
    
}



void LANA::construct_p_g(vector<int> dep_super_nodes, vector<int> non_essential_nodes){
    
    for(int i = 0; i < dep_super_nodes.size(); i++){
        int node_var = dep_super_nodes[i];
        cout << endl << " dependent node car " << node_var << endl;
        for(int j = 0; j < mv; j++){
            if (A_vg.get(j,node_var) != 0){
                pg.set(node_var,0, A_vg.get(j,node_var));
                pg.set(node_var,0, V_source.get(j,0) * pg.get(node_var,0));
        }
    }
    
    }
    
    for(int i = 0; i < non_essential_nodes.size(); i++){
        int node_var = non_essential_nodes[i];
        cout << endl << " non essential node car " << node_var << endl;
        for(int j = 0; j < mv; j++){
            if (A_vg.get(j,node_var) != 0){
                pg.set(node_var,0, A_vg.get(j,node_var));
                pg.set(node_var,0, V_source.get(j,0) * pg.get(node_var,0));
        }
        
    }
    }
    

    
    cout << endl << "particular solution pg" << endl;
    for(int i = 0; i < pg.get_num_rows(); i++){
        cout << pg.get(i,0) << endl;
    }
    
}


void LANA::this_node_has_to_be_independent(vector<int>& must_be_independent){
    int frequency[ng];
    for(int i = 0; i < ng; i++) frequency[i] = 0;
    
    for(int i = 0; i < mv; i++){
        for(int j = 0; j < ng; j++){
            if(A_vg.get(i,j) != 0){
                frequency[j]++;
                if (frequency[j] > 1) must_be_independent.push_back(j);
            }
        }
    }
}


void LANA::node_classification(){
    
    vector<int>  ind_super_nodes;
    vector <int> dep_super_nodes;
    vector<int>  ordinary_nodes;
    vector<int>  non_essential_nodes;
    int generalized_node[2] = {-1, -1};
    vector<int> must_be_independent;
    this_node_has_to_be_independent(must_be_independent);
    
    for (int i = 0; i < mv; i++){
        int count = -1;
        for (int j = 0; j < ng; j++){
            
            if (A_vg.get(i,j) != 0){
                count++;
                generalized_node[count] = j;
            }
            
        }
        
        if (generalized_node[1] == ground_node){
            non_essential_nodes.push_back(generalized_node[0]);
            ordinary_nodes.push_back(generalized_node[1]);
        } else {
        // Case1: the first node has already been classed as independent
        // so the second one is automatically dependent
        if ((_in_(generalized_node[0], ind_super_nodes) != -1) || (_in_(generalized_node[0], must_be_independent) != -1)){
            ind_super_nodes.push_back(generalized_node[0]);
            dep_super_nodes.push_back(generalized_node[1]);
        }
        
        // Case2: the second node has already been classed as independent
        // so the first one is automatically dependent
        
        else {
            ind_super_nodes.push_back(generalized_node[1]);
            dep_super_nodes.push_back(generalized_node[0]);
        }
        }
    }
    
    for(int i = 0; i < ng; i++){
        if (_in_(i, ind_super_nodes) == -1 && _in_(i, dep_super_nodes) == -1)
            ordinary_nodes.push_back(i);
    }
    
    
    int ground_index_dependent = _in_(ground_node, dep_super_nodes);
    if ( ground_index_dependent != -1){
        non_essential_nodes.push_back(ind_super_nodes[ground_index_dependent]);
        dep_super_nodes.erase(dep_super_nodes.begin() + ground_index_dependent);
        ind_super_nodes.erase(ind_super_nodes.begin() + ground_index_dependent);
    }
    
    
    cout << "dependent nodes from supernode: " << endl;
    for(int i = 0; i < dep_super_nodes.size(); i++){
        cout << dep_super_nodes[i] << " ";
    }
    cout << endl;
    
    cout << "independent nodes from supernode: " << endl;
    for(int i = 0; i < ind_super_nodes.size(); i++){
        cout << ind_super_nodes[i] << " ";
    }
    cout << endl;
    
    cout << "ordinary nodes:  " << endl;
    for(int i = 0; i < ordinary_nodes.size(); i++){
        cout << ordinary_nodes[i] << " ";
    }
    cout << endl;
    
    cout << "non essential nodes:  " << endl;
    for(int i = 0; i < non_essential_nodes.size(); i++){
        cout << non_essential_nodes[i] << " ";
    }
    cout << endl;
    
    construct_Dvg(ind_super_nodes, dep_super_nodes, ordinary_nodes, non_essential_nodes);
    construct_p_g(dep_super_nodes, non_essential_nodes);
    cout << endl;
    
}


void LANA::solve() {
    node_classification();
    
    // u - the number of independent nodes
    int u = D_vg.get_num_cols();
    
    // set D to deflated voltage dependency matrix D_vg by deleting the ground node row
    Matrix D(ng -1, u ,0);
    D = D_vg;
    D.delete_row(ground_node);

    cout << endl << "D" << endl;
    D.print_matrix();
    
    // ground the blocks of the incidence matrix
    
    Matrix Ar(mr, u, 0);
    Ar = A_r0 * D;
    
    Matrix Av(mv, u, 0);
    Av = A_v0 * D;
    
    Matrix Ai(mi, u, 0);
    Ai = A_i0 * D;
    
    
    // construct matrix K = Ar^T*G*Ar
    
    Matrix G(mr, mr, 0);
    G = G.identity(mr);
    
    Matrix Ar_T(mr, u, 0);
    Ar_T = Ar;
    Ar_T = Ar_T.transpose(Ar_T);
    Matrix K(u, u, 0);
    
    K = (Ar_T*G) * Ar;
    
    cout << endl << "K " << endl;
    K.print_matrix();
    
    // construct matrix b  = Ar0 * p0
    
    Matrix b(ng - 1, 1, 0);
    Matrix p_0(ng, 1, 0);
    p_0 = pg;
    p_0.delete_row(ground_node);

    b = (A_r0  * p_0);
    b = b * -1;
    
    cout << endl << "b" << endl;
    b.print_matrix();
    
    //construct vector f = Ai^T * I_source
    Matrix f(u, 1, 0);
    
    Matrix Ai_T(mi, u, 0);
    Ai_T = Ai;
    Ai_T = Ai_T.transpose(Ai_T);
    
    f = (Ai_T * I_source) * -1;
    cout << endl << "f" << endl;
    f.print_matrix();
    A_i0.print_matrix();
    A_v0.print_matrix();
    
    //Ai.print_matrix();
    
    //RHS - Ar^T * G * b - f
    
    Matrix RHS(u, 1, 0);
    
    RHS = (Ar_T * G *b) + f;
    cout << endl << "RHS" << endl;
    RHS.print_matrix();
    
    A_r0.print_matrix();
    
    
    Matrix u_solution(u, 1, 0);
    Matrix L(u, u, 0);
    L = L.identity(u);
    u_solution = u_solution.solve_NSLP(K,L,RHS);
    
    Matrix u_0(ng -1, 1, 0);
    u_0 = D * u_solution + p_0;
    cout << endl << "u_0" << endl;
    u_0.print_matrix();
    
    
    
    
    
    
}


