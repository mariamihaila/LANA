#include <iostream>
#include <vector>
#include "Matrix.hpp"
using namespace std;



class LANA{
private:
    
    int ng;
    int mr;
    int mv;
    int mi;
    
    int ground_node; 
    Matrix A_rg;
    Matrix A_vg;
    Matrix A_ig;
    
    Matrix A_r0;
    Matrix A_v0;
    Matrix A_i0;
    
    Matrix D_vg;
    Matrix pg;
    
    Matrix V_source;
    Matrix Resistor;
    Matrix I_source;
    
    
public:
    
   
    
    LANA(int ng, int mr, int mv, int mi, int ground_node, const Matrix& V_source, const Matrix& Resistor, const Matrix& I_source,
         const Matrix& A_rg, const Matrix& A_vg, const Matrix& A_ig, const Matrix& A_r0, const Matrix& A_v0, const Matrix& A_i0,      Matrix& D_vg, Matrix&pg);

    void construct_Dvg( vector<int> ind_super_nodes, vector<int> dep_super_nodes,
                        vector<int>ordinary_nodes, vector<int> non_essential_nodes);
    
    void construct_p_g();
    void node_classification ( );
    void construct_K ();
    void solve();
    void this_node_has_to_be_independent(vector<int>& must_be_independent);
};
