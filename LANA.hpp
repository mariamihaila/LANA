#include <iostream>
#include <vector>
#include "Matrix.hpp"
using namespace std;



class LANA{
private:
    
    int ng;     // Total # of nodes in circuit
    int mr;     // Total # of resistors
    int mv;     // Total # of voltage - sources
    int mi;     // Total # of current - sources
    
    int ground_node;
    
    Matrix V_source; // Voltage drops across V-sources (given)    || dimensions: mv x 1
    Matrix Resistor; // Resistance strengths (given)              || dimensions: mr x 1
    Matrix I_source; // Current running through C-sources (given) || dimensions: mi x 1
    
    // Incidience sub-block matricies
    
    Matrix A_rg;    // resistor sub-block         || dimensions:  mr x ng
    Matrix A_vg;    // voltage - source sub-block || dimensions:  mv x ng
    Matrix A_ig;    // current - source sub_block || dimensions:  mi x ng
    
    
    
    // Grounded incidence sub-block matricies that arise from deleting
    // the column corresponding to the chosen ground node. (n0 = ng - 1)

    Matrix A_r0;    // grounded resistor sub-block           || dimensions:  mr x n0
    Matrix A_v0;    // grounded voltage - source sub-block   || dimensions:  mv x n0
    Matrix A_i0;    // grounded current - source sub-block   || dimensions:  mi x n0
    
    
    // Deflation matrix for entire set of voltage-source KVL equations
    Matrix D;
    
    // Particular solution to entire entire set of voltage-source KVL equations
    Matrix pg;
    
    
    
    // Solution vectors that describe all the circuit behavior:
    
    Matrix u_0;   // voltage drops across nodes(excluding ground) || dimensions: n0 x 1
    
    Matrix i_r;   // current running thorguh resistors            || dimensions: mr x 1
    Matrix i_v;   // current running through V-Sources            || dimensions: mv x 1
    
    Matrix v_i;   // voltage drops across C-Sources               || dimensions: mi x 1
    Matrix v_r;   // voltage drops across resistors               || dimensions: mr x 1

      
    
public:
    LANA(int ng, int mr, int mv, int mi, int ground_node, const Matrix& V_source,
         const Matrix& Resistor, const Matrix& I_source,const Matrix& A_rg, const Matrix& A_vg,
         const Matrix& A_ig, const Matrix& A_r0, const Matrix& A_v0, const Matrix& A_i0,
         Matrix& D, Matrix&pg, const Matrix& u_0, const Matrix& i_r, const Matrix& i_v,
         const Matrix& v_i, const Matrix& v_r);
    

    void voltage_source_deflation_matrix();
    void KVL_particular_solution(); 
    void node_voltage_potentials(int u, Matrix& K, Matrix& D, Matrix& RHS, Matrix & p_0);
    void voltage_drop_across_resistors();
    void voltage_drop_across_current_sources();
    void current_through_resistors(Matrix& G);
    void current_through_voltage_sources();
    void solve();
    
};
