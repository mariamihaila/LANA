//
//  main.cpp
//  getting+back+in+the+groove
//
//  Created by Maria  on 3/29/21.
//

#include <iostream>
#include <cmath>
#include <iterator>
#include <fstream>
#include <stdlib.h>
#include <sstream>
#include <string>

#include "Matrix.hpp"
#include "LANA.hpp"

#include <vector>

using namespace std;


int main() {
    
    int ng; // # of nodes
    int mv; // # of voltage sources
    int mr; // # of resistors
    int mi; // # of current sources
    int m;  // total # of circuit variables
    int ground_node;
    
    ng = mv = mr = mi = m = ground_node = 0;
  
    
    cout << "Enter total number of nodes: ";
    cin >> ng;
    
    
    cout << "Enter total number of circuit variables: ";
    cin >> m;
    
    cout << "Enter mv, mr, mi" << endl;
    cin >> mv >> mr >> mi;
    
    cout << "Enter the ground node (nodes are counted from 0) : ";
    cin >> ground_node;
    
    
    // declare vectors for known circuit element values
    
    Matrix V_source (mv, 1, 0); // magnitude of voltage drop across V - Source
    Matrix Resistor (mr, 1, 0); // resistance strength
    Matrix I_source (mi, 1, 0); // magnitude of current flowing through current source
    
    // declare incidence matricies that describe node connectivity of each element
   
    Matrix A_rg(mr, ng, 0); // Resistor subblock
    Matrix A_vg(mv, ng, 0); // Voltage source subblock
    Matrix A_ig(mi, ng, 0); // Current source subblock
    

    string type; // Recognized types: V, R, and I
    
    int node1;   // node exited
    int node2;   // node entered
    
    float value;
    
    int mr_count = -1;
    int mv_count = -1;
    int mi_count = -1;

    
    // read the netlist and construct all three incidence matricies
    
    for(int i = 0; i < m; i++){
        
        cin >> type >> node1 >> node2 >> value;
        
        
        if(type.front() == 'R'){
            mr_count++;

            Resistor.set(mr_count, 0, value);
            A_rg.set(mr_count, node1, 1);
            A_rg.set(mr_count, node2, -1);
        }
        
        if(type.front() == 'V'){
            mv_count++;

            V_source.set(mv_count, 0, value);
            A_vg.set(mv_count, node1, 1);
            A_vg.set(mv_count, node2, -1);
            
        }
        
        if(type.front() == 'I'){
            mi_count++;

            I_source.set(mi_count, 0, value);
            A_ig.set(mi_count, node1, 1);
            A_ig.set(mi_count, node2, -1);
        
        }
    }
    
    
    cout << endl << "A_rg" << endl;
    A_rg.print_matrix();
    cout << endl;
    
    cout << endl << "A_vg" << endl;
    A_vg.print_matrix();
    cout << endl;
    
    cout << endl << "A_ig" << endl;
    A_ig.print_matrix();
    cout << endl;
    
    
    cout << endl << "V_source" << endl;
    V_source.print_matrix();
    cout << endl;
  
    
    // declare grounded incidence matricies
    Matrix A_r0(mr, ng -1, 0);
    Matrix A_v0(mv, ng -1, 0);
    Matrix A_i0(mi, ng -1, 0);
    
    // initialize grounded matricies by deleting colum
    // corresponding to the ground node
    
    A_r0 = A_rg;
    A_r0.delete_column(ground_node);
    
    A_v0 = A_vg;
    A_v0.delete_column(ground_node);

    A_i0 = A_ig;
    A_i0.delete_column(ground_node);

    
    Matrix D_vg(ng, ng, 0);
    Matrix pg(ng,1, 0);
    
    
    Matrix u_0 (ng - 1,1,1);   // voltage drops across nodes(excluding ground) || dimensions: n0 x 1
    
    Matrix i_r(mr,1,1);   // current running thorguh resistors            || dimensions: mr x 1
    Matrix i_v(mv,1,1);   // current running through V-Sources            || dimensions: mv x 1
    
    Matrix v_i(mi,1,1);   // voltage drops across C-Sources               || dimensions: mi x 1
    Matrix v_r(mr,1,1);    // voltage drops across resistors               || dimensions: mr x 1
    
    
    
    
    LANA circuit(ng, mr, mv, mi,ground_node, V_source, Resistor, I_source,
                  A_rg, A_vg, A_ig, A_r0, A_v0, A_i0, D_vg, pg, u_0,
                  i_r, i_v, v_i, v_r );
    

    circuit.solve();
     
  

}

