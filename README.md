# LANA DC Circuit Solver

## What is LANA?

The Linear Algebraic Nodal Analysis Algorithm is a circuit analysis algorithm developed by Professor Jeff Anderson at Foothill College (appliedlinearalgebra.com) that bridges together classical nodal analysis and modified nodal analysis to use nonsingular linear system problem to model the electrical behavior of a circuit containing only resistors, DC voltage sources, and DC current sources. 

## The Code Implementation: C++ console based DC circuit solver


### What kind of circuits are compatible?

 Any electrical circuit with resistors, DC voltage sources, and DC current sources. 

### What does it solve for? 
     
 - u_0:  The node voltage potentials
 - i_r: The currents running through resistors
 - i_v: The currents running through the voltage sources
 - v_r: The voltage drops across resistors
 - v_i: The voltage drops across the current sources

### How do I input a circuit?

##### Data needed:
-  **n_g** : the number of nodes in your circuit
-  **m**: the total number of circuit variables 
-  **mr**: the number of resistors
-  **mr**: the number of voltage sources
-  **mi**: the number of current sources
-  A modified spice netlist describing the connectivity of the circuit


### Anatomy of simplified SPICE netlist

**Letter Codes for circuit elements** 

Letter Code | Circuit Elements
------------ | -------------
R | Resistor
V| DC Voltage Source
I | DC Current Source

 In each netlist, there are as many rows as there are circuit variables.
 Each row of the netlist provides exactly 4 pieces of information: 

     [Letter Code: R / V / I] [node exited] [node entered] [resistance/ voltage strength / current strength]
  
   Notes on convention: 
   
   - The polarity of the resistors can be chosen arbitrarily, while the polarity for the voltage and current sources should correspond to how it is defined by the real circuit element. 
  
   - ideal circuit schematics usually number nodes starting from 1, but my program follows the computer science convention and counts the first node as node 0. So, when transcribing the connectivity from schematic to netlist, make sure to subtract one from every node.

Netlist Example: 


       R 0 5 1
       R 1 5 1
       R 1 2 1
       R 5 8 1
       R 5 8 1   
       R 5 6 1
       R 3 6 1
       R 6 8 1
       R 4 3 1
       R 6 7 1
       V 2 5 5
       V 3 2 2.5
       V 7 4 5
       I 8 0 1
       I 7 8 2.5



## The linear algebra behind the scenes; a quick overview of the code implementation:


### (1) Represent the node connectivity of circuit elements as incidence matricies. 

Program reads modified SPICE netlist from the console, and generates three incidence matricies -  for the resistors, the voltage sources, and the current sources - that desribe the node connectivity for each circuit element.

   -  A_rg: the resistor incidence matrix
   -  A_vg:  the voltage source incidence matrix
   -  A_ig:  the current source incidence matrix

### (2) Create a general solution to the set of Kirchoff's Voltage Law equations: A_vg * u_g = Vv ; a general linear systems problem.
    
  - solve for the particular solution to the KVL equations: p_g
  - construct the voltage source dependency matrix D, whos columns form the basis of the nullspace of the voltage source   dependency matrix (A_vg).

### (3) Solve the Linear Algabraic Nodal Analysis Equation: 
<img width="869" alt="Screen Shot 2021-06-14 at 2 31 04 PM" src="https://user-images.githubusercontent.com/40043538/121962126-2a28e180-cd1d-11eb-975d-1bf6ad287327.png"> 


