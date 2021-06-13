# LANA : The Linear Algebraic Nodal Analysis Algorithm 

## What is LANA?

The Linear Algebraic Nodal Analysis Algorithm is a novel circuit analysis algorithm developed by Professor Jeff Anderson at Foothill College (appliedlinearalgebra.com) that bridges together classical nodal analysis and modified nodal analysis to use nonsingular linear system problem to model the electrical behavior of a circuit containing only resistors, dc voltage sources, and dc current sources. 

## The code implementation:

Motivations .. ect ect


## How to input a circuit

##### Data needed:
- [ ] n_g: the number of nodes in your circuit
- [ ] m: the total number of circuit variables 
- [ ] mr: the number of resistors
- [ ] mr: the number of voltage sources
- [ ] mi: the number of current sources

##### Basic Rules: 

**Letter Codes for circuit components** 

Letter Code | Circuit Component
------------ | -------------
R | Resistor
V| DC Voltage Source
I | DC Current Source




#### Anatomy of simplified SPICE netlist:

*note: numbering the components (as shown in example below) is optional, but reccomended for checking your work later.

- Each row provides exactly 4 pieces of information 

      R1 0 5 1
      R2 1 5 1
      R3 1 2 1
      R4 5 8 1
      R5 5 8 1   
      R6 5 6 1
      R7 3 6 1
      R8 6 8 1
      R9 4 3 1
      R10 6 7 1
      V1 2 5 5
      V2 3 2 2.5
      V3 7 4 5
      I1 8 0 1
      I2 7 8 2.5

## What you get 

#### ⚡ Three incidence matricies that completley define the connectivity of the circuit

 -  A_rg: the resistor incidence matrix
 -  A_vg:  the voltage source incidence matrix
 -  A_ig:  the current source incidence matrix
     


#### ⚡ u_0:  The node voltage potentials
#### ⚡ i_r: The currents running through resistors
#### ⚡ i_v: The currents running through the voltage sources
#### ⚡ v_r: The voltage drops across resistors
#### ⚡ v_i: The voltage drops across the current sources

## The math behind the scenes: 

