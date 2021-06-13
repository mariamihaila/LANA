# LANA : The Linear Algebraic Nodal Analysis Algorithm 

DRAFT! 

## What is LANA?

The Linear Algebraic Nodal Analysis Algorithm is a circuit analysis algorithm developed by Professor Jeff Anderson at Foothill College (appliedlinearalgebra.com) that bridges together classical nodal analysis and modified nodal analysis to use nonsingular linear system problem to model the electrical behavior of a circuit containing only resistors, dc voltage sources, and dc current sources. 

## The code implementation:

What is it:  
- C++ console based DC circuit solver
- uses Matrix library that I wrote

Motivations .. ect ect

My motivation behind this project was to develop an easy to use program to solve DC circuits analysis problems efficiently, both from the computer's standpoint and from a human standpoint. Easy on the CPU, easy on the human - that's my motto. 


## How to input a circuit

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
   - The polarity (the direction for the resistors can be chosen arbitrarily, while the polarity for the voltage and current sources should correspond to how it is defined by the real circuit element. 
   - ideal circuit schematics usually number nodes starting from 1, but my program follows the computer science convention and counts the first node as node 0. So, when transcribing the connectivity from schematic to netlist, make sure to subtract one from every node.

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


## What you get 

####  Three incidence matricies that completley define the connectivity of the circuit

 -  A_rg: the resistor incidence matrix
 -  A_vg:  the voltage source incidence matrix
 -  A_ig:  the current source incidence matrix
     


####  u_0:  The node voltage potentials
####  i_r: The currents running through resistors
####  i_v: The currents running through the voltage sources
####  v_r: The voltage drops across resistors
####  v_i: The voltage drops across the current sources

## The math behind the scenes: 
-
-
-
-
-

