# LANA DC Circuit Solver


## What is LANA?

The Linear Algebraic Nodal Analysis Algorithm is a circuit analysis algorithm developed by Professor Jeff Anderson at Foothill College (appliedlinearalgebra.com) that bridges together classical nodal analysis and modified nodal analysis to use nonsingular linear system problem to model the electrical behavior of a circuit containing only resistors, DC voltage sources, and DC current sources. 

## My Independent Study Project: C++ console based DC circuit solver



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

## Example: 

### Ideal Circuit Schematic: 
<img width="805" alt="Screen Shot 2021-06-14 at 3 04 09 PM" src="https://user-images.githubusercontent.com/40043538/121965825-7f1b2680-cd22-11eb-9498-590e5b8ac319.png">


### Input: 

      Enter total number of nodes: 9
      Enter total number of circuit variables: 15
      Enter mv, mr, mi
      3 10 2
      Enter the ground node (nodes are counted from 0) : 8
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
      
      
### Output: 

      A_rg

      1 0 0 0 0 -1 0 0 0 
      0 1 0 0 0 -1 0 0 0 
      0 1 -1 0 0 0 0 0 0 
      0 0 0 0 0 1 0 0 -1 
      0 0 0 0 0 1 0 0 -1 
      0 0 0 0 0 1 -1 0 0 
      0 0 0 1 0 0 -1 0 0 
      0 0 0 0 0 0 1 0 -1 
      0 0 0 -1 1 0 0 0 0 
      0 0 0 0 0 0 1 -1 0 


      A_vg

      0 0 1 0 0 -1 0 0 0 
      0 0 -1 1 0 0 0 0 0 
      0 0 0 0 -1 0 0 1 0 


      A_ig

      -1 0 0 0 0 0 0 0 1 
      0 0 0 0 0 0 0 1 -1 


      The voltage drop across all nodes: -------------

      u_0 = 

      -0.868421 
      0.631579 
      3.13158 
      5.63158 
      0.18421 
      -1.86842 
      2.23684 
      5.18421 

      The voltage drop across the resistors: -------------

      v_r = 

      1 
      2.5 
      -2.5 
      -1.86842 
      -1.86842 
      -4.10526 
      3.39474 
      2.23684 
      -5.44737 
      -2.94737 

      The voltage drop across the current sources:  -------------

      v_i = 

      0.868421 
      5.18421 

      The current running through resistors:  -------------

      i_r = 

      1 
      2.5 
      -2.5 
      -1.86842 
      -1.86842 
      -4.10526 
      3.39474 
      2.23684 
      -5.44737 
      -2.94737 

      The current running through the voltage sources :  -------------

      i_v = 

      -11.3421 
      -8.84211 
      -5.44737 
      Program ended with exit code: 0
    
    ```

