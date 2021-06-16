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



## The linear algebra behind the scenes; a quick overview of the code implementation:


### (1) Represent the node connectivity of circuit elements as incidence matricies.

Program reads modified SPICE netlist from the console, and generates three incidence matricies -  for the resistors, the voltage sources, and the current sources - that desribe the node connectivity for each circuit element.

   -  A_rg: the resistor incidence matrix
   -  A_vg:  the voltage source incidence matrix
   -  A_ig:  the current source incidence matrix

### (2) Create a general solution to the set of Kirchoff's Voltage Law equations ( a general linear systems problem)




  - solve for the particular solution to the KVL equations: p_g
  - construct the voltage source dependency matrix D, whos columns form the basis of the nullspace of the voltage source   dependency matrix (A_vg).

### (4) Solve the Linear Algebraic Nodal Analysis Equation: 
- add equation 


### (3) Solve for the node voltage potentials
-  add equations
-  
### (4) Solve for the rest of the variables: 
- add equations


 For my independent study project, I took it upon myself to write a C++ program that can solve the complete electrical behavior of *any* DC circuit in a matter of milliseconds, saving time and energy of pen and paper analysis. I was motivated by efficiency - both for humans and for the computer. 
 
I began implementing the LANA algorithm in C++, writing a Matrix class to handle many of the major linear algebraic functions (including the ability to solve nonsingular and general linear systems problems). At first, I tried to directly transcribe the steps of the human algorithm into code. 
 
An essential part of the LANA algorithm is classifying all the nodes as a ground node, essential, nonessential, or supernode (and picking an independent and dependent node). These classifications aim to eliminate any voltage source node dependencies and end up with the minimum number of independent nodes. The next step is to construct a node voltage source dependency matrix D and a particular solution p_g, which will then be used to solve the final linear algebraic nodal analysis problem and then for all the variables. 

It’s undoubtedly an important step and embodies much of what LANA is all about: trimming the fat. 
Yet, my attempts to write a function to classify nodes like a human quickly tumbled into spaghetti code territory. The problem was that there was too much choice involved. When solving circuits on paper, I would approach this step with a hefty dose of inference and thinking ahead - cognitive skills that I could not code into a computer (at least, not yet). 

 Then, in the spirit of LANA, I trimmed the fat. I realized that matrix D and the particular solution p_g formed the general solution to Kirchoff's voltage law equations -  a general linear systems problem. 
 
 <img width="157" alt="Screen_Shot_2021-06-14_at_2 58 22_PM-removebg-preview" src="https://user-images.githubusercontent.com/40043538/121965022-429afb00-cd21-11eb-8450-99ebd090f76b.png">

 Of course, I knew this - *theoretically* - for a while, but it had yet to really click. Alongside working on the LANA circuit solver, I was also self-studying more  advanced topics in Linear Algebra. I really dove into the concepts of vector spaces, subspaces, nullspaces, images - to name a few. I finally understood that the columns of the voltage source dependency matrix D formed the basis of the null space of the voltage source dependency matrix.  
 
When I finally connected the concepts, I had the opportunity to apply them right away. I added a function to transform matrices to reduced row echelon form and find particular solutions to my Matrix class - all the functions necessary to solve general linear systems problems. I then deleted all my old node classification and particular solution functions. 

The differences are shown here: 







