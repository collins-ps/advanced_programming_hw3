# Homework 3

Implement a serial, two-dimensional version of the Barnes-Hut n-body algorithm as described in recorded lecture and course notes. Note that even the serial, 2D version can be tricky to implement efficiently, so I have provided some pseudocode and guidance in the recorded lectures in addition to a high-level description of the algorithm. There are three main parts and each requires care to get the details right. Step 1 is to build/add particles to the quadtree; Step 2 is to "summarize" the leaf nodes -- i.e., compute the center of mass at each parent; Step 3 is to scan the quadtree depth first left to right using the distance threshold to calculate the forces.

# Submission Requirements

Include a README with your name, the assignment, a list of any references you used, and a discussion of any shortcomings your code may have. You should also provide a Makefile and instructions for compiling and running your code.
