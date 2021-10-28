# Homework 3

Implement a serial, two-dimensional version of the Barnes-Hut n-body algorithm as described in recorded lecture and course notes. Note that even the serial, 2D version can be tricky to implement efficiently, so I have provided some pseudocode and guidance in the recorded lectures in addition to a high-level description of the algorithm. There are three main parts and each requires care to get the details right. Step 1 is to build/add particles to the quadtree; Step 2 is to "summarize" the leaf nodes -- i.e., compute the center of mass at each parent; Step 3 is to scan the quadtree depth first left to right using the distance threshold to calculate the forces.

# Submission Requirements

Include a README with your name, the assignment, a list of any references you used, and a discussion of any shortcomings your code may have. You should also provide a Makefile and instructions for compiling and running your code.

# MPCS 51100 HW3
**Name**: Phoebe Collins, UCID: 12277438

## References
I did not use any references other than class materials.

## Discussion on program and performance
I did not observe any compilation or run-time errors. The Valgrind report is clean.

With theta set to 0,

x and y values converge to a precision of 0.0001 up to approximately 920 particles for one iteration.

x and y values converge to a precision of 0.001 up to approximately 1350 particles for one iteration.

x and y values converge to a precision of 0.01 up to approximately 4400 particles for one iteration.

The algorithm's accuracy declines as the number of particles increases because of the non-associativity of floating point arithmetic and the chaos of the n-body problem which means that even small differences in the values could grow exponentially with time. 

## Compiling and running
`make nbody_bh` and then `./nbody_bh` which uses the default 3000 particles or `./nbody_bh x` where x is an integer representing the number of particles. Finally, run `make clean`.
