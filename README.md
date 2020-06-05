# SchurUtilities
Various utility scripts to calculate Schur polynomials
These were mainly used to experiment with the time complexity of the different methods of calculating Schur polynomials

![Length of partition vs time](/plots/length_plot.png)
![Sum of partition vs time](/plots/n_plot.png)

Two main methods were compared: the use of constraint programming to generate all semistandard Young tableaux. 
##### METHOD 1: Constraint Programming
Since Schur polynomials can be defined as a sum over all semistandard Young tableaux of a certain partition, it seems natural to generate all such tableau and then taking the sum. We use [Or Tools from Google](https://developers.google.com/optimization) to generate all possible tableaux.
Note that the number of statements needed to describe a tableaux of a partition of n is
  * at minimum, 0 (for a partition of the form (1,1,....,1))
  * at maximum, n-1 (for a partition of the form (n))
  * on average, seems to increase linearly in n
  U
