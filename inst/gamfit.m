## -*- texinfo -*-
## @deftypefn {Function File} {} [A B] = gamfit (@var{R})
## Finds the maximumlikelihood estimator for the Gamma distribution for R
## @seealso{gampdf, gaminv, gamrnd, gamlike}
## @end deftypefn

## This function works by minimizing the value of gamlike for the vector R.
## Just about any minimization function will work, all it has to do a
## minimize for one variable. Although the gamma distribution has two
## parameters, their product is the mean of the data. so a helper function
## for the search takes one parameter, calculates the other and then returns
## the value of gamlike.

## Note: Octave uses the inverse scale parameter, which is the opposite of
## Matlab. To work for Matlab, value of b needs to be inverted in a few
## places (marked with **)

## Written by Martijn van Oosterhout <kleptog@svana.org> (Nov 2006)
## This code is public domain

function res = gamfit(R)

avg = mean(R);

# This can be just about any search function. I choose this because it
# seemed to be the only one that might work in this situaition...
a=nmsmax( @gamfit_search, 1, [], [], avg, R );

b=a/avg;      # **

res=[a 1/b];

# Helper function so we only have to minimize for one variable. Also to
# inverting the output of gamlike, incase the optimisation function wants to
# maximize rather than minimize.

function res = gamfit_search( a, avg, R )

b=a/avg;      # **

res = -gamlike([a 1/b], R);

 	  	 
