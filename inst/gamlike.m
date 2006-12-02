## -*- texinfo -*-
## @deftypefn {Function File} {} X = gamlike ([@var{A} @var{B}], @var{R})
## Calculates the negative log-likelihood function for the Gamma
## distribution over vector R, with the given parameters A and B.
## @seealso{gampdf, gaminv, gamrnd, gamfit}
## @end deftypefn

## Written by Martijn van Oosterhout <kleptog@svana.org> (Nov 2006)
## This code is public domain

function res = gamlike(P,K)

a=P(1);
b=P(2);

res = -sum( log( gampdf(K, a, b) ) )

 	  	 
