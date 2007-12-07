## ff2n(n)
##   full-factor design with n binary terms.
##
## see also: fullfact

## This program is public domain

function A=ff2n(n) 
  A = fullfact (2 * ones (1,n)) - 1;

