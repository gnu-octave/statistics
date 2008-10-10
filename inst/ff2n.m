## -*- texinfo -*-
## @deftypefn {Function File} ff2n (@var{n})
## full-factor design with n binary terms.
##
## @seealso {fullfact}
## @end deftypefn

## This program is public domain

function A=ff2n(n) 
  A = fullfact (2 * ones (1,n)) - 1;

