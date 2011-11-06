## -*- texinfo -*-
## @deftypefn {Function File} ff2n (@var{n})
## Full-factor design with n binary terms.
##
## @seealso {fullfact}
## @end deftypefn

## Author: Paul Kienzle <pkienzle@users.sf.net>
## This program is public domain

function A=ff2n(n)
  A = fullfact (2 * ones (1,n)) - 1;
endfunction
