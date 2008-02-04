## Copyright (C) 2006 Frederick (Rick) A Niles
##
## This file is intended to be used with this software.
##
## This is free software; you can redistribute it and/or modify it
## under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 2, or (at your option)
## any later version.

## -*- texinfo -*-
## @deftypefn {Function File} {} jsucdf (@var{x}, @var{alpha1}, @var{alpha2})
## For each element of @var{x}, compute the cumulative distribution
## function (CDF) at @var{x} of the Johnson SU distribution with shape parameters
## @var{alpha1} and @var{alpha2}.
##
## Default values are @var{alpha1} = 1, @var{alpha2} = 1.
## @end deftypefn

## Author: Frederick (Rick) A Niles <niles@rickniles.com>
## Description: CDF of the Johnson SU distribution

## This function is derived from normcdf.m

## This is the TeX equation of this function:
## 
## \[ F(x) = \Phi\left(\alpha_1 + \alpha_2 
##    \log\left(x + \sqrt{x^2 + 1} \right)\right) \]
## 
## where \[ -\infty < x < \infty ; \alpha_2 > 0 \] and $\Phi$ is the
## standard normal cumulative distribution function.  $\alpha_1$ and
## $\alpha_2$ are shape parameters.


function cdf = jsucdf (x, alpha1, alpha2)

  if (! ((nargin == 1) || (nargin == 3)))
    usage ("jsucdf (x, alpha1, alpha2)");
  endif

  if (nargin == 1)
    m = 0;
    v = 1;
  endif

  if (!isscalar (alpha1) || !isscalar(alpha2))
    [retval, x, alpha1, alpha2] = common_size (x, alpha1, alpha2);
    if (retval > 0)
      error ("normcdf: x, alpha1 and alpha2 must be of common size or scalar");
    endif
  endif

  one = ones (size (x));
  cdf = stdnormal_cdf (alpha1 .* one + alpha2 .* log (x + sqrt(x.*x + one)));

endfunction
