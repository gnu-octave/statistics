## Copyright (C) 2019 Nir Krakauer <mail@nirkrakauer.net>
## Based on previous work by Martijn van Oosterhout <kleptog@svana.org>
## originally granted to the public domain.
##
## This file is part of the statistics package for GNU Octave.
##
## This program is free software: you can redistribute it and/or
## modify it under the terms of the GNU General Public License as
## published by the Free Software Foundation, either version 3 of the
## License, or (at your option) any later version.
##
## This program is distributed in the hope that it will be useful, but
## WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
## General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program; see the file COPYING.  If not, see
## <http://www.gnu.org/licenses/>.

## -*- texinfo -*-
## @deftypefn  {statistics} {@var{MLE} =} gamfit (@var{data})
##
## Calculate gamma distribution parameters.
##
## Find the maximum likelihood estimate parameters of the Gamma distribution
## of @var{data}.  @var{MLE} is a two element vector with shape parameter
## @var{A} and scale @var{B}.
##
## This function works by minimizing the value of gamlike for the vector R.
## Just about any minimization function will work, all it has to do is
## minimize for one variable.  Although the gamma distribution has two
## parameters, their product is the mean of the data.  So a helper function
## for the search takes one parameter, calculates the other and then returns
## the value of gamlike.
##
## @seealso{gamcdf, gampdf, gaminv, gamrnd, gamlike}
## @end deftypefn

function res = gamfit (R)

  if (nargin != 1)
    print_usage;
  endif

  avg = mean (R);

  ## Optimize with respect to log(a), since both a and b must be positive
  x = fminsearch (@(x) gamfit_search (x, avg, R), 0);
  a = exp(x);

  b = avg/a;

  res = [a b];

endfunction

## Helper function so we only have to minimize for one variable.
function res = gamfit_search (x, avg, R)
  a = exp (x);
  b = avg / a;
  res = gamlike ([a, b], R);
endfunction


## example data from https://www.real-statistics.com/distribution-fitting/distribution-fitting-via-maximum-likelihood/fitting-gamma-parameters-mle/
%!shared v, res
%! v = [1.2 1.6 1.7 1.8 1.9 2.0 2.2 2.6 3.0 3.5 4.0 4.8 5.6 6.6 7.6];
%! res = gamfit (v);
%!assert (res(1), 3.425, 1E-3);
%!assert (res(2), 0.975, 1E-3);
