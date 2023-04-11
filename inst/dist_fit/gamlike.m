## Copyright (C) 2022 Andreas Bertsatos <abertsatos@biol.uoa.gr>
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
## @deftypefn  {statistics} {@var{nLogL} =} gamlike (@var{params}, @var{R})
##
## Calculates the negative log-likelihood function for the Gamma
## distribution over vector @var{R}, with the given parameters @var{A} and
## @var{B} in a 2-element vector @var{params}.
##
## @seealso{gamcdf, gampdf, gaminv, gamrnd, gamfit}
## @end deftypefn

function nLogL = gamlike (params, R)

  if (nargin != 2)
    print_usage;
  endif

  a = params(1);
  b = params(2);

  nLogL = -sum (log (gampdf (R, a, b)));

endfunction

## Tests
%!error gamlike (1);
%!error gamlike (1, 2, 3);
%!test
%! [nLogL] = gamlike([2, 3], [2, 3, 4, 5, 6, 7, 8, 9]);
%! assert (nLogL, 19.4426, 1e-4);
