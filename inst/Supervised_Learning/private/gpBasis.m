## Copyright (C) 2026 Andreas Bertsatos <abertsatos@biol.uoa.gr>
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
## @deftypefn {Private Function} gpBasis (@var{X}, @var{basis})
##
## Build the explicit basis matrix of a Gaussian process regression.
##
## @var{X} is @math{NxD} and the returned @var{H} is @math{NxP}, one column per
## basis term, so that the deterministic part of the model is @qcode{H * Beta}.
## @var{basis} is @qcode{'None'}, @qcode{'Constant'}, @qcode{'Linear'},
## @qcode{'PureQuadratic'} as the model stores it, or a function handle taking
## @var{X} and returning the matrix itself.
##
## @qcode{'None'} gives a matrix with no columns, so the model has no
## deterministic part and @qcode{Beta} is empty.
##
## @end deftypefn

function H = gpBasis (X, basis)

  n = rows (X);

  if (is_function_handle (basis))
    H = basis (X);
    if (! isnumeric (H) || ! ismatrix (H) || rows (H) != n)
      error (strcat ("gpBasis: the basis function must return a matrix", ...
                     " with one row per observation."));
    endif
    return;
  endif

  switch (lower (basis))
    case 'none'
      H = zeros (n, 0);
    case 'constant'
      H = ones (n, 1);
    case 'linear'
      H = [ones(n, 1), X];
    case 'purequadratic'
      H = [ones(n, 1), X, X .^ 2];
    otherwise
      error ("gpBasis: unrecognized basis function '%s'.", basis);
  endswitch

endfunction
