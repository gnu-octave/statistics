## Copyright (C) 2025-26 Jayant Chauhan <0001jayant@gmail.com>
## Copyright (C) 2026 Andreas Bertsatos <abertsatos@biol.uoa.gr>
##
## This file is part of the statistics package for GNU Octave.
##
## This program is free software; you can redistribute it and/or modify it under
## the terms of the GNU General Public License as published by the Free Software
## Foundation; either version 3 of the License, or (at your option) any later
## version.
##
## This program is distributed in the hope that it will be useful, but WITHOUT
## ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
## FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
## details.
##
## You should have received a copy of the GNU General Public License along with
## this program; if not, see <http://www.gnu.org/licenses/>.

## -*- texinfo -*-
## @deftypefn {statistics} {@var{D} =} dummyvar (@var{group})
##
## Create dummy variables.
##
## @code{@var{D} = dummyvar (@var{group})} returns a matrix @var{D} containing
## the dummy variables associated with the grouping variables in @var{group}.
## Each row in @var{D} corresponds to the same observation across all variables
## in @var{group} and each column in @var{D} corresponds to a separate dummy
## variable.  @var{D} is a numeric matrix of @qcode{double} data type containing
## ones and zeros.
##
## The grouping variable in @var{group} can be specified in one of the following
## options:
##
## @itemize
## @item a positive integer vector representing the different group levels in
## the ordered range @code{1:max (@var{group})}.
##
## @item a positive integer matrix with each column corresponding to a separate
## grouping variable and the integer values representing the group levels within
## that grouping variable in the ordered range @code{1:max (@var{group})}.
##
## @item a @qcode{categorical} column vector, in which case the number and order
## of columns in @var{D} correspond to the categories returned by
## @code{categories (@var{group})}.  Categories that are defined but not present
## in @var{group} produce columns of zeros.  Elements of @var{group} that are
## @code{<undefined>} result in rows of @code{NaN} values in @var{D}.
##
## @item a @qcode{cell} array with its elements containing grouping variables
## specified as any of the above options.  Note that all grouping variables in
## the cell array must have the same number of observations.
## @end itemize
##
## @seealso{tabulate, grp2idx, grpstats}
## @end deftypefn

function D = dummyvar (g)

  if (nargin < 1)
    print_usage;
  endif

  if (!isnumeric (g) && !isa (g, 'categorical'))
    error ("Number of elements must not change. Use [] as one of the size inputs to automatically calculate the appropriate size for that dimension.");
  endif

  [nr, nc] = size (g);

  ## --- CATEGORICAL branch ---
  if (isa (g, 'categorical'))

    if (nc != 1 || ndims(g) > 2)
      error ("Categorical grouping variable must have one column.");
    endif

    [idx, gn] = grp2idx (g);
    K = numel (gn);
    D = zeros (nr, K);

    for i = 1:nr
      if (isnan (idx(i)))
        D(i,:) = NaN;
      else
        D(i, idx(i)) = 1;
      endif
    endfor

  ## --- NUMERIC branch ---
  elseif (isnumeric (g))

    if (ndims (g) != 2)
      error ("dummyvar: numeric grouping variable must be either a vector or a matrix.");
    endif
    if (!isempty(g) && (any (g(:) <= 0) || any (g(:) != fix (g(:)))))
      error ("dummyvar: numeric grouping variable must explicitly contain positive integers.");
    endif

    ## Force vector to column vector
    if (isvector (g) && nc > 1)
      g = g(:);
      nr = nc;
      nc = 1;
    endif

    K = max (g, [], 1);
    ## Handle empty input properly
    if (isempty(K))
        D = [];
        return;
    endif
    D = zeros (nr, sum (K));

    ij = 0;
    for i = 1:nc
      tmp = g(:,i);
      for j = 1:K(i)
        ij++;
        D(tmp == j, ij) = 1;
      endfor
    endfor

  endif

endfunction

## Test output
%!test
%! ## numeric grouping vector
%! g = [1; 2; 1; 3; 2];
%! D = dummyvar (g);
%! assert (D, [1, 0, 0; 0, 1, 0; 1, 0, 0; 0, 0, 1; 0, 1, 0]);
%!test
%! g = categorical ({'a'; 'b'; 'a'}, {'a', 'b', 'c'});
%! D = dummyvar (g);
%! cats = categories (g);
%! g_str = cellstr (g);
%! for k = 1:numel (cats)
%!   mask = strcmp (g_str, cats{k});
%!   assert (all (D(mask, k) == 1), true);
%!   assert (all (D(!mask, k) == 0), true);
%! endfor
%!test
%! g = categorical ({'a'; ''; 'b'}, {'a', 'b', 'c'});
%! D = dummyvar (g);
%! assert (D, [1, 0, 0; NaN, NaN, NaN; 0, 1, 0]);
%!test
%! colors = categorical ({'Red'; 'Blue'; 'Green'; 'Red'; 'Green'; 'Blue'});
%! D = dummyvar (colors);
%! assert (D, [0, 0, 1; 1, 0, 0; 0, 1, 0; 0, 0, 1; 0, 1, 0; 1, 0, 0]);
%!test
%! g1 = [1; 1; 1; 1; 2; 2; 2; 2];
%! g2 = [1; 2; 3; 1; 2; 3; 1; 2];
%! D = dummyvar ([g1, g2]);
%! D1 = [1, 0, 1, 0, 0; 1, 0, 0, 1, 0; 1, 0, 0, 0, 1; 1, 0, 1, 0, 0; ...
%!       0, 1, 0, 1, 0; 0, 1, 0, 0, 1; 0, 1, 1, 0, 0; 0, 1, 0, 1, 0];
%! assert (D, D1);
%!test
%! colors = {'red'; 'blue'; 'red'; 'green'; 'yellow'; 'blue'};
%! D = dummyvar (categorical (colors));
%! D1 = [0, 0, 1, 0; 1, 0, 0, 0; 0, 0, 1, 0; 0, 1, 0, 0; 0, 0, 0, 1; 1, 0, 0, 0];
%! assert (D, D1);
%!test
%! g = [1, 2, 1, 2, 1, 3, 2, 1];
%! D = dummyvar (g);
%! D1 = [1, 0, 0; 0, 1, 0; 1, 0, 0; 0, 1, 0; 1, 0, 0; 0, 0, 1; 0, 1, 0; 1, 0, 0];
%! assert (D, D1);
%!test
%! g = categorical ({'a'; 'b'; 'a'}, {'b', 'a', 'c'}, {'b', 'a', 'c'});
%! D = dummyvar (g);
%! assert (D, [0, 1, 0; 1, 0, 0; 0, 1, 0]);
%!test
%! g = categorical ({''; ''; ''}, {'x', 'y'});
%! D = dummyvar (g);
%! assert (D, [NaN, NaN; NaN, NaN; NaN, NaN]);
%!test
%! g = categorical ({'a'; 'a'; 'b'}, {'a', 'b', 'c'});
%! D = dummyvar (g);
%! assert (D, [1, 0, 0; 1, 0, 0; 0, 1, 0]);
%! assert (size (D, 2), 3);
%!error dummyvar ()
%!error<too many inputs> dummyvar (1, 2)
%!error<Categorical grouping variable must have one column> dummyvar (categorical ({'a', 'b'}))
%!error<Number of elements must not change> dummyvar ({'a', 'b'})
%!error<dummyvar: numeric grouping variable must be either a vector or a matrix> dummyvar (ones (3, 3, 3))
%!error<dummyvar: numeric grouping variable must explicitly contain positive integers> dummyvar ([2, 4, 0, 8, 1])
%!error<Number of elements must not change> dummyvar ([true; false])
