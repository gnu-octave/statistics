## Copyright (C) 2026 Andreas Bertsatos <abertsatos@biol.uoa.gr>
##
## This file is part of the statistics package for GNU Octave.
##
## This program is free software; you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 3 of the License, or
## (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program. If not, see <http://www.gnu.org/licenses/>.

## -*- texinfo -*-
## @deftypefn  {statistics} {@var{rho} =} partialcorr (@var{x})
## @deftypefnx {statistics} {@var{rho} =} partialcorr (@var{x}, @var{z})
## @deftypefnx {statistics} {@var{rho} =} partialcorr (@var{x}, @var{y}, @var{z})
## @deftypefnx {statistics} {[@var{rho}, @var{pval}] =} partialcorr (@dots{})
## @deftypefnx {statistics} {[@dots{}] =} partialcorr (@dots{}, @var{Name}, @var{Value})
##
## Linear or rank partial correlation coefficients.
##
## @code{@var{rho} = partialcorr (@var{x})} returns the sample linear partial
## correlation coefficients between pairs of variables in the @math{n}-by-@math{p}
## matrix @var{x}, controlling for the remaining columns of @var{x}.  Each element
## @code{@var{rho}(i,j)} is the partial correlation between @code{@var{x}(:,i)}
## and @code{@var{x}(:,j)}, adjusted for the other @math{p-2} columns.  @var{rho}
## is a symmetric @math{p}-by-@math{p} matrix with ones on the diagonal.
##
## @code{@var{rho} = partialcorr (@var{x}, @var{z})} controls instead for the
## variables in the @math{n}-by-@math{q} matrix @var{z}, returning the
## @math{p}-by-@math{p} partial correlations among the columns of @var{x}.
##
## @code{@var{rho} = partialcorr (@var{x}, @var{y}, @var{z})} returns the
## @math{p1}-by-@math{p2} matrix of partial correlations between the columns of
## the @math{n}-by-@math{p1} matrix @var{x} and the @math{n}-by-@math{p2} matrix
## @var{y}, controlling for @var{z}.  Element @code{@var{rho}(i,j)} is the partial
## correlation between @code{@var{x}(:,i)} and @code{@var{y}(:,j)}.
##
## @code{[@var{rho}, @var{pval}] = partialcorr (@dots{})} also returns @var{pval},
## a matrix of p-values for testing the hypothesis of no partial correlation
## against the alternative selected by @qcode{'Tail'}.
##
## The following @var{Name}/@var{Value} pairs are accepted:
##
## @table @asis
## @item @qcode{'Type'}
## @qcode{'Pearson'} (default) for linear partial correlation, or
## @qcode{'Spearman'} for rank partial correlation (computed on the ranks of the
## data).  @qcode{'Kendall'} is @emph{not} supported and raises an error, as in
## @sc{matlab}.
##
## @item @qcode{'Rows'}
## @qcode{'all'} (default) uses all rows regardless of missing values (any
## @code{NaN} yields a @code{NaN} result); @qcode{'complete'} uses only the rows
## with no missing values across all supplied variables; @qcode{'pairwise'} uses,
## for each computed coefficient, the rows with no missing values among just the
## variables involved in that coefficient.
##
## @item @qcode{'Tail'}
## The alternative hypothesis for @var{pval}: @qcode{'both'} (default, nonzero
## correlation), @qcode{'right'} (greater than zero), or @qcode{'left'} (less than
## zero).
## @end table
##
## The partial correlation is computed by regressing each of the two variables on
## the controlling variables (with an intercept) and correlating the residuals.
## The p-value uses a Student's @math{t} statistic with @math{n - 2 - k} degrees
## of freedom, where @math{k} is the number of controlling variables and @math{n}
## the number of observations used.
##
## @seealso{partialcorri, corr, corrcoef, tiedrank}
## @end deftypefn

function [rho, pval] = partialcorr (varargin)

  if (nargin < 1)
    print_usage ();
  endif

  ## Separate the leading numeric matrices from the Name/Value options.
  nmat = 0;
  while (nmat < numel (varargin) && isnumeric (varargin{nmat+1}))
    nmat += 1;
  endwhile

  if (nmat < 1 || nmat > 3)
    error ("partialcorr: invalid number of input matrices.");
  endif

  args = varargin(1:nmat);
  [Type, Rows, Tail, rem] = parsePairedArguments ({'Type', 'Rows', 'Tail'}, ...
                              {'pearson', 'all', 'both'}, varargin(nmat+1:end));
  if (! isempty (rem))
    error ("partialcorr: unknown or unpaired optional argument.");
  endif
  if (! (ischar (Type) && ischar (Rows) && ischar (Tail)))
    error ("partialcorr: 'Type', 'Rows', and 'Tail' values must be strings.");
  endif
  Type = lower (Type);
  Rows = lower (Rows);
  Tail = lower (Tail);
  if (strcmp (Type, 'kendall'))
    error ("partialcorr: cannot compute Kendall's partial rank correlation.");
  elseif (! any (strcmp (Type, {'pearson', 'spearman'})))
    error ("partialcorr: '%s' is not a valid 'Type'.", Type);
  endif
  if (! any (strcmp (Rows, {'all', 'complete', 'pairwise'})))
    error ("partialcorr: '%s' is not a valid 'Rows' option.", Rows);
  endif
  if (! any (strcmp (Tail, {'both', 'right', 'left'})))
    error ("partialcorr: '%s' is not a valid 'Tail' option.", Tail);
  endif

  ## Assign the roles of the input matrices.
  ##   1 matrix : pairs of columns of X, controlling for the remaining columns.
  ##   2 matrix : pairs of columns of X, controlling for Z.
  ##   3 matrix : columns of X vs columns of Y, controlling for Z.
  A = args{1};
  useRemaining = (nmat == 1);
  square = (nmat <= 2);
  if (nmat == 3)
    B = args{2};
  else
    B = A;
  endif
  if (nmat >= 2)
    Zc = args{end};
  else
    Zc = [];
  endif

  check_matrix ("partialcorr", A, "X");
  n = rows (A);
  if (nmat == 3)
    check_matrix ("partialcorr", B, "Y");
    if (rows (B) != n)
      error ("partialcorr: X and Y must have the same number of rows.");
    endif
  endif
  if (! isempty (Zc))
    check_matrix ("partialcorr", Zc, "Z");
    if (rows (Zc) != n)
      error ("partialcorr: Z must have the same number of rows as X.");
    endif
  endif

  p1 = columns (A);
  p2 = columns (B);

  ## Rows with no missing value across every supplied variable.
  completerows = ! any (isnan ([A, B, Zc]), 2);

  rho = NaN (p1, p2);
  pval = NaN (p1, p2);

  for i = 1:p1
    for j = 1:p2

      if (square && j < i)
        continue;                     # filled below by symmetry
      endif

      a = A(:,i);
      b = B(:,j);
      if (useRemaining)
        C = A(:, setdiff (1:p1, [i, j]));
      else
        C = Zc;
      endif

      [r, nu] = resid_partial (a, b, C, Type, Rows, completerows);
      rho(i,j) = r;
      pval(i,j) = student_pval (r, nu - 2 - columns (C), Tail);

      if (square)
        rho(j,i) = rho(i,j);
        pval(j,i) = pval(i,j);
      endif

    endfor
  endfor

endfunction

## Validate that a matrix input is real, numeric, and 2-D.
function check_matrix (fname, x, label)
  if (! (isnumeric (x) && isreal (x) && ismatrix (x) && ndims (x) == 2))
    error ("%s: %s must be a real numeric matrix.", fname, label);
  endif
endfunction

## Partial correlation of a and b controlling for C, with NaN handling.
function [r, nu] = resid_partial (a, b, C, Type, Rows, completerows)

  n = numel (a);
  switch (Rows)
    case 'all'
      idx = true (n, 1);
    case 'complete'
      idx = completerows;
    case 'pairwise'
      idx = ! any (isnan ([a, b, C]), 2);
  endswitch

  nu = sum (idx);

  ## With 'all', any missing value propagates to a NaN coefficient.
  if (strcmp (Rows, 'all') && any (isnan ([a, b, C])(:)))
    r = NaN;
    nu = n;
    return;
  endif

  a = a(idx);
  b = b(idx);
  C = C(idx,:);

  if (strcmp (Type, 'spearman'))
    a = tiedrank (a);
    b = tiedrank (b);
    for c = 1:columns (C)
      C(:,c) = tiedrank (C(:,c));
    endfor
  endif

  m = numel (a);
  M = [ones(m, 1), C];
  ra = a - M * (M \ a);
  rb = b - M * (M \ b);
  sa = sqrt (sum (ra .^ 2));
  sb = sqrt (sum (rb .^ 2));
  if (sa == 0 || sb == 0)
    r = NaN;
  else
    r = sum (ra .* rb) / (sa * sb);
    r = max (-1, min (1, r));
  endif

endfunction

## p-value from a Student's t statistic for a correlation r with DOF degrees
## of freedom, under the requested tail.
function p = student_pval (r, dof, Tail)
  if (dof <= 0 || isnan (r))
    p = NaN;
    return;
  endif
  t = r .* sqrt (dof ./ (1 - r .^ 2));
  switch (Tail)
    case 'both'
      p = 2 * tcdf (-abs (t), dof);
    case 'right'
      p = tcdf (-t, dof);
    case 'left'
      p = tcdf (t, dof);
  endswitch
endfunction

%!demo
%! ## Partial correlations among four variables, each pair adjusted for the
%! ## other two.
%! x = [0.42 1.30 -0.85 0.11; 1.15 -0.47 0.33 1.82; -0.98 0.55 1.21 -0.34; ...
%!      0.63 2.10 -0.19 0.48; 1.88 -1.02 0.74 0.05; -0.31 0.86 -1.44 1.29; ...
%!      0.77 0.14 0.58 -0.71; -1.52 1.77 0.02 0.94; 0.29 -0.63 1.36 0.37];
%! rho = partialcorr (x)

## shared test data
%!shared D, X, Y, Z
%! D = [ 0.42  1.30 -0.85  0.11  2.04
%!       1.15 -0.47  0.33  1.82 -0.62
%!      -0.98  0.55  1.21 -0.34  0.77
%!       0.63  2.10 -0.19  0.48 -1.15
%!       1.88 -1.02  0.74  0.05  0.39
%!      -0.31  0.86 -1.44  1.29  0.92
%!       0.77  0.14  0.58 -0.71  1.63
%!      -1.52  1.77  0.02  0.94 -0.28
%!       0.29 -0.63  1.36  0.37  0.51
%!       2.01  0.48 -0.77 -1.08  0.14
%!      -0.44  1.05  0.91  0.66 -0.83
%!       0.90 -0.29 -0.36  1.47  1.22];
%! X = D(:,1:2);
%! Y = D(:,3);
%! Z = D(:,4:5);

## single-matrix form: symmetric, unit diagonal, controls for remaining columns
%!test
%! rho = partialcorr (D);
%! assert_equal (diag (rho), ones (5, 1), 1e-12);
%! assert_equal (rho, rho', 1e-12);
%! assert_equal (rho(1,2), -0.8399, 1e-4);
%! assert_equal (rho(1,5), -0.6172, 1e-4);
%! assert_equal (rho(3,4), -0.6763, 1e-4);

## single-matrix form p-values: zero on the diagonal, symmetric off-diagonal
%!test
%! [rho, p] = partialcorr (D);
%! assert_equal (diag (p), zeros (5, 1), 1e-12);
%! assert_equal (p, p', 1e-12);
%! assert_equal (p(1,2), 0.0046, 1e-4);
%! assert_equal (p(1,5), 0.0766, 1e-4);

## two-matrix form: pairs of X controlling for Z, with p-values
%!test
%! [rho, p] = partialcorr (X, Z);
%! assert_equal (rho, [1 -0.5888; -0.5888 1], 1e-4);
%! assert_equal (p, [0 0.0733; 0.0733 0], 1e-4);

## three-matrix form: X vs Y controlling for Z, with p-values
%!test
%! [r, p] = partialcorr (X, Y, Z);
%! assert_equal (r, [-0.2073; -0.5273], 1e-4);
%! assert_equal (p, [0.5656; 0.1173], 1e-4);

## Spearman rank partial correlation
%!test
%! r = partialcorr (X, Y, Z, 'Type', 'Spearman');
%! assert_equal (r, [-0.3079; -0.4984], 1e-4);

## one-sided p-values
%!test
%! [~, pr] = partialcorr (X, Y, Z, 'Tail', 'right');
%! [~, pl] = partialcorr (X, Y, Z, 'Tail', 'left');
%! assert_equal (pr, [0.7172; 0.9414], 1e-4);
%! assert_equal (pl, [0.2828; 0.0586], 1e-4);

## NaN handling: 'all' propagates, 'complete'/'pairwise' delete rows
%!test
%! DN = D;
%! DN(3,2) = NaN;
%! assert (all (isnan (partialcorr (DN))(:)));
%! rc = partialcorr (DN, 'Rows', 'complete');
%! assert_equal (rc(1,2), -0.8311, 1e-4);
%! assert_equal (rc(4,5), -0.6163, 1e-4);

## 'pairwise' differs from 'complete' when entries use different columns
%!test
%! XN = X;
%! XN(3,1) = NaN;
%! rc = partialcorr (XN, Y, Z, 'Rows', 'complete');
%! [rp, pp] = partialcorr (XN, Y, Z, 'Rows', 'pairwise');
%! assert_equal (rc, [-0.0116; -0.5849], 1e-4);
%! assert_equal (rp, [-0.0116; -0.5273], 1e-4);
%! assert_equal (pp, [0.9764; 0.1173], 1e-4);

## input validation
%!error <partialcorr: invalid number of input matrices.> ...
%! partialcorr (ones (5), ones (5), ones (5), ones (5))
%!error <cannot compute Kendall's partial rank correlation.> ...
%! partialcorr (ones (10, 3), 'Type', 'Kendall')
%!error <is not a valid 'Type'.> partialcorr (ones (10, 3), 'Type', 'foo')
%!error <is not a valid 'Rows' option.> partialcorr (ones (10, 3), 'Rows', 'foo')
%!error <X and Y must have the same number of rows.> ...
%! partialcorr (ones (10, 2), ones (8, 1), ones (10, 2))
