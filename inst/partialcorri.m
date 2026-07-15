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
## @deftypefn  {statistics} {@var{rho} =} partialcorri (@var{y}, @var{x})
## @deftypefnx {statistics} {@var{rho} =} partialcorri (@var{y}, @var{x}, @var{z})
## @deftypefnx {statistics} {[@var{rho}, @var{pval}] =} partialcorri (@dots{})
## @deftypefnx {statistics} {[@dots{}] =} partialcorri (@dots{}, @var{Name}, @var{Value})
##
## Partial correlation of each response with each predictor, adjusting for the
## remaining predictors.
##
## @code{@var{rho} = partialcorri (@var{y}, @var{x})} returns the sample partial
## correlation coefficients between the columns of the @math{n}-by-@math{p}
## response matrix @var{y} and the columns of the @math{n}-by-@math{q} predictor
## matrix @var{x}.  Element @code{@var{rho}(i,j)} is the partial correlation
## between @code{@var{y}(:,i)} and @code{@var{x}(:,j)}, adjusted for the other
## columns of @var{x} (that is, all columns of @var{x} except the @math{j}-th).
## @var{rho} is a @math{p}-by-@math{q} matrix.
##
## @code{@var{rho} = partialcorri (@var{y}, @var{x}, @var{z})} additionally
## controls for the variables in the @math{n}-by-@math{r} matrix @var{z}, so that
## @code{@var{rho}(i,j)} is adjusted for both the other columns of @var{x} and all
## columns of @var{z}.
##
## @code{[@var{rho}, @var{pval}] = partialcorri (@dots{})} also returns @var{pval},
## a matrix of p-values for testing the hypothesis of no partial correlation
## against the alternative selected by @qcode{'Tail'}.
##
## The @qcode{'Type'}, @qcode{'Rows'}, and @qcode{'Tail'} @var{Name}/@var{Value}
## options are accepted with the same meaning as in @code{partialcorr}.
## @qcode{'Kendall'} is @emph{not} supported and raises an error, as in
## @sc{matlab}.
##
## @seealso{partialcorr, corr, corrcoef, tiedrank}
## @end deftypefn

function [rho, pval] = partialcorri (varargin)

  if (nargin < 2)
    print_usage ();
  endif

  ## Separate the leading numeric matrices from the Name/Value options.
  nmat = 0;
  while (nmat < numel (varargin) && isnumeric (varargin{nmat+1}))
    nmat += 1;
  endwhile

  if (nmat < 2 || nmat > 3)
    error ("partialcorri: invalid number of input matrices.");
  endif

  Yv = varargin{1};
  Xv = varargin{2};
  if (nmat == 3)
    Zc = varargin{3};
  else
    Zc = [];
  endif

  ## Parse and validate options (errors on 'Kendall' under this function's name).
  [Type, Rows, Tail, rem] = parsePairedArguments ({'Type', 'Rows', 'Tail'}, ...
                              {'pearson', 'all', 'both'}, varargin(nmat+1:end));
  if (! isempty (rem))
    error ("partialcorri: unknown or unpaired optional argument.");
  endif
  if (! (ischar (Type) && ischar (Rows) && ischar (Tail)))
    error ("partialcorri: 'Type', 'Rows', and 'Tail' values must be strings.");
  endif
  Type = lower (Type);
  Rows = lower (Rows);
  Tail = lower (Tail);
  if (strcmp (Type, 'kendall'))
    error ("partialcorri: cannot compute Kendall's partial rank correlation.");
  elseif (! any (strcmp (Type, {'pearson', 'spearman'})))
    error ("partialcorri: '%s' is not a valid 'Type'.", Type);
  endif
  if (! any (strcmp (Rows, {'all', 'complete', 'pairwise'})))
    error ("partialcorri: '%s' is not a valid 'Rows' option.", Rows);
  endif
  if (! any (strcmp (Tail, {'both', 'right', 'left'})))
    error ("partialcorri: '%s' is not a valid 'Tail' option.", Tail);
  endif

  if (! (isnumeric (Yv) && isreal (Yv) && ismatrix (Yv)))
    error ("partialcorri: Y must be a real numeric matrix.");
  endif
  if (! (isnumeric (Xv) && isreal (Xv) && ismatrix (Xv)))
    error ("partialcorri: X must be a real numeric matrix.");
  endif
  n = rows (Yv);
  if (rows (Xv) != n)
    error ("partialcorri: X and Y must have the same number of rows.");
  endif
  if (! isempty (Zc) && rows (Zc) != n)
    error ("partialcorri: Z must have the same number of rows as Y.");
  endif

  qx = columns (Xv);
  py = columns (Yv);
  rho = NaN (py, qx);
  pval = NaN (py, qx);

  ## Each predictor column in turn: adjust for the other predictors (and Z), then
  ## reuse partialcorr's cross form (Y vs the single predictor, controlling for
  ## the rest).  This shares the ranking, NaN, and p-value machinery.
  for j = 1:qx
    C = [Xv(:, [1:j-1, j+1:qx]), Zc];
    [rj, pj] = partialcorr (Yv, Xv(:,j), C, ...
                            'Type', Type, 'Rows', Rows, 'Tail', Tail);
    rho(:,j) = rj;
    pval(:,j) = pj;
  endfor

endfunction

%!demo
%! ## Partial correlation of a response with each of two predictors, each
%! ## adjusted for the other predictor.
%! y = [-0.85; 0.33; 1.21; -0.19; 0.74; -1.44; 0.58; 0.02; 1.36];
%! x = [0.42 1.30; 1.15 -0.47; -0.98 0.55; 0.63 2.10; 1.88 -1.02; ...
%!      -0.31 0.86; 0.77 0.14; -1.52 1.77; 0.29 -0.63];
%! rho = partialcorri (y, x)

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

## single response, two predictors, controlling for Z, with p-values
%!test
%! [r, p] = partialcorri (Y, X, Z);
%! assert_equal (r, [-0.7539, -0.8212], 1e-4);
%! assert_equal (p, [0.0189, 0.0066], 1e-4);

## matches the corresponding entries of the symmetric partialcorr matrix
%!test
%! rho = partialcorr (D);
%! assert_equal (partialcorri (Y, X, Z), rho(3,1:2), 1e-12);

## Spearman rank partial correlation
%!test
%! r = partialcorri (Y, X, Z, 'Type', 'Spearman');
%! assert_equal (r, [-0.8381, -0.8677], 1e-4);

## multiple responses, with and without an extra control set
%!test
%! r1 = partialcorri ([D(:,3), D(:,4)], X, D(:,5));
%! assert_equal (r1, [-0.5444, -0.6714; -0.3517, -0.2830], 1e-4);
%! r2 = partialcorri ([D(:,3), D(:,4)], X);
%! assert_equal (r2, [-0.4740, -0.5784; -0.3098, -0.1934], 1e-4);

## input validation
%!error <partialcorri: invalid number of input matrices.> ...
%! partialcorri (ones (5), ones (5), ones (5), ones (5))
%!error <cannot compute Kendall's partial rank correlation.> ...
%! partialcorri (ones (10, 1), ones (10, 2), 'Type', 'Kendall')
%!error <X and Y must have the same number of rows.> ...
%! partialcorri (ones (10, 1), ones (8, 2))
