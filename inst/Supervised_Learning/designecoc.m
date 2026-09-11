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
## @deftypefn  {statistics} {@var{M} =} designecoc (@var{K}, @var{name})
## @deftypefnx {statistics} {@var{M} =} designecoc (@dots{}, @qcode{'NumTrials'}, @var{n})
##
## Coding design matrix for an error correcting output codes model.
##
## @code{@var{M} = designecoc (@var{K}, @var{name})} returns the coding design
## for @var{K} classes named by @var{name}.  @var{M} is a @math{KxL} matrix of
## -1, 0 and +1 with one row per class and one column per binary learner: a
## learner is trained to tell the classes marked +1 in its column from those
## marked -1, and a class marked 0 takes no part in it.
##
## @var{K} must be an integer of at least 2.  @var{name} must be one of:
##
## @multitable @columnfractions 0.23 0.17 0.60
## @headitem Design @tab Columns @tab Description
##
## @item @qcode{'onevsone'} @tab @math{K(K-1)/2} @tab One column per pair of
## classes, the earlier class +1 and the later -1, pairs taken in order.
##
## @item @qcode{'onevsall'} @tab @math{K} @tab One column per class, that
## class +1 and every other -1.
##
## @item @qcode{'binarycomplete'} @tab @math{2^(K-1)-1} @tab Every way of
## splitting the classes into two non-empty groups, with the first class
## always +1.
##
## @item @qcode{'ternarycomplete'} @tab @math{(3^K-2^(K+1)+1)/2} @tab Every
## way of splitting into two non-empty groups while leaving any classes out.
## It grows fast: 28501 columns at @math{K = 10}.
##
## @item @qcode{'ordinal'} @tab @math{K-1} @tab Column @math{j} separates the
## first @math{j} classes from the rest, for classes that are ordered.
##
## @item @qcode{'denserandom'} @tab about @math{10log_2 K} @tab Random -1 and
## +1, no class left out.
##
## @item @qcode{'sparserandom'} @tab about @math{15log_2 K} @tab Random -1, 0
## and +1, a class left out of a column with probability @math{0.5}.
## @end multitable
##
## @code{@var{M} = designecoc (@dots{}, @qcode{'NumTrials'}, @var{n})} draws
## @var{n} random designs and keeps the one whose rows are furthest apart,
## which is what makes a random design correct errors.  The default is 10000.
## It is accepted but does nothing for the five designs that are not random.
##
## @math{K = 2} gives the single column @code{[-1; 1]} whatever the design,
## there being only one way to tell two classes apart.
##
## @subheading Deviation from MATLAB
##
## The two random designs cannot be reproduced from MATLAB and neither can
## their width.  Measured on R2024a: five runs at @math{K = 10} gave
## @qcode{'denserandom'} 38, 38, 38, 40 and 40 columns and
## @qcode{'sparserandom'} 57, 55, 55 and 54, so the number of columns varies
## between runs of MATLAB itself.  The five other designs are exact.
##
## @seealso{fitcecoc, ClassificationECOC}
## @end deftypefn

function M = designecoc (K, name, varargin)

  ## Input validation
  if (nargin < 2)
    error ("designecoc: too few input arguments.");
  endif
  if (mod (numel (varargin), 2) != 0)
    error ("designecoc: name-value arguments must be in pairs.");
  endif
  if (! (isnumeric (K) && isscalar (K) && isreal (K) && K == fix (K)
         && K >= 2))
    error ("designecoc: K must be an integer of at least 2.");
  endif
  if (! (ischar (name) && isrow (name)))
    error ("designecoc: NAME must be a character vector.");
  endif

  NumTrials = 10000;
  for i = 1:2:numel (varargin)
    if (! (ischar (varargin{i}) && isrow (varargin{i})))
      error ("designecoc: invalid name-value argument.");
    endif
    switch (tolower (varargin{i}))
      case 'numtrials'
        NumTrials = varargin{i+1};
        if (! (isnumeric (NumTrials) && isscalar (NumTrials)
               && isreal (NumTrials) && NumTrials == fix (NumTrials)
               && NumTrials >= 1))
          error (strcat ("designecoc: 'NumTrials' must be a positive", ...
                         " integer."));
        endif
      otherwise
        error ("designecoc: invalid name-value argument.");
    endswitch
  endfor

  designs = {'onevsone', 'onevsall', 'binarycomplete', 'ternarycomplete', ...
             'ordinal', 'denserandom', 'sparserandom'};
  design = tolower (name);
  if (! any (strcmp (design, designs)))
    error ("designecoc: '%s' is not a coding design.", name);
  endif

  ## Two classes are told apart one way however the design is named, and
  ## every rule below would otherwise give a wider matrix saying the same
  ## thing twice.  Measured on R2024a.
  if (K == 2)
    M = [-1; 1];
    return;
  endif

  switch (design)

    case 'onevsone'
      L = K * (K - 1) / 2;
      M = zeros (K, L);
      c = 0;
      for i = 1:K-1
        for j = i+1:K
          c++;
          M(i, c) = 1;
          M(j, c) = -1;
        endfor
      endfor

    case 'onevsall'
      M = 2 * eye (K) - 1;

    case 'ordinal'
      M = ones (K, K - 1);
      for j = 1:K-1
        M(1:j, j) = -1;
      endfor

    case 'binarycomplete'
      ## The first class is +1 throughout, so a column is a subset of the
      ## other K-1 and the columns run from all but the last down to none,
      ## the second class being the most significant bit.
      L = 2 ^ (K - 1) - 1;
      M = ones (K, L);
      for c = 1:L
        m = L - c;
        for r = 2:K
          M(r, c) = 2 * bitget (m, K - r + 1) - 1;
        endfor
      endfor

    case 'ternarycomplete'
      ## A column is taken only one way round, so the last class taking part
      ## is always +1 and everything after it sits out.  For each position of
      ## that class, the classes before it run through -1, 0 and +1 with the
      ## first moving fastest, and a column with no -1 in it is not a split.
      cols = {};
      for p = 2:K
        digits = p - 1;
        for t = 0:(3 ^ digits - 1)
          v = zeros (digits, 1);
          rem_ = t;
          for d = 1:digits
            v(d) = mod (rem_, 3) - 1;
            rem_ = floor (rem_ / 3);
          endfor
          if (! any (v == -1))
            continue;
          endif
          cols{end+1} = [v; 1; zeros(K - p, 1)];
        endfor
      endfor
      M = cell2mat (cols);

    otherwise
      M = randomDesign (K, design, NumTrials);

  endswitch

endfunction

## The two random designs.  Columns are drawn at random, those that say
## nothing or repeat another are dropped, and the draw is repeated so that the
## design kept is the one whose rows stand furthest apart, which is what lets
## a random design correct an error.
function M = randomDesign (K, design, NumTrials)

  if (strcmp (design, 'denserandom'))
    if (K <= 5)
      warning (strcat ("designecoc: use 'binarycomplete' instead of", ...
                       " 'denserandom' for 5 or fewer classes."));
    endif
    L = ceil (10 * log2 (K));
    p0 = 0;
  else
    if (K <= 4)
      warning (strcat ("designecoc: use 'ternarycomplete' instead of", ...
                       " 'sparserandom' for 4 or fewer classes."));
    endif
    L = ceil (15 * log2 (K));
    p0 = 0.5;
  endif

  M = [];
  best = -Inf;
  for t = 1:NumTrials
    C = sign (rand (K, L) - 0.5);
    if (p0 > 0)
      C(rand (K, L) < p0) = 0;
    endif
    C = keepUsableColumns (C);
    if (isempty (C) || any (all (C >= 0, 2)) || any (all (C <= 0, 2)))
      continue;
    endif
    d = minRowDistance (C);
    if (d > best)
      best = d;
      M = C;
    endif
  endfor

  if (isempty (M))
    error (strcat ("designecoc: could not draw a usable '%s' design for", ...
                   " %d classes."), design, K);
  endif

endfunction

## A column with no +1 or no -1 trains nothing, and one that repeats another
## up to sign trains the same learner twice.  Taking each column with its
## first class on the +1 side makes a column and its mirror the same row, so
## the repeats fall out of one call rather than a comparison of every pair.
function C = keepUsableColumns (C)

  C = C(:, any (C > 0, 1) & any (C < 0, 1));
  if (isempty (C))
    return;
  endif
  [~, first] = max (C != 0, [], 1);
  lead = C(sub2ind (size (C), first, 1:columns (C)));
  [~, keep] = unique ((C .* lead)', 'rows', 'first');
  C = C(:, sort (keep));

endfunction

## How far apart the two closest rows are, counting a position where one row
## sits out as half a step, which is how an error correcting code is scored.
## Written as one cross product: the distance between rows a and b is
## (L - C(a,:) * C(b,:)') / 2, so the whole table is (L - C*C') / 2.
function d = minRowDistance (C)

  D = (columns (C) - C * C.') / 2;
  D(logical (eye (rows (C)))) = Inf;
  d = min (D(:));

endfunction

## Tests
%!test  # MATLAB parity: the pairwise design, three and four classes
%! assert_equal (designecoc (3, 'onevsone'), [1, 1, 0; -1, 0, 1; 0, -1, -1]);
%! assert_equal (designecoc (4, 'onevsone'), ...
%!               [1, 1, 1, 0, 0, 0; -1, 0, 0, 1, 1, 0; ...
%!                0, -1, 0, -1, 0, 1; 0, 0, -1, 0, -1, -1]);

%!test  # MATLAB parity: one class against the rest
%! assert_equal (designecoc (3, 'onevsall'), [1, -1, -1; -1, 1, -1; -1, -1, 1]);
%! assert_equal (designecoc (5, 'onevsall'), 2 * eye (5) - 1);

%!test  # MATLAB parity: every split in two, the first class always +1
%! assert_equal (designecoc (3, 'binarycomplete'), ...
%!               [1, 1, 1; 1, -1, -1; -1, 1, -1]);
%! assert_equal (designecoc (4, 'binarycomplete'), ...
%!               [1, 1, 1, 1, 1, 1, 1; 1, 1, 1, -1, -1, -1, -1; ...
%!                1, -1, -1, 1, 1, -1, -1; -1, 1, -1, 1, -1, 1, -1]);

%!test  # MATLAB parity: ordered classes, cut after each in turn
%! assert_equal (designecoc (3, 'ordinal'), [-1, -1; 1, -1; 1, 1]);
%! assert_equal (designecoc (4, 'ordinal'), ...
%!               [-1, -1, -1; 1, -1, -1; 1, 1, -1; 1, 1, 1]);

%!test  # MATLAB parity: every split with classes left out
%! assert_equal (designecoc (3, 'ternarycomplete'), ...
%!               [-1, -1, 0, 1, -1, -1; 1, -1, -1, -1, 0, 1; ...
%!                0, 1, 1, 1, 1, 1]);

%!test  # MATLAB parity: the widths of the five exact designs
%! for K = [4, 5, 6, 10]
%!   assert_equal (columns (designecoc (K, 'onevsone')), K * (K - 1) / 2);
%!   assert_equal (columns (designecoc (K, 'onevsall')), K);
%!   assert_equal (columns (designecoc (K, 'ordinal')), K - 1);
%!   assert_equal (columns (designecoc (K, 'binarycomplete')), 2 ^ (K - 1) - 1);
%! endfor

%!test  # MATLAB parity: the ternary design grows as its closed form says
%! assert_equal (columns (designecoc (6, 'ternarycomplete')), 301);

%!test  # MATLAB parity: two classes are one column whatever the design
%! for d = {'onevsone', 'onevsall', 'binarycomplete', 'ternarycomplete', ...
%!          'ordinal'}
%!   assert_equal (designecoc (2, d{1}), [-1; 1]);
%! endfor

%!test  # a dense random design leaves no class out
%! M = designecoc (10, 'denserandom', 'NumTrials', 20);
%! assert_equal (rows (M), 10);
%! assert_equal (any (M(:) == 0), false);

%!test  # a sparse random design does leave classes out
%! M = designecoc (10, 'sparserandom', 'NumTrials', 20);
%! assert_equal (rows (M), 10);
%! assert_equal (any (M(:) == 0), true);

%!test  # every column of a random design trains a real learner
%! M = designecoc (8, 'denserandom', 'NumTrials', 20);
%! assert_equal (all (any (M > 0, 1) & any (M < 0, 1)), true);

%!test  # a design is named without regard to case
%! assert_equal (designecoc (4, 'OneVsAll'), designecoc (4, 'onevsall'));

%!warning<use 'binarycomplete' instead of 'denserandom' for 5 or fewer classes.> ...
%! designecoc (4, 'denserandom', 'NumTrials', 5);

%!warning<use 'ternarycomplete' instead of 'sparserandom' for 4 or fewer classes.> ...
%! designecoc (4, 'sparserandom', 'NumTrials', 5);

## Test input validation
%!error<designecoc: too few input arguments.> designecoc (3)
%!error<designecoc: name-value arguments must be in pairs.> ...
%! designecoc (3, 'onevsone', 'NumTrials')
%!error<designecoc: K must be an integer of at least 2.> designecoc (1, 'onevsone')
%!error<designecoc: K must be an integer of at least 2.> designecoc (2.5, 'onevsone')
%!error<designecoc: K must be an integer of at least 2.> designecoc ('a', 'onevsone')
%!error<designecoc: NAME must be a character vector.> designecoc (3, 42)
%!error<designecoc: 'nosuch' is not a coding design.> designecoc (3, 'nosuch')
%!error<designecoc: invalid name-value argument.> ...
%! designecoc (3, 'onevsone', 'NoSuch', 1)
%!error<designecoc: 'NumTrials' must be a positive integer.> ...
%! designecoc (3, 'denserandom', 'NumTrials', 0)
