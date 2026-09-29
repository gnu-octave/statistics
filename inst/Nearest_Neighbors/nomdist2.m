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
## FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
## more details.
##
## You should have received a copy of the GNU General Public License along with
## this program; if not, see <http://www.gnu.org/licenses/>.

## -*- texinfo -*-
## @deftypefn  {statistics} {@var{D} =} nomdist2 (@var{X}, @var{Y})
## @deftypefnx {statistics} {@var{D} =} nomdist2 (@var{X}, @var{Y}, @var{measure})
## @deftypefnx {statistics} {@var{D} =} nomdist2 (@dots{}, @var{Name}, @var{Value})
##
## Dissimilarity between the rows of two samples of nominal data.
##
## @code{@var{D} = nomdist2 (@var{X}, @var{Y})} returns the Goodall 3
## dissimilarity of each row of @var{X} to each row of @var{Y}, as the
## @math{N*M} matrix @var{D}, where @var{X} has @math{N} rows and @var{Y}
## @math{M}@.  Their columns are the same nominal variables and their values
## the levels, of any type @code{nomdist} takes; where @var{X} is a table,
## @var{Y} must be one holding every variable of @var{X}, read by name.  A
## @var{Y} of no rows gives an @math{N*0} @var{D}.
##
## Level frequencies are counted over @var{X}, the sample the rows of @var{Y}
## are measured against, so @var{D} is not symmetric in the two:
## @code{nomdist2 (@var{X}, @var{Y})} is not in general
## @code{nomdist2 (@var{Y}, @var{X})'}.  @code{nomdist2 (@var{X}, @var{X})}
## equals @code{squareform (nomdist (@var{X}))} off the diagonal; on it,
## several measures give two identical rows a nonzero dissimilarity.
##
## @code{@var{D} = nomdist2 (@var{X}, @var{Y}, @var{measure})} uses the
## measure named, one of @qcode{'anderberg'}, @qcode{'burnaby'},
## @qcode{'eskin'}, @qcode{'gambaryan'}, @qcode{'goodall1'},
## @qcode{'goodall2'}, @qcode{'goodall3'}, @qcode{'goodall4'}, @qcode{'iof'},
## @qcode{'lin'}, @qcode{'lin1'}, @qcode{'of'}, @qcode{'sm'},
## @qcode{'smirnov'}, @qcode{'ve'} or @qcode{'vm'}, in any case, each defined
## as in @code{nomdist}.
##
## @code{@var{D} = nomdist2 (@dots{}, @var{Name}, @var{Value})} takes the
## options @code{nomdist} takes.  @qcode{'Weights'} gives one weight per
## variable.  @qcode{'Reference'} names a sample to count level frequencies
## over in place of @var{X}, which must hold every level @var{X} holds;
## @code{[@var{X}; @var{Y}]} counts them over both samples.
##
## @var{Y} may hold a level the reference sample does not.  Each measure then
## answers with the limit of its formula as that level's count falls to 0.
## @qcode{'sm'}, @qcode{'eskin'}, @qcode{'gambaryan'}, the four Goodall
## measures, @qcode{'smirnov'}, @qcode{'ve'} and @qcode{'vm'} never read that
## count and answer as usual.  @qcode{'of'} and @qcode{'burnaby'} score the
## variable 0, @qcode{'anderberg'} places the pair at 1, and @qcode{'lin'} and
## @qcode{'lin1'} at @code{Inf}.  Refused, where the limit does not exist:
## @qcode{'iof'}, whose similarity rises as a level grows rarer, and
## @qcode{'of'} and @qcode{'burnaby'} where the reference holds that variable
## at a single level.  A variable of weight 0 is left out, and so never
## refused.
##
## @seealso{nomdist, pdist2}
## @end deftypefn

function D = nomdist2 (X, Y, varargin)

  ## Input validation
  if (nargin < 2)
    error ("nomdist2: too few input arguments.");
  endif
  [C, measure, w, errmsg] = __nomprep__ (X, Y, true, varargin);
  if (! isempty (errmsg))
    error ("nomdist2: %s", errmsg);
  endif
  [Xc, Yc, Rc] = C{:};
  if (isempty (Yc))
    D = zeros (rows (Xc), 0);
    return;
  endif

  ## A level of Y the reference does not hold, where no limit exists
  if (any (strcmp (measure, {'iof', 'of', 'burnaby'})))
    if (isempty (Rc))
      Rc = Xc;
    endif
    for k = find (w > 0)
      if (all (ismember (Yc(:,k), Rc(:,k))))
        continue;
      endif
      if (strcmp (measure, 'iof'))
        error (strcat ("nomdist2: 'iof' cannot score a level of Y that", ...
                       " the reference sample does not hold."));
      endif
      if (all (Rc(:,k) == Rc(1,k)))
        error (strcat ("nomdist2: '%s' cannot score a level of Y that", ...
                       " the reference sample does not hold, where that", ...
                       " sample holds the variable at one level."), measure);
      endif
    endfor
  endif

  D = __nomdist__ (Xc, Yc, C{3}, measure, w);

endfunction

%!shared X
%! X = [1, 1, 1; 1, 2, 1; 1, 1, 2; 2, 2, 2; 2, 1, 1; 3, 2, 2; 4, 1, 3];
%!test
%! D = squareform (nomdist (X));
%! D2 = nomdist2 (X, X);
%! assert_equal (D2(! eye (7)), D(! eye (7)));
%!test
%! D = squareform (nomdist (X, 'lin1', 'Weights', [0.7, 1, 0.4]));
%! D2 = nomdist2 (X, X, 'lin1', 'Weights', [0.7, 1, 0.4]);
%! assert_equal (D2(! eye (7)), D(! eye (7)));
%!test
%! ## Identical rows need not be at 0: 1 - (3 - 24/42) / 3
%! D2 = nomdist2 (X, X);
%! assert_equal (D2(1,1), 4 / 21, -1e-14);
%!assert_equal (size (nomdist2 (X, X(1:2,:))), [7, 2])
%!assert_equal (nomdist2 (X, X(1:2,:), 'of'), nomdist2 (X, X, 'of')(:,1:2))
%!assert_equal (nomdist2 ([1; 1; 2], [1; 2; 2]), [1/3, 1, 1; 1/3, 1, 1; ...
%!                                                1, 0, 0], -1e-14)
%!assert_equal (nomdist2 ([1; 2; 2], [1; 1; 2]), [0, 0, 1; 1, 1, 1/3; ...
%!                                                1, 1, 1/3], -1e-14)
%!test
%! A = [1, 1; 2, 1; 2, 2];
%! B = [1, 2; 3, 1];
%! D = squareform (nomdist ([A; B], 'lin'));
%! assert_equal (nomdist2 (A, B, 'lin', 'Reference', [A; B]), D(1:3,4:5), ...
%!               -1e-14);
%!assert_equal (nomdist2 (X, zeros (0, 3)), zeros (7, 0))
%!test
%! T = table ([1; 2; 2], {'a'; 'a'; 'b'}, 'VariableNames', {'A', 'B'});
%! U = table ({'b'; 'a'}, [2; 1], 'VariableNames', {'B', 'A'});
%! assert_equal (nomdist2 (T, U, 'sm'), [1, 0; 0.5, 0.5; 0, 1]);
%!test
%! ## A new level of Y: OF scores its variable 0
%! s = 1 / (1 + log (3) * log (3 / 2));
%! D = nomdist2 ([1, 1; 2, 1; 2, 2], [3, 1; 3, 2], 'of');
%! assert_equal (D, [1, 2/s-1; 1, 2/s-1; 2/s-1, 1], -1e-14);
%!assert_equal (nomdist2 ({'a'; 'b'}, {'c'; 'a'}, 'burnaby'), [1, 0; 1, 0])
%!assert_equal (nomdist2 ([1; 2], [3; 1], 'anderberg'), [1, 0; 1, 1])
%!assert_equal (nomdist2 ([1; 2], [3; 1], 'lin'), [Inf, 0; Inf, Inf])
%!assert_equal (nomdist2 ([1; 2], [3; 1], 'goodall4'), [1, 1; 1, 1])
%!assert_equal (nomdist2 ([1, 1; 1, 1; 2, 1; 2, 1], [1, 3], 'iof', ...
%!                        'Weights', [1, 0]), [0; 0; 1; 1] * log (2) ^ 2, ...
%!              -1e-14)
%!assert_equal (nomdist2 ([1, 1; 1, 2], [2, 1], 'of', 'Weights', [0, 1]), ...
%!              [0; log(2) ^ 2], -1e-14)
%!assert_equal (nomdist2 ([1; 2], [3; 1], 'iof', 'Reference', ...
%!                        [1; 1; 2; 2; 3; 3]), [1, 0; 1, 1] * log (2) ^ 2, ...
%!              -1e-14)

%!error<nomdist2: too few input arguments.> nomdist2 ([1; 2])
%!error<nomdist2: Y must be a numeric, logical, categorical, string or cellstr matrix, or a table.> ...
%! nomdist2 ([1; 2], {1; 2})
%!error<nomdist2: Y must be a matrix, as X is.> ...
%! nomdist2 ([1; 2], table ([1; 2]))
%!error<nomdist2: Y must be a table, as X is.> ...
%! nomdist2 (table ([1; 2]), [1; 2])
%!error<nomdist2: Y must hold every variable of X.> ...
%! nomdist2 (table ([1; 2]), table ([1; 2], 'VariableNames', {'Z'}))
%!error<nomdist2: every variable of Y must be one column of levels.> ...
%! nomdist2 (table ([1; 2], 'VariableNames', {'A'}), ...
%!           table ([1, 2], 'VariableNames', {'A'}))
%!error<nomdist2: Y must have as many columns as X.> ...
%! nomdist2 ([1; 2], [1, 2])
%!error<nomdist2: Y must hold the same type as X in every variable.> ...
%! nomdist2 ([1; 2], {'a'})
%!error<nomdist2: Y must not hold missing values.> ...
%! nomdist2 ([1; 2], [1; NaN])
%!error<nomdist2: 'iof' cannot score a level of Y that the reference sample does not hold.> ...
%! nomdist2 ([1; 2], [3; 1], 'iof')
%!error<nomdist2: 'of' cannot score a level of Y that the reference sample does not hold, where that sample holds the variable at one level.> ...
%! nomdist2 ([1; 1], 2, 'of')
%!error<nomdist2: 'burnaby' cannot score a level of Y that the reference sample does not hold, where that sample holds the variable at one level.> ...
%! nomdist2 ([1; 1], 2, 'burnaby')
