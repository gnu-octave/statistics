## Copyright (C) 2023 Andreas Bertsatos <abertsatos@biol.uoa.gr>
##
## This file is part of the statistics package for GNU Octave.
##
## This program is free software: you can redistribute it and/OR
## modify it under the terms of the GNU General Public License as
## published by the Free Software Foundation, either version 3 of the
## License, OR (at your option) any later version.
##
## This program is distributed in the hope that it will be useful, but
## WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
## General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program; see the file COPYING.  If not, see
## <http://www.gnu.org/licenses/>.

## -*- texinfo -*-
## @deftypefn  {statistics} {@var{h} =} fishertest (@var{x})
## @deftypefnx {statistics} {@var{h} =} fishertest (@var{x}, @var{param1}, @var{value1}, @dots{})
## @deftypefnx {statistics} {[@var{h}, @var{pval}] =} fishertest (@dots{})
## @deftypefnx {statistics} {[@var{h}, @var{pval}, @var{stats}] =} fishertest (@dots{})
##
## Fisher's exact test.
##
## @code{@var{h} = fishertest (@var{x})} performs Fisher's exact test on a
## @math{2x2} contingency table given in matrix @var{x}. This is a test of the
## hypothesis that there are no non-random associations between the two 2-level
## categorical variables in @var{x}.  @code{fishertest} returns the result of
## the tested hypothesis in @var{h}.  @var{h} = 0 indicates that the null
## hypothesis (of no association) cannot be rejected at the 5% significance
## level.  @var{h} = 1 indicates that the null hypothesis can be rejected at the
## 5% level.  @var{x} must contain only non-negative integers.  Use the
## @code{crosstab} function to generate the contingency table from samples of two
## categorical variables.  Fisher's exact test is not suitable when all integers
## in @var{x} are very large.  User can use the Chi-square test in this case.
##
## @code{[@var{h}, @var{pval}] = fishertest (@var{x})} returns the p-value in
## @var{pval}.  That is the probability of observing the given result, or one
## more extreme, by chance if the null hypothesis is true.  Small values of
## @var{pval} cast doubt on the validity of the null hypothesis.
##
## @code{[@var{p}, @var{pval}, @var{stats}] = fishertest (@dots{})} returns the
## structure @var{stats} with the following fields:
##
## @multitable @columnfractions 0.05 0.3 0.65
## @item @tab @qcode{OddsRatio} @tab -- the odds ratio
## @item @tab @qcode{ConfidenceInterval} @tab -- the asymptotic confidence
## interval for the odds ratio.  If any of the four entries in the contingency
## table @var{x} is zero, the confidence interval will not be computed, and
## @qcode{[-Inf Inf]} will be displayed.
## @end multitable
##
## @code{[@dots{}] = fishertest (@dots{}, @var{name}, @var{value}, @dots{})}
## specifies one or more of the following name/value pairs:
##
## @multitable @columnfractions 0.05 0.2 0.75
## @headitem @tab Name @tab Value
## @item @tab @qcode{"alpha"} @tab the significance level. Default is 0.05.
##
## @item @tab @qcode{"tail"} @tab a string specifying the alternative hypothesis
## @end multitable
## @multitable @columnfractions 0.1 0.25 0.65
## @item @tab @qcode{"both"} @tab odds ratio not equal to 1, indicating
## association between two variables (two-tailed test, default)
## @item @tab @qcode{"left"} @tab odds ratio greater than 1 (right-tailed test)
## @item @tab @qcode{"right"} @tab odds ratio is less than 1 (left-tailed test)
## @end multitable
##
## @seealso{crosstab, chi2test, mcnemar_test, ztest2}
## @end deftypefn

function [h, p, stats] = fishertest (x, varargin)

  if (nargin < 1)
    error ("fishertest: contingency table is missing.");
  endif

  if (nargin > 5)
    error ("fishertest: too many input parameters.");
  endif

  ## Check contingency table
  if (! ismatrix (x) || ndims (x) != 2)
    error ("fishertest: X must be a 2-dimensional matrix.");
  endif
  if (any (x(:) < 0) || any (isnan (x(:))) || any (isinf (x(:))) || ...
                        iscomplex (x) || any (fix (x(:)) != x(:)))
    error ("fishertest: X must contain only non-negative real integers.");
  endif
  if (all (x(:) >= 1e7))
    error ("fishertest: cannot handle large entries (>=1e7).");
  endif

  ## Add defaults and parse optional arguments
  alpha = 0.05;
  tail = "both";
  if (nargin > 1)
    params = numel (varargin);
    if ((params / 2) != fix (params / 2))
      error ("fishertest: optional arguments must be in Name-Value pairs.")
    endif
    for idx = 1:2:params
      name = varargin{idx};
      value = varargin{idx+1};
      switch (lower (name))
        case "alpha"
          alpha = value;
          if (! isscalar (alpha) || ! isnumeric (alpha) || ...
                alpha <= 0 || alpha >= 1)
            error ("fishertest: invalid value for alpha.");
          endif
        case "tail"
          tail = value;
          if (! any (strcmpi (tail, {"both", "left", "right"})))
            error ("fishertest: invalid value for tail.");
          endif
        otherwise
          error ("fishertest: invalid name for optional arguments.");
      endswitch
    endfor
  endif

  ## For 2x2 contingency table apply Fisher's exact test
  ## For larger tables apply the Fisher-Freeman-Halton variance

  if (all (size (x) == 2))

    ## Get margin sums
    r1 = sum (x(1,:));
    r2 = sum (x(2,:));
    c1 = sum (x(:,1));
    c2 = sum (x(:,2));
    sz = sum (x(:));

    ## Use try_catch block to avoid memory overflow for large numbers
    try
      if (strcmp (tail, "left"))
          p = hygecdf (x(1,1), sz, r1, c1);
      else
          if (min (r1, c1) <=  min (r2, c2))
              x11 = (0 : min (r1, c1))';
          else
              x22 = (0 : min (r2, c2))';
              x12 = c2 - x22;
              x11 = r1 - x12;
          endif
          switch tail
            case "both"
              p1 = hygepdf (x11, sz, r1, c1);
              p2 = hygepdf (x(1,1), sz, r1, c1);
              p = sum (p1(p1 < p2 + 10 * eps (p2)));
            case "right"
              xr = x11(x11 >= x(1,1));
              p = sum(hygepdf(xr,sz,r1,c1));
          endswitch
      endif
    catch
      error ("fishertest: cannot handle large entries.");
    end_try_catch

    ## Return test decision
    h = (p <= alpha);

    ## Calculate extra output arguments (if necessary)
    if (nargout > 2)
      OR = x(1,1) * x(2,2) / x(1,2) / x(2,1);
      if (any (x(:) == 0))
          CI = [-Inf, Inf];
      else
          SE = sqrt (1 / x(1,1) + 1 / x(1,2) + 1 / x(2,1) + 1 / x(2,2));
          LB = OR * exp (-norminv (1 - alpha / 2) * SE);
          UB = OR * exp (norminv (1 - alpha / 2) * SE);
          CI = [LB, UB];
      endif
      stats = struct ("OddsRatio", OR, "ConfidenceInterval", CI);
    endif

  else
    error ("fishertest: the Fisher-Freeman-Halton test is not implemented yet.");
  endif

endfunction

%!demo
%! ## A Fisher's exact test example
%!
%! x = [3, 1; 1, 3]
%! [h, p, stats] = fishertest(x)

## Test output against MATLAB R2018
%!assert (fishertest ([3, 4; 5, 7]), false);
%!assert (isa (fishertest ([3, 4; 5, 7]), "logical"), true);
%!test
%! [h, pval, stats] = fishertest ([3, 4; 5, 7]);
%! assert (pval, 1, 1e-14);
%! assert (stats.OddsRatio, 1.05);
%! CI = [0.159222057151289, 6.92429189601808];
%! assert (stats.ConfidenceInterval, CI, 1e-14)
%!test
%! [h, pval, stats] = fishertest ([3, 4; 5, 0]);
%! assert (pval, 0.08080808080808080, 1e-14);
%! assert (stats.OddsRatio, 0);
%! assert (stats.ConfidenceInterval, [-Inf, Inf])


## Test input validation
%!error fishertest ();
%!error fishertest (1, 2, 3, 4, 5, 6);
%!error<fishertest: X must be a 2-dimensional matrix.> ...
%! fishertest (ones (2, 2, 2));
%!error<fishertest: X must contain only non-negative real integers.> ...
%! fishertest ([1, 2; -3, 4]);
%!error<fishertest: X must contain only non-negative real integers.> ...
%! fishertest ([1, 2; 3, 4+i]);
%!error<fishertest: X must contain only non-negative real integers.> ...
%! fishertest ([1, 2; 3, 4.2]);
%!error<fishertest: X must contain only non-negative real integers.> ...
%! fishertest ([NaN, 2; 3, 4]);
%!error<fishertest: X must contain only non-negative real integers.> ...
%! fishertest ([1, Inf; 3, 4]);
%!error<fishertest: cannot handle large entries.> ...
%! fishertest (ones (2) * 1e8);
%!error<fishertest: invalid value for alpha.> ...
%! fishertest ([1, 2; 3, 4], "alpha", 0);
%!error<fishertest: invalid value for alpha.> ...
%! fishertest ([1, 2; 3, 4], "alpha", 1.2);
%!error<fishertest: invalid value for alpha.> ...
%! fishertest ([1, 2; 3, 4], "alpha", "val");
%!error<fishertest: invalid value for tail.>  ...
%! fishertest ([1, 2; 3, 4], "tail", "val");
%!error<fishertest: invalid value for tail.>  ...
%! fishertest ([1, 2; 3, 4], "alpha", 0.01, "tail", "val");
%!error<fishertest: invalid name for optional arguments.> ...
%! fishertest ([1, 2; 3, 4], "alpha", 0.01, "badoption", 3);
