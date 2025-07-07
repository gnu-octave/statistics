## Copyright (C) 2014 Tony Richardson <richardson.tony@gmail.com>
## Copyright (C) 2022 Andrew Penn <A.C.Penn@sussex.ac.uk>
## Copyright (C) 2022 Andreas Bertsatos <abertsatos@biol.uoa.gr>
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
## @deftypefn  {statistics} {[@var{h}, @var{pval}, @var{ci}, @var{stats}] =} ttest (@var{x})
## @deftypefnx {statistics} {[@var{h}, @var{pval}, @var{ci}, @var{stats}] =} ttest (@var{x}, @var{m})
## @deftypefnx {statistics} {[@var{h}, @var{pval}, @var{ci}, @var{stats}] =} ttest (@var{x}, @var{y})
## @deftypefnx {statistics} {[@var{h}, @var{pval}, @var{ci}, @var{stats}] =} ttest (@var{x}, @var{m}, @var{Name}, @var{Value})
## @deftypefnx {statistics} {[@var{h}, @var{pval}, @var{ci}, @var{stats}] =} ttest (@var{x}, @var{y}, @var{Name}, @var{Value})
##
## Test for mean of a normal sample with unknown variance.
##
## Perform a t-test of the null hypothesis @code{mean (@var{x}) ==
## @var{m}} for a sample @var{x} from a normal distribution with unknown
## mean and unknown standard deviation.  Under the null, the test statistic
## @var{t} has a Student's t distribution.  The default value of
## @var{m} is 0.
##
## If the second argument @var{y} is a vector, a paired-t test of the
## hypothesis @code{mean (@var{x}) = mean (@var{y})} is performed. If @var{x}
## and @var{y} are vectors, they must have the same size and dimensions.
##
## @var{x} (and @var{y}) can also be matrices. For matrices, @qcode{ttest}
## performs separate t-tests along each column, and returns a vector of results.
## @var{x} and @var{y} must have the same number of columns. The Type I error
## rate of the resulting vector of @var{pval} can be controlled by entering
## @var{pval} as input to the function @qcode{multcompare}.
##
## @qcode{ttest} treats NaNs as missing values, and ignores them.
##
## Name-Value pair arguments can be used to set various options.
## @qcode{"alpha"} can be used to specify the significance level
## of the test (the default value is 0.05).  @qcode{"tail"}, can be used
## to select the desired alternative hypotheses.  If the value is
## @qcode{"both"} (default) the null is tested against the two-sided
## alternative @code{mean (@var{x}) != @var{m}}.
## If it is @qcode{"right"} the one-sided alternative @code{mean (@var{x})
## > @var{m}} is considered.  Similarly for @qcode{"left"}, the one-sided
## alternative @code{mean (@var{x}) < @var{m}} is considered.
## When argument @var{x} is a matrix, @qcode{"dim"} can be used to select
## the dimension over which to perform the test.  (The default is the
## first non-singleton dimension).
##
## If @var{h} is 1 the null hypothesis is rejected, meaning that the tested
## sample does not come from a Student's t distribution.  If @var{h} is 0, then
## the null hypothesis cannot be rejected and it can be assumed that @var{x}
## follows a Student's t distribution.  The p-value of the test is returned in
## @var{pval}.  A 100(1-alpha)% confidence interval is returned in @var{ci}.
##
## @var{stats} is a structure containing the value of the test statistic
## (@var{tstat}), the degrees of freedom (@var{df}) and the sample's standard
## deviation (@var{sd}).
##
## @seealso{hotelling_t2test, ttest2, hotelling_t2test2}
## @end deftypefn

function [h, p, ci, stats] = ttest (x, my, varargin)

  ## Set default arguments
  my_default = 0;
  alpha = 0.05;
  tail = "both";

  ## Find the first non-singleton dimension of x
  dim = min (find (size (x) != 1));
  if (isempty (dim))
    dim = 1;
  endif

  if (nargin == 1)
    my = my_default;
  endif

  i = 1;
  while (i <= length (varargin))
    switch lower (varargin{i})
      case "alpha"
        i = i + 1;
        alpha = varargin{i};
      case "tail"
        i = i + 1;
        tail = varargin{i};
      case "dim"
        i = i + 1;
        dim = varargin{i};
      otherwise
        error ("ttest: Invalid Name argument.");
    endswitch
    i = i + 1;
  endwhile

  if (! isa (tail, "char"))
    error ("ttest: tail argument must be a string.");
  endif

  if (any (and (! isscalar (my), size (x) != size (my))))
    error ("ttest: Arrays in paired test must be the same size.");
  endif

  ## Set default values if arguments are present but empty
  if (isempty (my))
    my = my_default;
  endif

  ## This adjustment allows everything else to remain the
  ## same for both the one-sample t test and paired tests.
  x = x - my;
  if (! isscalar (my))
    my = 0;
  endif

  ## Calculate the test statistic value (tval)
  n = sum (!isnan (x), dim);
  x_bar = mean (x, dim, "omitnan");
  stats.tstat = [];
  stats.df = n - 1;
  stats.sd = std (x, 0, dim, "omitnan");
  x_bar_std = stats.sd ./ sqrt(n);
  tval = (x_bar) ./ x_bar_std;
  stats.tstat = tval;

  ## Based on the "tail" argument determine the P-value, the critical values,
  ## and the confidence interval.
  switch lower (tail)
    case "both"
      p = 2 * (1 - tcdf (abs (tval), n - 1));
      tcrit = - tinv (alpha / 2, n - 1);
      ci = [x_bar-tcrit.*x_bar_std; x_bar+tcrit.*x_bar_std] + my;
    case "left"
      p = tcdf (tval, n - 1);
      tcrit = - tinv (alpha, n - 1);
      ci = [-inf*ones(size(x_bar)); my+x_bar+tcrit.*x_bar_std];
    case "right"
      p = 1 - tcdf (tval, n - 1);
      tcrit = - tinv (alpha, n - 1);
      ci = [my+x_bar-tcrit.*x_bar_std; inf*ones(size(x_bar))];
    otherwise
      error ("ttest: Invalid value for tail argument.");
  endswitch

  ## Reshape the ci array to match MATLAB shaping
  if (isscalar (x_bar) && dim == 2)
    ci = ci(:)';
  elseif (size (x_bar, 2) < size (x_bar, 1))
    ci = reshape (ci(:), length (x_bar), 2);
  endif

  ## Determine the test outcome
  ## MATLAB returns this a double instead of a logical array
  h = double (p < alpha);
endfunction

%!test
%! x = 8:0.1:12;
%! [h, pval, ci] = ttest (x, 10);
%! assert (h, 0)
%! assert (pval, 1, 10*eps)
%! assert (ci, [9.6219 10.3781], 1E-5)
%! [h, pval, ci0] = ttest (x, 0);
%! assert (h, 1)
%! assert (pval, 0)
%! assert (ci0, ci, 2e-15)
%! [h, pval, ci] = ttest (x, 10, "tail", "right", "dim", 2, "alpha", 0.05);
%! assert (h, 0)
%! assert (pval, 0.5, 10*eps)
%! assert (ci, [9.68498 Inf], 1E-5)
%!error ttest ([8:0.1:12], 10, "tail", "invalid");
%!error ttest ([8:0.1:12], 10, "tail", 25);
