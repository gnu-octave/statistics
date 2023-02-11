## Copyright (C) 1996-2017 Kurt Hornik
## Copyright (C) 2023 Andreas Bertsatos <abertsatos@biol.uoa.gr>
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
## @deftypefn  {statistics} {[@var{h}, @var{pval}, @var{stats}] =} hotelling_t2test (@var{x})
## @deftypefnx {statistics} {[@dots{}] =} hotelling_t2test (@var{x}, @var{m})
## @deftypefnx {statistics} {[@dots{}] =} hotelling_t2test (@var{x}, @var{y})
## @deftypefnx {statistics} {[@dots{}] =} hotelling_t2test (@var{x}, @var{m}, @var{Name}, @var{Value})
## @deftypefnx {statistics} {[@dots{}] =} hotelling_t2test (@var{x}, @var{y}, @var{Name}, @var{Value})
##
## Compute Hotelling's T^2 ("T-squared") test for a single sample or two
## dependent samples (paired-samples).
##
## For a sample @var{x} from a multivariate normal distribution with unknown
## mean and covariance matrix, test the null hypothesis that
## @code{mean (@var{x}) == @var{m}}.
##
## For two dependent samples @var{x} and @var{y} from a multivariate normal
## distributions with unknown means and covariance matrices, test the null
## hypothesis that @code{mean (@var{x} - @var{y}) == 0}.
##
## @qcode{hotelling_t2test} treats NaNs as missing values, and ignores the
## corresponding rows.
##
## Name-Value pair arguments can be used to set statistical significance.
## @qcode{"alpha"} can be used to specify the significance level of the test
## (the default value is 0.05).
##
## If @var{h} is 1 the null hypothesis is rejected, meaning that the tested
## sample does not come from a multivariate distribution with mean @var{m}, or
## in case of two dependent samples that they do not come from the same
## multivariate distribution.  If @var{h} is 0, then the null hypothesis cannot
## be rejected and it can be assumed that it holds true.
##
## The p-value of the test is returned in @var{pval}.
##
## @var{stats} is a structure containing the value of the Hotelling's @math{T^2}
## test statistic in the field "Tsq", and the degrees of freedom of the F
## distribution in the fields "df1" and "df2".  Under the null hypothesis,
## @math{(n-p) T^2 / (p(n-1))} has an F distribution with @math{p} and
## @math{n-p} degrees of freedom, where @math{n} and @math{p} are the
## numbers of samples and variables, respectively.
##
## @seealso{hotelling_t2test2}
## @end deftypefn

function [h, pval, stats] = hotelling_t2test (x, my, varargin)

  ## Check for minimum number of input arguments
  if (nargin < 1)
    print_usage ();
  endif

  ## Check X being a valid data set
  if (isscalar (x) || ndims (x) > 2)
    error ("hotelling_t2test: X must be a vector or a 2D matrix.");
  endif

  ## Set default arguments
  alpha = 0.05;

  ## Fix MY when X is a single input argument
  if (nargin == 1)
    if (isvector (x))
      my = 0;
    elseif (ismatrix (x))
      [n, p] = size (x);
      my = zeros (1, p);
    endif
  endif

  ## When X and MY are of equal size, then assume paired-sample
  if (isequal (size (x), size(my)))
    x = x - my;
    if (isvector (x))
      my = 0;
    elseif (ismatrix (x))
      [n, p] = size (x);
      my = zeros (1, p);
    endif
  endif

  ## Remove rows containing any NaNs
  x = rmmissing (x);

  ## Check additional options
  i = 1;
  while (i <= length (varargin))
    switch lower (varargin{i})
      case "alpha"
        i = i + 1;
        alpha = varargin{i};
        ## Check for valid alpha
        if (! isscalar (alpha) || ! isnumeric (alpha) || ...
                    alpha <= 0 || alpha >= 1)
          error ("hotelling_t2test: invalid value for alpha.");
        endif
      otherwise
        error ("hotelling_t2test: invalid Name argument.");
    endswitch
    i = i + 1;
  endwhile

  ## Conditional error checking for X being a vector or matrix
  if (isvector (x))
    if (! isscalar (my))
      error ("hotelling_t2test: if X is a vector, M must be a scalar.");
    endif
    n = length (x);
    p = 1;
  elseif (ismatrix (x))
    [n, p] = size (x);
    if (n <= p)
      error ("hotelling_t2test: X must have more rows than columns.");
    endif
    if (isvector (my) && length (my) == p)
      my = reshape (my, 1, p);
    else
      error (strcat (["hotelling_t2test: if X is a matrix, M must be a"], ...
                     [" vector of length equal to the columns of X."]));
    endif
  endif

  ## Calculate the necessary statistics
  d = mean (x) - my;
  stats.Tsq = n * d * (cov (x) \ d');
  stats.df1 = p;
  stats.df2 = n - p;
  pval = 1 - fcdf ((n-p) * stats.Tsq / (p * (n-1)), stats.df1, stats.df2);

  ## Determine the test outcome
  ## MATLAB returns this a double instead of a logical array
  h = double (pval < alpha);

endfunction

## Test input validation
%!error<Invalid call to hotelling_t2test.  Correct usage> hotelling_t2test ();
%!error<hotelling_t2test: X must be a vector or a 2D matrix.> ...
%! hotelling_t2test (1);
%!error<hotelling_t2test: X must be a vector or a 2D matrix.> ...
%! hotelling_t2test (ones(2,2,2));
%!error<hotelling_t2test: invalid value for alpha.> ...
%! hotelling_t2test (ones(20,2), [0, 0], "alpha", 1);
%!error<hotelling_t2test: invalid value for alpha.> ...
%! hotelling_t2test (ones(20,2), [0, 0], "alpha", -0.2);
%!error<hotelling_t2test: invalid value for alpha.> ...
%! hotelling_t2test (ones(20,2), [0, 0], "alpha", "a");
%!error<hotelling_t2test: invalid value for alpha.> ...
%! hotelling_t2test (ones(20,2), [0, 0], "alpha", [0.01, 0.05]);
%!error<hotelling_t2test: invalid Name argument.> ...
%! hotelling_t2test (ones(20,2), [0, 0], "name", 0.01);
%!error<hotelling_t2test: if X is a vector, M must be a scalar.> ...
%! hotelling_t2test (ones(20,1), [0, 0]);
%!error<hotelling_t2test: X must have more rows than columns.> ...
%! hotelling_t2test (ones(4,5), [0, 0, 0, 0, 0]);
%!error<hotelling_t2test: if X is a matrix, M must be a vector of length> ...
%! hotelling_t2test (ones(20,5), [0, 0, 0, 0]);

## Test results
%!test
%! [h, pval, stats] = hotelling_t2test (randn(5000,5));
%! assert (h, 0);
%! assert (stats.df1, 5);
%! assert (stats.df2, 4995);
