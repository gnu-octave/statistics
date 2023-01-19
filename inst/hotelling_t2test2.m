## Copyright (C) 1996-2017 Kurt Hornik
## Copyright (C) 2023 Andreas Bertsatos <abertsatos@biol.uoa.gr>
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
## @deftypefn  {statistics} [@var{h}, @var{pval}, @var{stats}] = hotelling_t2test2 (@var{x}, @var{y})
## @deftypefnx {statistics} [@dots{}] = hotelling_t2test2 (@var{x}, @var{y}, @var{Name}, @var{Value})
##
## Compute Hotelling's T^2 ("T-squared") test for two independent samples.
##
## For two samples @var{x} from multivariate normal distributions with
## the same number of variables (columns), unknown means and unknown
## equal covariance matrices, test the null hypothesis
## @code{mean (@var{x}) == mean (@var{y})}.
##
## @qcode{hotelling_t2test2} treats NaNs as missing values, and ignores the
## corresponding rows for each sample independently.
##
## Name-Value pair arguments can be used to set statistical significance.
## @qcode{"alpha"} can be used to specify the significance level of the test
## (the default value is 0.05).
##
## If @var{h} is 1 the null hypothesis is rejected, meaning that the tested
## samples do not come from the same multivariate distribution.  If @var{h} is
## 0, then the null hypothesis cannot be rejected and it can be assumed that
## both samples come from the same multivariate distribution.
##
## The p-value of the test is returned in @var{pval}.
##
## @var{stats} is a structure containing the value of the Hotelling's @math{T^2}
## test statistic in the field "Tsq", and the degrees of freedom of the F
## distribution in the fields "df1" and "df2".  Under the null hypothesis,
## @tex
## $$
## {(n_x+n_y-p-1) T^2 \over p(n_x+n_y-2)}
## $$
## @end tex
## @ifnottex
##
## @example
## (n_x+n_y-p-1) T^2 / (p(n_x+n_y-2))
## @end example
##
## @end ifnottex
## @noindent
## has an F distribution with @math{p} and @math{n_x+n_y-p-1} degrees of
## freedom, where @math{n_x} and @math{n_y} are the sample sizes and
## @math{p} is the number of variables.
##
## @seealso{hotelling_t2test}
## @end deftypefn

function [h, pval, stats] = hotelling_t2test2 (x, y, varargin)

  ## Check for minimum number of input arguments
  if (nargin < 2)
    print_usage ();
  endif
  ## Check X being a valid data set
  if (isscalar (x) || ndims (x) > 2)
    error ("hotelling_t2test2: X must be a vector or a 2D matrix.");
  endif
  ## Check Y being a valid data set
  if (isscalar (y) || ndims (y) > 2)
    error ("hotelling_t2test2: Y must be a vector or a 2D matrix.");
  endif
  
  ## Set default arguments
  alpha = 0.05;
  
  ## Remove rows containing any NaNs
  x = rmmissing (x);
  y = rmmissing (y);
  
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
          error ("hotelling_t2test2: invalid value for alpha.");
        endif
      otherwise
        error ("hotelling_t2test2: invalid Name argument.");
    endswitch
    i = i + 1;
  endwhile
  
  ## Conditional error checking for X being a vector or matrix
  if (isvector (x))
    n_x = length (x);
    if (! isvector (y))
      error ("hotelling_t2test2: if X is a vector, Y must also be a vector.");
    else
      n_y = length (y);
      p = 1;
    endif
  elseif (ismatrix (x))
    [n_x, p] = size (x);
    [n_y, q] = size (y);
    if (p != q)
      error (strcat (["hotelling_t2test2: X and Y must have the same"], ...
                     [" number of columns."]));
    endif
  endif

  ## Calculate the necessary statistics
  d = mean (x) - mean (y);
  S = ((n_x - 1) * cov (x) + (n_y - 1) * cov (y)) / (n_x + n_y - 2);
  stats.Tsq  = (n_x * n_y / (n_x + n_y)) * d * (S \ d');
  stats.df1 = p;
  stats.df2 = n_x + n_y - p - 1;
  pval = 1 - fcdf ((n_x + n_y - p - 1) * stats.Tsq / (p * (n_x + n_y - 2)), ...
                    stats.df1, stats.df2);

  ## Determine the test outcome
  ## MATLAB returns this a double instead of a logical array
  h = double (pval < alpha);

endfunction

## Test input validation
%!error<Invalid call to hotelling_t2test2.  Correct usage> hotelling_t2test2 ();
%!error<Invalid call to hotelling_t2test2.  Correct usage> ...
%! hotelling_t2test2 ([2, 3, 4, 5, 6]);
%!error<hotelling_t2test2: X must be a vector or a 2D matrix.> ...
%! hotelling_t2test2 (1, [2, 3, 4, 5, 6]);
%!error<hotelling_t2test2: X must be a vector or a 2D matrix.> ...
%! hotelling_t2test2 (ones (2,2,2), [2, 3, 4, 5, 6]);
%!error<hotelling_t2test2: Y must be a vector or a 2D matrix.> ...
%! hotelling_t2test2 ([2, 3, 4, 5, 6], 2);
%!error<hotelling_t2test2: Y must be a vector or a 2D matrix.> ...
%! hotelling_t2test2 ([2, 3, 4, 5, 6], ones (2,2,2));
%!error<hotelling_t2test2: invalid value for alpha.> ...
%! hotelling_t2test2 (ones (20,2), ones (20,2), "alpha", 1);
%!error<hotelling_t2test2: invalid value for alpha.> ...
%! hotelling_t2test2 (ones (20,2), ones (20,2), "alpha", -0.2);
%!error<hotelling_t2test2: invalid value for alpha.> ...
%! hotelling_t2test2 (ones (20,2), ones (20,2), "alpha", "a");
%!error<hotelling_t2test2: invalid value for alpha.> ...
%! hotelling_t2test2 (ones (20,2), ones (20,2), "alpha", [0.01, 0.05]);
%!error<hotelling_t2test2: invalid Name argument.> ...
%! hotelling_t2test2 (ones (20,2), ones (20,2), "name", 0.01);
%!error<hotelling_t2test2: if X is a vector, Y must also be a vector.> ...
%! hotelling_t2test2 (ones (20,1), ones (20,2));
%!error<hotelling_t2test2: X and Y must have the same number of columns.> ...
%! hotelling_t2test2 (ones (20,2), ones (25,3));

## Test results
%!test
%! [h, pval, stats] = hotelling_t2test2 (randn(5000,5), randn(3000,5));
%! assert (h, 0);
%! assert (stats.df1, 5);
%! assert (stats.df2, 7994);