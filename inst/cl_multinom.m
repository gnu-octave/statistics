## Copyright (C) 2009 Levente Torok <TorokLev@gmail.com>
## Copyright (C) 2023 Andreas Bertsatos <abertsatos@biol.uoa.gr>
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
## @deftypefn  {statistics} {@var{CL} =} cl_multinom (@var{X}, @var{N}, @var{b})
## @deftypefnx {statistics} {@var{CL} =} cl_multinom (@var{X}, @var{N}, @var{b}, @var{method})
##
## Confidence level of multinomial portions.
##
## @code{cl_multinom} returns confidence level of multinomial parameters
## estimated as @math{p = X / sum(X)} with predefined confidence interval
## @var{b}.  Finite population is also considered.
##
## This function calculates the level of confidence at which the samples
## represent the true distribution given that there is a predefined tolerance
## (confidence interval). This is the upside down case of the typical excercises
## at which we want to get the confidence interval given the confidence level
## (and the estimated parameters of the underlying distribution).
## But once we accept (lets say at elections) that we have a standard predefined
## maximal acceptable error rate (e.g. @var{b}=0.02 ) in the estimation and we
## just want to know that how sure we can be that the measured proportions are
## the same as in the entire population (ie. the expected value and mean of the
## samples are roughly the same) we need to use this function.
##
## @subheading Arguments
## @multitable @columnfractions 0.1 0.01 0.10 0.01 0.78
## @headitem Variable @tab @tab Type @tab @tab Description
## @item @var{X} @tab @tab int vector @tab @tab sample frequencies bins.
## @item @var{N} @tab @tab int scalar @tab @tab Population size that was sampled
## by @var{X}.  If @qcode{N < sum (@var{X})}, infinite number assumed.
## @item @var{b} @tab @tab real vector @tab @tab confidence interval. If vector,
## it should be the size of @var{X} containing confence interval for each cells.
## If scalar, each cell will have the same value of b unless it is zero or -1.
## If value is 0, @var{b} = 0.02 is assumed which is standard choice at
## elections otherwise it is calculated in a way that one sample in a cell
## alteration defines the confidence interval.
## @item @var{method} @tab @tab string @tab @tab An optional argument
## for defining the calculation method.  Available choices are
## @qcode{"bromaghin"} (default), @qcode{"cochran"}, and @qcode{agresti_cull}.
## @end multitable
##
## Note!  The @qcode{agresti_cull} method is not exactly the solution at
## reference given below but an adjustment of the solutions above.
##
## @subheading Returns
## Confidence level.
##
## @subheading Example
## CL = cl_multinom ([27; 43; 19; 11], 10000, 0.05)
## returns 0.69 confidence level.
##
## @subheading References
## @enumerate
## @item
## "bromaghin" calculation type (default) is based on the article:
##
## Jeffrey F. Bromaghin, "Sample Size Determination for Interval Estimation
## of Multinomial Probabilities", The American Statistician  vol 47, 1993,
## pp 203-206.
##
## @item
## "cochran" calculation type is based on article:
##
## Robert T. Tortora, "A Note on Sample Size Estimation for Multinomial
## Populations", The American Statistician, , Vol 32. 1978,  pp 100-102.
##
## @item
## "agresti_cull" calculation type is based on article:
##
## A. Agresti and B.A. Coull, "Approximate is better than 'exact' for
## interval estimation of binomial portions", The American Statistician,
## Vol. 52, 1998, pp 119-126
## @end enumerate
##
## @end deftypefn

function CL = cl_multinom (X, N, b = 0.05, method = "bromaghin")

  if (nargin < 2 || nargin > 4)
    print_usage;
  elseif (! ischar (method))
    error ("cl_multinom: argument method must be a string.");
  endif

  k = rows (X);
  nn = sum (X);
  p = X / nn;

  if (isscalar (b))
    if (b==0)
      b=0.02;
    endif
    b = ones (rows (X), 1 ) * b;
    if (b<0)
      b = 1 ./ max (X, 1);
    endif
  endif
  bb = b .* b;

  if (N == nn)
    CL = 1;
    return;
  endif

  if (N < nn)
    fpc = 1;
  else
    fpc = (N - 1) / (N - nn); # finite population correction tag
  endif

  beta = p .* (1 - p);

  switch lower (method)
    case "cochran"
      t = sqrt (fpc * nn * bb ./ beta);
      alpha = (1 - normcdf (t)) * 2;

    case "bromaghin"
      t = sqrt (fpc * (nn * 2 * bb ) ./ ...
               (beta - 2 * bb + sqrt (beta .* beta - bb .* (4 * beta - 1))));
      alpha = (1 - normcdf (t)) * 2;

    case "agresti_cull"
      ts = fpc * nn * bb ./ beta ;
      if (k <= 2)
        alpha = 1 - chi2cdf (ts, k - 1); # adjusted Wilson interval
      else
        alpha = 1 - chi2cdf (ts / k, 1); # Goodman interval with Bonferroni arg.
      endif
    otherwise
      error ("cl_multinom: unknown calculation type '%s'.", method);
  endswitch

  CL = 1 - max( alpha );

endfunction

%!demo
%! CL = cl_multinom ([27; 43; 19; 11], 10000, 0.05)

## Test input validation
%!error<Invalid call to cl_multinom.  Correct usage> cl_multinom ();
%!error cl_multinom (1, 2, 3, 4, 5);
%!error<cl_multinom: argument method must be a string.> ...
%! cl_multinom (1, 2, 3, 4);
%!error<cl_multinom: unknown calculation type.> ...
%! cl_multinom (1, 2, 3, "some string");
