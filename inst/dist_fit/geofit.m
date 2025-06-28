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
## @deftypefn  {statistics} {@var{pshat} =} geofit (@var{x})
## @deftypefnx {statistics} {[@var{pshat}, @var{psci}] =} geofit (@var{x})
## @deftypefnx {statistics} {[@var{pshat}, @var{psci}] =} geofit (@var{x}, @var{alpha})
## @deftypefnx {statistics} {[@var{pshat}, @var{psci}] =} geofit (@var{x}, @var{alpha}, @var{freq})
##
## Estimate parameter and confidence intervals for the geometric distribution.
##
## @code{@var{pshat} = geofit (@var{x})} returns the maximum likelihood estimate
## (MLE) of the probability of success for the geometric distribution.  @var{x}
## must be a vector.
##
## @code{[@var{pshat}, @var{psci}] = geofit (@var{x}, @var{alpha})} also returns
## the @qcode{100 * (1 - @var{alpha})} percent confidence intervals of the
## estimated parameter.  By default, the optional argument @var{alpha} is 0.05
## corresponding to 95% confidence intervals.  Pass in @qcode{[]} for
## @var{alpha} to use the default values.
##
## @code{[@dots{}] = geofit (@var{x}, @var{alpha}, @var{freq})} accepts a
## frequency vector, @var{freq}, of the same size as @var{x}.  @var{freq}
## typically contains integer frequencies for the corresponding elements in
## @var{x}, but it can contain any non-integer non-negative values.  By default,
## or if left empty, @qcode{@var{freq} = ones (size (@var{x}))}.
##
## The geometric distribution models the number of failures (@var{x}) of a
## Bernoulli trial with probability @var{ps} before the first success.
##
## Further information about the geometric distribution can be found at
## @url{https://en.wikipedia.org/wiki/Geometric_distribution}
##
## @seealso{geocdf, geoinv, geopdf, geornd, geostat}
## @end deftypefn

function [pshat, psci] = geofit (x, alpha, freq)

  ## Check input arguments
  if (nargin < 1)
    error ("geofit: function called with too few input arguments.");
  endif

  ## Check data in X
  if (any (x < 0))
    error ("geofit: X cannot have negative values.");
  endif
  if (! isvector (x))
    error ("geofit: X must be a vector.");
  endif

  ## Check ALPHA
  if (nargin < 2 || isempty (alpha))
    alpha = 0.05;
  elseif (! isscalar (alpha) || ! isreal (alpha) || alpha <= 0 || alpha >= 1)
    error ("geofit: wrong value for ALPHA.");
  endif

  ## Check frequency vector
  if (nargin < 3 || isempty (freq))
    freq = ones (size (x));
  elseif (! isequal (size (x), size (freq)))
    error ("geofit: X and FREQ vector mismatch.");
  endif

  ## Expand frequency and censor vectors (if necessary)
  if (! all (freq == 1))
    xf = [];
    for i = 1:numel (freq)
      xf = [xf, repmat(x(i), 1, freq(i))];
    endfor
    x = xf;
    freq = ones (size (x));
  endif

  ## Compute PS estimate
  pshat = 1 ./ (1 + mean (x));

  ## Compute confidence interval of PS
  if (nargout > 1)
    sz = numel (x);
    serr = pshat .* sqrt ((1 - pshat) ./ sz);
    psci = norminv ([alpha/2; 1-alpha/2], [pshat; pshat], [serr; serr]);
  endif

endfunction

%!demo
%! ## Sample 2 populations from different geometric distributions
%! rande ("seed", 1);    # for reproducibility
%! r1 = geornd (0.15, 1000, 1);
%! rande ("seed", 2);    # for reproducibility
%! r2 = geornd (0.5, 1000, 1);
%! r = [r1, r2];
%!
%! ## Plot them normalized and fix their colors
%! hist (r, 0:0.5:20.5, 1);
%! h = findobj (gca, "Type", "patch");
%! set (h(1), "facecolor", "c");
%! set (h(2), "facecolor", "g");
%! hold on
%!
%! ## Estimate their probability of success
%! pshatA = geofit (r(:,1));
%! pshatB = geofit (r(:,2));
%!
%! ## Plot their estimated PDFs
%! x = [0:15];
%! y = geopdf (x, pshatA);
%! plot (x, y, "-pg");
%! y = geopdf (x, pshatB);
%! plot (x, y, "-sc");
%! xlim ([0, 15])
%! ylim ([0, 0.6])
%! legend ({"Normalized HIST of sample 1 with ps=0.15", ...
%!          "Normalized HIST of sample 2 with ps=0.50", ...
%!          sprintf("PDF for sample 1 with estimated ps=%0.2f", ...
%!                  mean (pshatA)), ...
%!          sprintf("PDF for sample 2 with estimated ps=%0.2f", ...
%!                  mean (pshatB))})
%! title ("Two population samples from different geometric distributions")
%! hold off

## Test output
%!test
%! x = 0:5;
%! [pshat, psci] = geofit (x);
%! assert (pshat, 0.2857, 1e-4);
%! assert (psci, [0.092499; 0.478929], 1e-5);
%!test
%! x = 0:5;
%! [pshat, psci] = geofit (x, [], [1 1 1 1 1 1]);
%! assert (pshat, 0.2857, 1e-4);
%! assert (psci, [0.092499; 0.478929], 1e-5);
%!assert (geofit ([1 1 2 3]), geofit ([1 2 3], [] ,[2 1 1]))

## Test input validation
%!error<geofit: function called with too few input arguments.> geofit ()
%!error<geofit: X cannot have negative values.> geofit (-1, [1 2 3 3])
%!error<geofit: wrong value for ALPHA.> geofit (1, 0)
%!error<geofit: wrong value for ALPHA.> geofit (1, 1.2)
%!error<geofit: wrong value for ALPHA.> geofit (1, [0.02 0.05])
%!error<geofit: X and FREQ vector mismatch.> ...
%! geofit ([1.5, 0.2], [], [0, 0, 0, 0, 0])
%!error<geofit: X and FREQ vector mismatch.> ...
%! geofit ([1.5, 0.2], [], [1, 1, 1])
