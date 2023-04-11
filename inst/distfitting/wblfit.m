## Copyright (C) 2023 Andreas Bertsatos <abertsatos@biol.uoa.gr>
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
## @deftypefn  {statistics} {@var{paramhat} =} wblfit (@var{x})
## @deftypefnx {statistics} {[@var{paramhat}, @var{paramci}] =} wblfit (@var{x})
## @deftypefnx {statistics} {[@var{paramhat}, @var{paramci}] =} wblfit (@var{x}, @var{alpha})
## @deftypefnx {statistics} {[@dots{}] =} wblfit (@var{x}, @var{alpha}, @var{censor})
## @deftypefnx {statistics} {[@dots{}] =} wblfit (@var{x}, @var{alpha}, @var{censor}, @var{freq})
## @deftypefnx {statistics} {[@dots{}] =} wblfit (@var{x}, @var{alpha}, @var{censor}, @var{freq}, @var{options})
##
## Estimate mean and confidence intervals for the Weibull distribution.
##
## @code{@var{muhat} = wblfit (@var{x})} returns maximum likelihood estimates of
## the parameters of the Weibull distribution given the data in @var{x}.
## @qcode{@var{paramhat}(1)} is the scale parameter, @math{λ}, and
## @qcode{@var{paramhat}(2)} is the shape parameter, @math{k}.
##
## @code{[@var{paramhat}, @var{paramci}] = wblfit (@var{x})} returns the 95%
## confidence intervals for the parameter estimates.
##
## @code{[@dots{}] = wblfit (@var{x}, @var{alpha})} returns
## @qcode{100 * (1 - @var{alpha})} percent confidence intervals for the
## parameter estimates.  By default, @math{@var{alpha} = 0.05} corresponding
## fidence intervals.  Pass in @qcode{[]} for @var{alpha} to use the default
## values.
##
## @code{[@dots{}] = wblfit (@var{x}, @var{alpha}, @var{censor})}  accepts a
## boolean vector of the same size as @var{x} with 1 for observations that
## are right-censored and 0 for observations that are observed exactly.  By
## default, or if left empty, @qcode{@var{censor} = zeros (size (@var{x}))}.
##
## @code{[@dots{}] = wblfit (@var{x}, @var{alpha}, @var{censor}, @var{freq})}
## accepts a frequency vector of the same size as @var{x}.  @var{freq} typically
## contains integer frequencies for the corresponding elements in @var{x}, but
## may contain any non-integer non-negative values.  By default, or if left
## empty, @qcode{@var{freq} = ones (size (@var{x}))}.
##
## @code{[@dots{}] = evfit (@dots{}, @var{options})} specifies control
## parameters for the iterative algorithm used to compute ML estimates with the
## @code{fminsearch} function.  @var{options} is a structure with the following
## fields and their default values:
## @itemize
## @item @var{options}@qcode{.Display = "off"};
## @item @var{options}@qcode{.TolX = 1e-6};
## @end itemize
##
## @seealso{wblcdf, wblinv, wblpdf, wblrnd, wbllike, wblstat}
## @end deftypefn

function [parmhat, parmci] = wblfit (x, alpha, censor, freq, options)

  ## Check input arguments
  if (! isvector (x))
    error ("wblfit: X must be a vector.");
  elseif (any (x <= 0))
    error ("wblfit: X must contain only positive values.");
  endif

  ## Check alpha
  if (nargin < 2 || isempty (alpha))
    alpha = 0.05;
  else
    if (! isscalar (alpha) || ! isreal (alpha) || alpha <= 0 || alpha >= 1)
      error ("wblfit: Wrong value for ALPHA.");
    endif
  endif

  ## Check censor vector
  if (nargin < 3 || isempty (censor))
    censor = zeros (size (x));
  elseif (! isequal (size (x), size (censor)))
    error ("wblfit: X and CENSOR vector mismatch.");
  endif

  ## Check frequency vector
  if (nargin < 4 || isempty (freq))
    freq = ones (size (x));
  elseif (! isequal (size (x), size (freq)))
    error ("wblfit: X and FREQ vector mismatch.");
  endif

  ## Get options structure or add defaults
  if (nargin < 5)
    options.Display = "off";
    options.TolX = 1e-6;
  else
    if (! isstruct (options) || ! isfield (options, "Display") || ...
                                ! isfield (options, "TolX"))
      error (strcat (["wblfit: 'options' 5th argument must be a structure"], ...
                     [" with 'Display' and 'TolX' fields present."]));
    endif
  endif

  ## Fit an extreme value distribution to the logged data, then transform to
  ## the Weibull parameter scales.
  [parmhatEV, parmciEV] = evfit (log (x), alpha, censor, freq, options);
  parmhat = [exp(parmhatEV(1)), 1./parmhatEV(2)];
  if (nargout > 1)
    parmci = [exp(parmciEV(:,1)) 1./parmciEV([2 1],2)];
  endif

endfunction

%!demo
%! ## Sample 3 populations from 3 different Weibull distibutions
%! r = [wblrnd(2, 4, 2000, 1), wblrnd(5, 2, 2000, 1), wblrnd(1, 5, 2000, 1)];
%!
%! ## Plot them normalized and fix their colors
%! hist (r, 30, 2);
%! h = findobj(gca,'Type','patch');
%! set(h(1),'facecolor',"c");
%! set(h(2),'facecolor',"g");
%! set(h(3),'facecolor',"r");
%! ylim ([0, 2]);
%! xlim ([0, 15])
%! hold on
%!
%! ## Estimate their lambda parameter
%! lambda_kA = wblfit (r(:,1));
%! lambda_kB = wblfit (r(:,2));
%! lambda_kC = wblfit (r(:,3));
%!
%! ## Plot their estimated PDFs
%! x = [0:0.1:15];
%! y = wblpdf (x, lambda_kA(1), lambda_kA(2));
%! plot (x, y, "-pr");
%! y = wblpdf (x, lambda_kB(1), lambda_kB(2));
%! plot (x, y, "-sg");
%! y = wblpdf (x, lambda_kC(1), lambda_kC(2));
%! plot (x, y, "-^c");
%! hold off
%! legend ({"Normalized HIST of sample 1 with λ=2 and k=4", ...
%!          "Normalized HIST of sample 2 with λ=5 and k=2", ...
%!          "Normalized HIST of sample 3 with λ=1 and k=5", ...
%!          sprintf("PDF for sample 1 with estimated λ=%0.2f and k=%0.2f", ...
%!                  lambda_kA(1), lambda_kA(2)), ...
%!          sprintf("PDF for sample 2 with estimated λ=%0.2f and k=%0.2f", ...
%!                  lambda_kB(1), lambda_kB(2)), ...
%!          sprintf("PDF for sample 3 with estimated λ=%0.2f and k=%0.2f", ...
%!                  lambda_kC(1), lambda_kC(2))})
%! title ("Three population samples from different Weibull distibutions")
%! hold off

## test results
%!test
%! x = 1:50;
%! [paramhat, paramci] = wblfit (x);
%! paramhat_out = [28.3636, 1.7130];
%! paramci_out = [23.9531, 1.3551; 33.5861, 2.1655];
%! assert (paramhat, paramhat_out, 1e-4);
%! assert (paramci, paramci_out, 1e-4);
%!test
%! x = 1:50;
%! [paramhat, paramci] = wblfit (x, 0.01);
%! paramci_out = [22.7143, 1.2589; 35.4179, 2.3310];
%! assert (paramci, paramci_out, 1e-4);

## test input validation
%!error<wblfit: X must be a vector.> wblfit (ones (2,5));
%!error<wblfit: X must contain only positive values.> wblfit ([-1 2 3 4]);
%!error<wblfit: Wrong value for ALPHA.> wblfit ([1, 2, 3, 4, 5], 1.2);
%!error<wblfit: Wrong value for ALPHA.> wblfit ([1, 2, 3, 4, 5], 0);
%!error<wblfit: Wrong value for ALPHA.> wblfit ([1, 2, 3, 4, 5], "alpha");
%!error<wblfit: X and CENSOR vector mismatch.> ...
%! wblfit ([1, 2, 3, 4, 5], 0.05, [1 1 0]);
%!error<wblfit: X and CENSOR vector mismatch.> ...
%! wblfit ([1, 2, 3, 4, 5], [], [1 1 0 1 1]');
%!error<wblfit: X and FREQ vector mismatch.> ...
%! wblfit ([1, 2, 3, 4, 5], 0.05, zeros (1,5), [1 1 0]);
%!error<wblfit: X and FREQ vector mismatch.> ...
%! wblfit ([1, 2, 3, 4, 5], [], [], [1 1 0 1 1]');
%!error<wblfit: 'options' 5th argument must be a structure> ...
%! wblfit ([1, 2, 3, 4, 5], 0.05, [], [], 2);
