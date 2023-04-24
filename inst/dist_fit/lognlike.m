## Copyright (C) 2023 Andreas Bertsatos <abertsatos@biol.uoa.gr>
##
## This file is part of the statistics package for GNU Octave.
##
## Octave is free software; you can redistribute it and/or modify it
## under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 3 of the License, or (at
## your option) any later version.
##
## Octave is distributed in the hope that it will be useful, but
## WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
## General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with Octave; see the file COPYING.  If not, see
## <http://www.gnu.org/licenses/>.

## -*- texinfo -*-
## @deftypefn  {statistics} {@var{nlogL} =} lognlike (@var{params}, @var{x})
## @deftypefnx {statistics} {[@var{nlogL}, @var{avar}] =} lognlike (@var{params}, @var{x})
## @deftypefnx {statistics} {[@dots{}] =} lognlike (@var{params}, @var{x}, @var{censor})
## @deftypefnx {statistics} {[@dots{}] =} lognlike (@var{params}, @var{x}, @var{censor}, @var{freq})
##
## Negative log-likelihood for the normal distribution.
##
## @code{@var{nlogL} = lognlike (@var{params}, @var{x})} returns the negative
## log-likelihood of the data in @var{x} for the log-normal distribution given
## parameters @qcode{@var{paramhat}([1, 2])}, which correspond to the mean and
## standard deviation of the associated normal distribution.  Missing values,
## @qcode{NaNs}, are ignored.  Negative values of @var{x} are treated as missing
## values.  @var{nlogL} is a scalar.
##
## @code{[@var{nlogL}, @var{avar}] = lognlike (@var{params}, @var{x})}
## returns the inverse of Fisher's information matrix, @var{avar}.  If the input
## parameter values in @var{params} are the maximum likelihood estimates, the
## diagonal elements of @var{avar} are their asymptotic variances.  @var{avar}
## is based on the observed Fisher's information, not the expected information.
##
## @code{[@dots{}] = lognlike (@var{params}, @var{x}, @var{censor})} accepts
## a boolean vector of the same size as @var{x} that is 1 for observations
## that are right-censored and 0 for observations that are observed exactly.
##
## @code{[@dots{}] = lognlike (@var{params}, @var{x}, @var{censor},
## @var{freq})} accepts a frequency vector of the same size as @var{x}.
## @var{freq} typically contains integer frequencies for the corresponding
## elements in @var{x}, but it may contain any non-integer non-negative
## values.  Pass in [] for @var{censor} to use its default value.
##
## @seealso{logncdf, logninv, lognpdf, lognrnd, lognfit, lognstat}
## @end deftypefn

function [nlogL, avar] = lognlike (params, x, censor, freq)

  ## Check input arguments
  if (nargin < 2)
    error ("lognlike: too few input arguments.");
  endif
  if (! isvector (x))
    error ("lognlike: X must be a vector.");
  endif
  if (numel (params) != 2)
    error ("lognlike: PARAMS must be a two-element vector.");
  endif
  if (nargin < 3 || isempty (censor))
    censor = [];
  elseif (! isequal (size (x), size (censor)))
    error ("lognlike: X and CENSOR vectors mismatch.");
  endif
  if nargin < 4 || isempty (freq)
    freq = [];
  elseif (isequal (size (x), size (freq)))
    nulls = find (freq == 0);
    if (numel (nulls) > 0)
      x(nulls) = [];
      if (numel (censor) == numel (freq))
        censor(nulls) = [];
      endif
      freq(nulls) = [];
    endif
  else
    error ("lognlike: X and FREQ vectors mismatch.");
  endif

  ## Treat negative data in X as missing values
  x(x < 0) = NaN;

  ## Calculate on log data
  logx = log(x);
  if (nargout <= 1)
    nlogL = normlike (params, logx, censor, freq);
  else
    [nlogL, avar] = normlike (params, logx, censor, freq);
  endif

  ## Compute censored and frequency
  if (isempty (freq))
    freq = 1;
  endif
  if (isempty (censor))
    censor = 0;
  endif
  nlogL = nlogL + sum (freq .* logx .* (1 - censor));

endfunction

## Test input validation
%!error<lognlike: too few input arguments.> lognlike ([12, 15]);
%!error<lognlike: X must be a vector.> lognlike ([12, 15], ones (2));
%!error<lognlike: PARAMS must be a two-element vector.> ...
%! lognlike ([12, 15, 3], [1:50]);
%!error<lognlike: X and CENSOR vectors mismatch.> ...
%! lognlike ([12, 15], [1:50], [1, 2, 3]);
%!error<lognlike: X and FREQ vectors mismatch.> ...
%! lognlike ([12, 15], [1:50], [], [1, 2, 3]);

## Results compared with Matlab
%!test
%! x = 1:50;
%! [nlogL, avar] = lognlike ([0, 0.25], x);
%! avar_out = [-5.4749e-03, 2.8308e-04; 2.8308e-04, -1.1916e-05];
%! assert (nlogL, 3962.330333301793, 1e-10);
%! assert (avar, avar_out, 1e-7);
%!test
%! x = 1:50;
%! [nlogL, avar] = lognlike ([0, 0.25], x * 0.5);
%! avar_out = [-7.6229e-03, 4.8722e-04; 4.8722e-04, -2.6754e-05];
%! assert (nlogL, 2473.183051225747, 1e-10);
%! assert (avar, avar_out, 1e-7);
%!test
%! x = 1:50;
%! [nlogL, avar] = lognlike ([0, 0.5], x);
%! avar_out = [-2.1152e-02, 2.2017e-03; 2.2017e-03, -1.8535e-04];
%! assert (nlogL, 1119.072424020455, 1e-12);
%! assert (avar, avar_out, 1e-6);
%!test
%! x = 1:50;
%! censor = ones (1, 50);
%! censor([2, 4, 6, 8, 12, 14]) = 0;
%! [nlogL, avar] = lognlike ([0, 0.5], x, censor);
%! avar_out = [-1.9823e-02, 2.0370e-03; 2.0370e-03, -1.6618e-04];
%! assert (nlogL, 1091.746371145497, 1e-12);
%! assert (avar, avar_out, 1e-6);
%!test
%! x = 1:50;
%! censor = ones (1, 50);
%! censor([2, 4, 6, 8, 12, 14]) = 0;
%! [nlogL, avar] = lognlike ([0, 1], x, censor);
%! avar_out = [-6.8634e-02, 1.3968e-02; 1.3968e-02, -2.1664e-03];
%! assert (nlogL, 349.3969104144271, 1e-12);
%! assert (avar, avar_out, 1e-6);
