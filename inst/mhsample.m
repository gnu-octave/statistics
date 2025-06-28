## Copyright (C) 1995-2022 The Octave Project Developers
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
## @deftypefn  {statistics} {[@var{smpl}, @var{accept}] =} mhsample (@var{start}, @var{nsamples}, @var{property}, @var{value}, @dots{})
##
## Draws @var{nsamples} samples from a target stationary distribution @var{pdf}
## using Metropolis-Hastings algorithm.
##
## Inputs:
##
## @itemize
## @item
## @var{start} is a @var{nchain} by @var{dim} matrix of starting points for each
## Markov chain.  Each row is the starting point of a different chain and each
## column corresponds to a different dimension.
##
## @item
## @var{nsamples} is the number of samples, the length of each Markov chain.
## @end itemize
##
## Some property-value pairs can or must be specified, they are:
##
## (Required) One of:
##
## @itemize
## @item
## "pdf" @var{pdf}: a function handle of the target stationary distribution to
## be sampled.  The function should accept different locations in each row and
## each column corresponds to a different dimension.
##
## or
##
## @item
## "logpdf" @var{logpdf}: a function handle of the log of the target stationary
## distribution to be sampled.  The function should accept different locations
## in each row and each column corresponds to a different dimension.
## @end itemize
##
## In case optional argument @var{symmetric} is set to false (the default), one
## of:
##
## @itemize
## @item
## "proppdf" @var{proppdf}: a function handle of the proposal distribution that
## is sampled from with @var{proprnd} to give the next point in the chain.  The
## function should accept two inputs, the random variable and the current
## location each input should accept different locations in each row and each
## column corresponds to a different dimension.
##
## or
##
## @item
## "logproppdf" @var{logproppdf}: the log of "proppdf".
## @end itemize
##
## The following input property/pair values may be needed depending on the
## desired output:
##
## @itemize
## @item
## "proprnd" @var{proprnd}: (Required) a function handle which generates random
## numbers from @var{proppdf}.  The function should accept different locations
## in each row and each column corresponds to a different dimension
## corresponding with the current location.
##
## @item
## "symmetric" @var{symmetric}: true or false based on whether @var{proppdf} is
## a symmetric distribution.  If true, @var{proppdf} (or @var{logproppdf}) need
## not be specified. The default is false.
##
## @item
## "burnin" @var{burnin} the number of points to discard at the beginning, the
## default is 0.
##
## @item
## "thin" @var{thin}: omits @var{thin}-1 of every @var{thin} points in the
## generated Markov chain.  The default is 1.
##
## @item
## "nchain" @var{nchain}: the number of Markov chains to generate.  The default
## is 1.
## @end itemize
##
## Outputs:
##
## @itemize
## @item
## @var{smpl}: a @var{nsamples} x @var{dim} x @var{nchain} tensor of random
## values drawn from @var{pdf}, where the rows are different random values, the
## columns correspond to the dimensions of @var{pdf}, and the third dimension
## corresponds to different Markov chains.
##
## @item
## @var{accept} is a vector of the acceptance rate for each chain.
## @end itemize
##
## Example : Sampling from a normal distribution
##
## @example
## @group
## start = 1;
## nsamples = 1e3;
## pdf = @@(x) exp (-.5 * x .^ 2) / (pi ^ .5 * 2 ^ .5);
## proppdf = @@(x,y) 1 / 6;
## proprnd = @@(x) 6 * (rand (size (x)) - .5) + x;
## [smpl, accept] = mhsample (start, nsamples, "pdf", pdf, "proppdf", ...
## proppdf, "proprnd", proprnd, "thin", 4);
## histfit (smpl);
## @end group
## @end example
##
## @seealso{rand, slicesample}
## @end deftypefn

function [smpl, accept] = mhsample (start, nsamples, varargin)

  if (nargin < 6)
    print_usage ();
  endif

  sizestart  = size (start);
  pdf        = [];
  proppdf    = [];
  logpdf     = [];
  logproppdf = [];
  proprnd    = [];
  sym        = false;
  K          = 0;                       # burnin
  m          = 1;                       # thin
  nchain     = 1;

  for k = 1:2:length (varargin)
    if (ischar (varargin{k}))
      switch lower(varargin{k})
        case "pdf"
          if (isa (varargin{k+1}, "function_handle"))
            pdf = varargin{k+1};
          else
            error ("mhsample: pdf must be a function handle");
          endif

        case "proppdf"
          if (isa (varargin{k+1}, "function_handle"))
            proppdf = varargin{k+1};
          else
            error ("mhsample: proppdf must be a function handle");
          endif

        case "logpdf"
          if (isa (varargin{k+1}, "function_handle"))
            pdf = varargin{k+1};
          else
            error ("mhsample: logpdf must be a function handle");
          endif

        case "logproppdf"
          if (isa (varargin{k+1}, "function_handle"))
            proppdf = varargin{k+1};
          else
            error ("mhsample: logproppdf must be a function handle");
          endif

        case "proprnd"
          if (isa (varargin{k+1}, "function_handle"))
            proprnd = varargin{k+1};
          else
            error ("mhsample: proprnd must be a function handle");
          endif

        case "symmetric"
          if (isa (varargin{k+1}, "logical"))
            sym = varargin{k+1};
          else
            error ("mhsample: sym must be true or false");
          endif

        case "burnin"
          if (varargin{k+1}>=0)
            K = varargin{k+1};
          else
            error ("mhsample: K must be greater than or equal to 0");
          endif

        case "thin"
          if (varargin{k+1} >= 1)
            m = varargin{k+1};
          else
            error ("mhsample: m must be greater than or equal to 1");
          endif

        case "nchain"
          if (varargin{k+1} >= 1)
            nchain = varargin{k+1};
          else
            error ("mhsample: nchain must be greater than or equal to 1");
          endif

        otherwise
          warning (["mhsample: Ignoring unknown option " varargin{k}]);
      endswitch
    else
      error (["mhsample: " varargin{k} " is not a valid property."]);
    endif
  endfor

  if (! isempty (pdf) && isempty (logpdf))
    logpdf=@(x) rloge (pdf (x));
  elseif (isempty (pdf) && isempty (logpdf))
    error ("mhsample: pdf or logpdf must be input.");
  endif
  if (! isempty (proppdf) && isempty (logproppdf))
    logproppdf = @(x, y) rloge (proppdf (x, y));
  elseif (isempty (proppdf) && isempty (logproppdf) && ! sym)
    error ("mhsample: proppdf or logproppdf must be input unless 'symmetrical' is true.");
  endif
  if (! isa (proprnd, "function_handle"))
    error ("mhsample: proprnd must be a function handle.");
  endif
  if (length (sizestart) == 2)
    sizestart = [sizestart 0];
  end
  smpl = zeros (nsamples, sizestart(2), nchain);

  if (all (sizestart([1 3]) == [1 nchain]))
    ## Could remove, not Matlab compatible but allows continuing chains
    smpl(1, :, :) = start;
  elseif (all (sizestart([1 3]) == [nchain 0]))
    smpl(1, :, :) = permute (start, [3, 2, 1]);
  elseif (all (sizestart([1 3]) == [1 0]))
    ## Could remove, not Matlab compatible but allows all chains to start
    ## at the same location
    smpl(1, :, :) = repmat (start,[1, 1, nchain]);
  else
    error ("mhsample: start must be a nchain by dim matrix.");
  endif
  cx = permute (smpl(1, :, :),[3, 2, 1]);
  accept = zeros (nchain, 1);
  i = 1;
  rnd = log (rand (nchain, nsamples*m+K));
  for k = 1:nsamples*m+K
    canacc = rem (k-K, m) == 0;
    px = proprnd (cx);
    if (sym)
      A = logpdf (px) - logpdf(cx);
    else
      A = (logpdf (px) + logproppdf (cx, px)) - (logpdf (cx) + logproppdf (px, cx));
    endif
    ac = rnd(:, k) < min (A, 0);
    cx(ac, :) = px(ac, :);
    accept(ac)++;
    if (canacc)
      smpl(i, :, :) = permute (cx, [3, 2, 1]);
    end
    if (k > K && canacc)
      i++;
    endif
  endfor
  accept ./= (nsamples * m + K);

endfunction


function y = rloge (x)

  y = -inf (size (x));
  xg0 = x > 0;
  y(xg0) = log (x(xg0));

endfunction

%!demo
%! ## Define function to sample
%! d = 2;
%! mu = [-1; 2];
%! rand ("seed", 5)  # for reproducibility
%! Sigma = rand (d);
%! Sigma = (Sigma + Sigma');
%! Sigma += eye (d) * abs (eigs (Sigma, 1, "sa")) * 1.1;
%! pdf = @(x)(2*pi)^(-d/2)*det(Sigma)^-.5*exp(-.5*sum((x.'-mu).*(Sigma\(x.'-mu)),1));
%! ## Inputs
%! start = ones (1, 2);
%! nsamples = 500;
%! sym = true;
%! K = 500;
%! m = 10;
%! rand ("seed", 8)  # for reproducibility
%! proprnd = @(x) (rand (size (x)) - .5) * 3 + x;
%! [smpl, accept] = mhsample (start, nsamples, "pdf", pdf, "proprnd", proprnd, ...
%!                            "symmetric", sym, "burnin", K, "thin", m);
%! figure;
%! hold on;
%! plot (smpl(:, 1), smpl(:, 2), 'x');
%! [x, y] = meshgrid (linspace (-6, 4), linspace(-3, 7));
%! z = reshape (pdf ([x(:), y(:)]), size(x));
%! mesh (x, y, z, "facecolor", "None");
%! ## Using sample points to find the volume of half a sphere with radius of .5
%! f = @(x) ((.25-(x(:,1)+1).^2-(x(:,2)-2).^2).^.5.*(((x(:,1)+1).^2+(x(:,2)-2).^2)<.25)).';
%! int = mean (f (smpl) ./ pdf (smpl));
%! errest = std (f (smpl) ./ pdf (smpl)) / nsamples ^ .5;
%! trueerr = abs (2 / 3 * pi * .25 ^ (3 / 2) - int);
%! printf ("Monte Carlo integral estimate int f(x) dx = %f\n", int);
%! printf ("Monte Carlo integral error estimate %f\n", errest);
%! printf ("The actual error %f\n", trueerr);
%! mesh (x, y, reshape (f([x(:), y(:)]), size(x)), "facecolor", "None");

%!demo
%! ## Integrate truncated normal distribution to find normalization constant
%! pdf = @(x) exp (-.5*x.^2)/(pi^.5*2^.5);
%! nsamples = 1e3;
%! rand ("seed", 5)  # for reproducibility
%! proprnd = @(x) (rand (size (x)) - .5) * 3 + x;
%! [smpl, accept] = mhsample (1, nsamples, "pdf", pdf, "proprnd", proprnd, ...
%!                            "symmetric", true, "thin", 4);
%! f = @(x) exp(-.5 * x .^ 2) .* (x >= -2 & x <= 2);
%! x = linspace (-3, 3, 1000);
%! area(x, f(x));
%! xlabel ('x');
%! ylabel ('f(x)');
%! int = mean (f (smpl) ./ pdf (smpl));
%! errest = std (f (smpl) ./ pdf (smpl)) / nsamples^ .5;
%! trueerr = abs (erf (2 ^ .5) * 2 ^ .5 * pi ^ .5 - int);
%! printf ("Monte Carlo integral estimate int f(x) dx = %f\n", int);
%! printf ("Monte Carlo integral error estimate %f\n", errest);
%! printf ("The actual error %f\n", trueerr);

## Test output
%!test
%! nchain = 1e4;
%! start = rand (nchain, 1);
%! nsamples = 1e3;
%! pdf = @(x) exp (-.5*(x-1).^2)/(2*pi)^.5;
%! proppdf = @(x, y) 1/3;
%! proprnd = @(x) 3 * (rand (size (x)) - .5) + x;
%! [smpl, accept] = mhsample (start, nsamples, "pdf", pdf, "proppdf", proppdf, ...
%!                            "proprnd", proprnd, "thin", 2, "nchain", nchain, ...
%!                            "burnin", 0);
%! assert (mean (mean (smpl, 1), 3), 1, .01);
%! assert (mean (var (smpl, 1), 3), 1, .01)

## Test input validation
%!error mhsample ();
%!error mhsample (1);
%!error mhsample (1, 1);
%!error mhsample (1, 1, "pdf", @(x)x);
%!error mhsample (1, 1, "pdf", @(x)x, "proprnd", @(x)x+rand(size(x)));

