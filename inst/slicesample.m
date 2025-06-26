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
## @deftypefn  {statistics} {[@var{smpl}, @var{neval}] =} slicesample (@var{start}, @var{nsamples}, @var{property}, @var{value}, @dots{})
##
## Draws @var{nsamples} samples from a target stationary distribution @var{pdf}
## using slice sampling of Radford M. Neal.
##
## Input:
## @itemize
## @item
## @var{start} is a 1 by @var{dim} vector of the starting point of the
## Markov chain. Each column corresponds to a different dimension.
##
## @item
## @var{nsamples} is the number of samples, the length of the Markov chain.
## @end itemize
##
## Next, several property-value pairs can or must be specified, they are:
##
## (Required properties) One of:
##
## @itemize
## @item
## @var{"pdf"}: the value is a function handle of the target stationary
## distribution to be sampled.  The function should accept different locations
## in each row and each column corresponds to a different dimension.
##
## or
##
## @item
## @var{logpdf}: the value is a function handle of the log of the target
## stationary distribution to be sampled. The function should accept different
## locations in each row and each column corresponds to a different dimension.
## @end itemize
##
## The following input property/pair values may be needed depending on the
## desired outut:
##
## @itemize
## @item
## "burnin" @var{burnin} the number of points to discard at the beginning, the default
## is 0.
##
## @item
## "thin" @var{thin} omitts @var{m}-1 of every @var{m} points in the generated
## Markov chain. The default is 1.
##
## @item
## "width" @var{width} the maximum Manhattan distance between two samples.
## The default is 10.
## @end itemize
##
## Outputs:
## @itemize
##
## @item
## @var{smpl} is a @var{nsamples} by @var{dim} matrix of random
## values drawn from @var{pdf} where the rows are different random values, the
## columns correspond to the dimensions of @var{pdf}.
##
## @item
## @var{neval} is the number of function evaluations per sample.
## @end itemize
## Example : Sampling from a normal distribution
##
## @example
## @group
## start = 1;
## nsamples = 1e3;
## pdf = @@(x) exp (-.5 * x .^ 2) / (pi ^ .5 * 2 ^ .5);
## [smpl, accept] = slicesample (start, nsamples, "pdf", pdf, "thin", 4);
## histfit (smpl);
## @end group
## @end example
##
## @seealso{rand, mhsample, randsample}
## @end deftypefn

function [smpl, neval] = slicesample (start, nsamples, varargin)

  if (nargin < 4)
    error ("slicesample: function called with too few input arguments.");
  endif

  sizestart = size (start);
  pdf = [];
  logpdf = [];
  width = 10;
  burnin = 0;
  thin = 1;
  for k = 1:2:length (varargin)
    if (ischar (varargin{k}))
      switch lower (varargin{k})
        case "pdf"
          if (isa (varargin{k+1}, "function_handle"))
            pdf = varargin{k+1};
          else
            error ("slicesample: pdf must be a function handle.");
          endif
        case "logpdf"
          if (isa (varargin{k+1}, "function_handle"))
            pdf = varargin{k+1};
          else
            error ("slicesample: logpdf must be a function handle.");
          endif
        case "width"
          if (numel (varargin{k+1}) == 1 || numel (varargin{k+1}) == sizestart(2))
            width = varargin{k+1}(:).';
          else
            error ("slicesample: width must be a scalar or 1 by dim vector.");
          endif
        case "burnin"
          if (varargin{k+1}>=0)
            burnin = varargin{k+1};
          else
            error ("slicesample: burnin must be greater than or equal to 0.");
          endif
        case "thin"
          if (varargin{k+1}>=1)
            thin = varargin{k+1};
          else
            error ("slicesample: thin must be greater than or equal to 1.");
          endif
        otherwise
          warning (["slicesample: Ignoring unknown option " varargin{k}]);
      endswitch
    else
      error (["slicesample: " varargin{k} " is not a valid property."]);
    endif
  endfor

  if (! isempty (pdf) && isempty (logpdf))
    logpdf = @(x) rloge (pdf (x));
  elseif (isempty (pdf) && isempty (logpdf))
    error ("slicesample: pdf or logpdf must be input.");
  endif
  dim = sizestart(2);
  smpl = zeros (nsamples, dim);

  if (all (sizestart == [1 dim]))
    smpl(1, :) = start;
  else
    error ("slicesample: start must be a 1 by dim vector.");
  endif

  maxit = 100;
  neval = 0;

  fgreaterthan = @(x, fxc) logpdf (x) >= fxc;

  ti = burnin + nsamples * thin;

  rndexp = rande (ti, 1);
  crand = rand (ti, dim);
  prand = rand (ti, dim);

  xc = smpl(1, :);
  for i = 1:ti
    neval++;
    sliceheight = logpdf (xc) - rndexp(i);
    c = width .* crand(i, :);
    lb = xc - c;
    ub = xc + width - c;
    #Only for single variable as bounds can not be found with point when dim > 1
    if (dim == 1)
      for k=1:maxit
        neval++;
        if (! fgreaterthan (lb, sliceheight))
          break
        endif
        lb -= width;
      end
      if (k == maxit)
        warning ("slicesample: Step out exceeded maximum iterations");
      endif
      for k = 1:maxit
        neval++;
        if (! fgreaterthan (ub, sliceheight))
          break
        endif
        ub += width;
      end
      if (k == maxit)
        warning ("slicesample: Step out exceeded maximum iterations");
      endif
    end
    xp = (ub - lb) .* prand(i, :) + lb;
    for k=1:maxit
      neval++;
      isgt = fgreaterthan (xp,sliceheight);
      if (all (isgt))
        break
      endif
      lc = ! isgt & xp < xc;
      uc = ! isgt & xp > xc;
      lb(lc) = xp(lc);
      ub(uc) = xp(uc);
      xp = (ub - lb) .* rand (1, dim) + lb;
    end
    if (k == maxit)
      warning ("slicesample: Step in exceeded maximum iterations");
    endif
    xc = xp;
    if (i > burnin)
      indx = (i - burnin) / thin;
      if rem (indx, 1) == 0
        smpl(indx, :) = xc;
      end
    end
  end
  neval = neval / (nsamples * thin + burnin);
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
%! Sigma += eye (d)*abs (eigs (Sigma, 1, "sa")) * 1.1;
%! pdf = @(x)(2*pi)^(-d/2)*det(Sigma)^-.5*exp(-.5*sum((x.'-mu).*(Sigma\(x.'-mu)),1));
%!
%! ## Inputs
%! start = ones (1,2);
%! nsamples = 500;
%! K = 500;
%! m = 10;
%! rande ("seed", 4);  rand ("seed", 5)  # for reproducibility
%! [smpl, accept] = slicesample (start, nsamples, "pdf", pdf, "burnin", K, "thin", m, "width", [20, 30]);
%! figure;
%! hold on;
%! plot (smpl(:,1), smpl(:,2), 'x');
%! [x, y] = meshgrid (linspace (-6,4), linspace(-3,7));
%! z = reshape (pdf ([x(:), y(:)]), size(x));
%! mesh (x, y, z, "facecolor", "None");
%!
%! ## Using sample points to find the volume of half a sphere with radius of .5
%! f = @(x) ((.25-(x(:,1)+1).^2-(x(:,2)-2).^2).^.5.*(((x(:,1)+1).^2+(x(:,2)-2).^2)<.25)).';
%! int = mean (f (smpl) ./ pdf (smpl));
%! errest = std (f (smpl) ./ pdf (smpl)) / nsamples^.5;
%! trueerr = abs (2/3*pi*.25^(3/2)-int);
%! fprintf ("Monte Carlo integral estimate int f(x) dx = %f\n", int);
%! fprintf ("Monte Carlo integral error estimate %f\n", errest);
%! fprintf ("The actual error %f\n", trueerr);
%! mesh (x,y,reshape (f([x(:), y(:)]), size(x)), "facecolor", "None");

%!demo
%! ## Integrate truncated normal distribution to find normilization constant
%! pdf = @(x) exp (-.5*x.^2)/(pi^.5*2^.5);
%! nsamples = 1e3;
%! rande ("seed", 4);  rand ("seed", 5)  # for reproducibility
%! [smpl, accept] = slicesample (1, nsamples, "pdf", pdf, "thin", 4);
%! f = @(x) exp (-.5 * x .^ 2) .* (x >= -2 & x <= 2);
%! x = linspace (-3, 3, 1000);
%! area (x, f(x));
%! xlabel ("x");
%! ylabel ("f(x)");
%! int = mean (f (smpl) ./ pdf (smpl));
%! errest = std (f (smpl) ./ pdf (smpl)) / nsamples ^ 0.5;
%! trueerr = abs (erf (2 ^ 0.5) * 2 ^ 0.5 * pi ^ 0.5 - int);
%! fprintf("Monte Carlo integral estimate int f(x) dx = %f\n", int);
%! fprintf("Monte Carlo integral error estimate %f\n", errest);
%! fprintf("The actual error %f\n", trueerr);

## Test output
%!test
%! start = 0.5;
%! nsamples = 1e3;
%! pdf = @(x) exp (-.5*(x-1).^2)/(2*pi)^.5;
%! [smpl, accept] = slicesample (start, nsamples, "pdf", pdf, "thin", 2, "burnin", 0, "width", 5);
%! assert (mean (smpl, 1), 1, .15);
%! assert (var (smpl, 1), 1, .25);

## Test input validation
%!error slicesample ();
%!error slicesample (1);
%!error slicesample (1, 1);

