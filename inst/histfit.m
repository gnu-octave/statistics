## Copyright (C) 2003 Alberto Terruzzi <t-albert@libero.it>
## Copyright (C) 2022-2024 Andreas Bertsatos <abertsatos@biol.uoa.gr>
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
## @deftypefn  {statistics} {} histfit (@var{x}, @var{nbins})
## @deftypefnx {statistics} {@var{h} =} histfit (@var{x}, @var{nbins})
##
## Plot histogram with superimposed fitted normal density.
##
## @code{histfit (@var{x}, @var{nbins})} plots a histogram of the values in
## the vector @var{x} using @var{nbins} bars in the histogram.  With one input
## argument, @var{nbins} is set to the square root of the number of elements in
## @var{x}.
##
## @code{@var{h} = histfit (@var{x}, @var{nbins})} returns the bins and fitted
## line handles of the plot in @var{h}.
##
## Example
##
## @example
## histfit (randn (100, 1))
## @end example
##
## @seealso{bar, hist, pareto}
## @end deftypefn

function [varargout] = histfit (varargin)

  ## Check for axes handle
  if (numel (varargin) < 1)
    error ("histfit: too few input arguments.");
  endif
  if (isaxes (varargin{1}))
    ax = varargin{1};
    varargin(1) = [];
  else
    ax = gca ();
  endif

  ## Get data
  if (numel (varargin) < 1)
    error ("histfit: too few input arguments.");
  else
    x = varargin{1};
    if (! isnumeric (x) || ! isreal (x) || ! isvector (x) || isscalar (x))
      error ("histfit: X must be a numeric vector of real numbers.");
    endif
    ## Remove missing values
    x(isnan (x)) = [];
    xsize = numel (x);
    ## Check for valid data
    if (xsize < 1)
      error ("histfit: no data in X.");
    endif
  endif

  ## Get nbins
  if (numel (varargin) > 1)
    nbins = varargin{2};
    if (! (isreal (nbins) && isscalar (nbins) && fix (nbins) == nbins))
      error ("histfit: NBINS must be a real scalar integer value.");
    endif
  else
    nbins = ceil (sqrt (xsize));
  endif

  ## Get distribution
  if (numel (varargin) > 2)
    distname = varargin{3};
    ## Add list of supported probability distribution objects
    PDO = {'Beta'; 'BirnbaumSaunders'; 'Burr'; 'Exponential'; ...
           'ExtremeValue'; 'Gamma'; 'GeneralizedExtremeValue'; ...
           'GeneralizedPareto'; 'InverseGaussian'; 'Kernel'; 'Logistic'; ...
           'Loglogistic'; 'Lognormal'; 'Nakagami'; 'NegativeBinomial'; ...
           'Normal'; 'Poisson'; 'Rayleigh'; 'Rician'; 'tLocationScale'; ...
           'Weibull'; "bisa"; "ev"; "gev"; "gp"; "invg"; "nbin"; "tls"; "wbl"};

    ## Check distribution name
    if (! (ischar (distname) && size (distname, 1) == 1))
      error ("histfit: DISTNAME must be a character vector.");
    elseif (! any (strcmpi (distname, PDO)))
      error ("histfit: unrecognized distribution name.");
    elseif (strcmpi (distname, "kernel"))
      error ("histfit: 'Kernel' distribution is not supported yet.");
    endif
  else
    distname = "normal";
  endif

  ## Plot the histogram
  if (strcmpi (distname, "poisson"))
    binwidth = 1;
    xmin = min (x) - 1;
    xmax = max (x) + 1;
    [binsize, bincenter] = hist (x, [xmin:xmax]);
  else
    [binsize, bincenter] = hist (x, nbins);
    binwidth = max (diff (bincenter));
    xmin = min (x) - binwidth / 2;
    xmax = max (x) + binwidth / 2;
  endif
  h = bar (ax, bincenter, binsize, 1, "facecolor", "b");

  ## Fit distibution to data
  pd = fitdist (x, distname);

  ## Compute density function
  if (strcmpi (distname, "poisson"))
    x = [min(x):max(x)]';
    y = pdf (pd, x);
  else
    x = [xmin:(xmax-xmin)/100:xmax]';
    y = pdf (pd, x);
  endif

  ## Normalize density line and overplot the histogram
  y = xsize * y * binwidth;
  hold on;
  if (strcmpi (distname, "poisson"))
    h(2) = plot (ax, x, y, ";;r-o");
  else
    h(2) = plot (ax, x, y, ";;r-");
  endif
  xlim ([xmin, xmax]);
  hold off;

  ## Return the plot's handle if requested
  if (nargout == 1)
    varargout{1} = h;
  endif
endfunction

%!demo
%! histfit (randn (100, 1))

%!demo
%! histfit (poissrnd (2, 1000, 1), 10, "Poisson")

%!demo
%! histfit (betarnd (3, 10, 1000, 1), 10, "beta")

## Test plotting
%!test
%! hf = figure ("visible", "off");
%! unwind_protect
%!   x = [2, 4, 3, 2, 4, 3, 2, 5, 6, 4, 7, 5, 9, 8, 10, 4, 11];
%!   histfit (x);
%! unwind_protect_cleanup
%!   close (hf);
%! end_unwind_protect
%!test
%! hf = figure ("visible", "off");
%! unwind_protect
%!   x = [2, 4, 3, 2, NaN, 3, 2, 5, 6, 4, 7, 5, 9, 8, 10, 4, 11];
%!   histfit (x);
%! unwind_protect_cleanup
%!   close (hf);
%! end_unwind_protect
%!test
%! hf = figure ("visible", "off");
%! unwind_protect
%!   x = [2, 4, 3, 2, NaN, 3, 2, 5, 6, 4, 7, 5, 9, 8, 10, 4, 11];
%!   histfit (x, 3);
%! unwind_protect_cleanup
%!   close (hf);
%! end_unwind_protect
%!test
%! hf = figure ("visible", "off");
%! unwind_protect
%!   histfit (randn (100, 1));
%! unwind_protect_cleanup
%!   close (hf);
%! end_unwind_protect
%!test
%! hf = figure ("visible", "off");
%! unwind_protect
%!   histfit (poissrnd (2, 1000, 1), 10, "Poisson");
%! unwind_protect_cleanup
%!   close (hf);
%! end_unwind_protect
%!test
%! hf = figure ("visible", "off");
%! unwind_protect
%!   histfit (betarnd (3, 10, 1000, 1), 10, "beta");
%! unwind_protect_cleanup
%!   close (hf);
%! end_unwind_protect
%!test
%! hf = figure ("visible", "off");
%! unwind_protect
%!   ax = gca ();
%!   histfit (ax, randn (100, 1));
%! unwind_protect_cleanup
%!   close (hf);
%! end_unwind_protect
%!test
%! hf = figure ("visible", "off");
%! unwind_protect
%!   ax = gca ();
%!   histfit (ax, poissrnd (2, 1000, 1), 10, "Poisson");
%! unwind_protect_cleanup
%!   close (hf);
%! end_unwind_protect
%!test
%! hf = figure ("visible", "off");
%! unwind_protect
%!   ax = gca ();
%!   histfit (ax, betarnd (3, 10, 1000, 1), 10, "beta");
%! unwind_protect_cleanup
%!   close (hf);
%! end_unwind_protect

## Test input validation
%!error<histfit: too few input arguments.> histfit ();
%!error<histfit: too few input arguments.> histfit (gca);
%!error<histfit: X must be a numeric vector of real numbers.> histfit ("wer");
%!error<histfit: no data in X.> histfit ([NaN, NaN, NaN]);
%!error<histfit: NBINS must be a real scalar integer value.> ...
%! histfit (randn (100, 1), 5.6);
%!error<histfit: DISTNAME must be a character vector.> ...
%! histfit (randn (100, 1), 8, 5);
%!error<histfit: DISTNAME must be a character vector.> ...
%! histfit (randn (100, 1), 8, {"normal"});
%!error<histfit: unrecognized distribution name.> ...
%! histfit (randn (100, 1), 8, "ASDASDASD");
%!error<histfit: 'Kernel' distribution is not supported yet.> ...
%! histfit (randn (100, 1), 8, "Kernel");
