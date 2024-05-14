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
## @deftypefn  {statistics} {} histfit (@var{x})
## @deftypefnx {statistics} {} histfit (@var{x}, @var{nbins})
## @deftypefnx {statistics} {} histfit (@var{x}, @var{nbins}, @var{distname})
## @deftypefnx {statistics} {} histfit (@var{ax}, @dots{})
## @deftypefnx {statistics} {@var{h} =} histfit (@dots{})
##
## Plot histogram with superimposed distribution fit.
##
## @code{histfit (@var{x)} plots a histogram of the values in the vector @var{x}
## using the number of bins equal to the square root of the number of nonmissing
## elements in @var{x} and superimposes a fitted normal density function.
##
## @code{histfit (@var{x}, @var{nbins})} plots a histogram of the values in the
## vector @var{x} using @var{nbins} number of bins in the histogram and
## superimposes a fitted normal density function.
##
## @code{histfit (@var{x}, @var{nbins}, @var{distname})} plots a histogram of
## the values in the vector @var{x} using @var{nbins} number of bins in the
## histogram and superimposes a fitted density function from the distribution
## specified by @var{distname}.
##
## @code{histfit (@var{ax}, @dots{})} uses the axes handle @var{ax} to plot the
## histogram and the fitted density function onto followed by any of the input
## argument combinations specified in the previous syntaxes.
##
## @code{@var{h} = histfit (@dots{})} returns a vector of handles @var{h}, where
## @qcode{@var{h}(1)} is the handle to the histogram and @qcode{@var{h}(1)} is
## the handle to the density curve.
##
## Note: calling @code{fitdist} without any input arguments will return a cell
## array of character vectors listing all supported distributions.
##
## @seealso{bar, hist, normplot, fitdist}
## @end deftypefn

function [varargout] = histfit (varargin)

  ## Add list of supported probability distribution objects
  PDO = {'Beta'; 'BirnbaumSaunders'; 'Burr'; 'Exponential'; 'ExtremeValue'; ...
         'Gamma'; 'GeneralizedExtremeValue'; 'GeneralizedPareto'; ...
         'InverseGaussian'; 'Logistic'; 'Loglogistic'; 'Lognormal'; ...
         'Nakagami'; 'NegativeBinomial'; 'Normal'; 'Poisson'; 'Rayleigh'; ...
         'Rician'; 'tLocationScale'; 'Weibull'};

  ABBR = {"bisa"; "ev"; "gev"; "gp"; "invg"; "nbin"; "tls"; "wbl"};

  ## Check for zero input arguments
  if (numel (varargin) < 1)
    varargout{1} = PDO;
    return
  endif

  ## Check for axes handle
  if (isaxes (varargin{1}))
    ax = varargin{1};
    varargin(1) = [];
    get_current_axes = false;
  else
    get_current_axes = true;
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
    ## Check distribution name
    if (! (ischar (distname) && size (distname, 1) == 1))
      error ("histfit: DISTNAME must be a character vector.");
    elseif (strcmpi (distname, "kernel"))
      error ("histfit: 'Kernel' distribution is not supported yet.");
    elseif (! (any (strcmpi (distname, PDO)) || any (strcmpi (distname, ABBR))))
      error ("histfit: unrecognized distribution name.");
    endif
  else
    distname = "normal";
  endif

  ## Create axes handle (if necessary)
  if (get_current_axes)
    ax = gca ();
  endif

  ## Plot the histogram
  if (any (strcmpi (distname, {"poisson", "NegativeBinomial", "nbin"})))
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
  if (any (strcmpi (distname, {"poisson", "NegativeBinomial", "nbin"})))
    x = [min(x):max(x)]';
    y = pdf (pd, x);
  else
    x = [xmin:(xmax-xmin)/100:xmax]';
    y = pdf (pd, x);
  endif

  ## Normalize density line and overplot the histogram
  y = xsize * y * binwidth;
  hold on;
  if (any (strcmpi (distname, {"poisson", "NegativeBinomial", "nbin"})))
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
%!test
%! hf = figure ("visible", "off");
%! unwind_protect
%!   ax = axes ("parent", hf);
%!   fail ("histfit (ax)", "histfit: too few input arguments.");
%! unwind_protect_cleanup
%!   close (hf);
%! end_unwind_protect
%!error<histfit: X must be a numeric vector of real numbers.> ...
%! histfit ('wer')
%!error<histfit: no data in X.> histfit ([NaN, NaN, NaN]);
%!error<histfit: NBINS must be a real scalar integer value.> ...
%! histfit (randn (100, 1), 5.6)
%!error<histfit: DISTNAME must be a character vector.> ...
%! histfit (randn (100, 1), 8, 5)
%!error<histfit: DISTNAME must be a character vector.> ...
%! histfit (randn (100, 1), 8, {'normal'})
%!error<histfit: 'Kernel' distribution is not supported yet.> ...
%! histfit (randn (100, 1), 8, 'Kernel')
%!error<histfit: unrecognized distribution name.> ...
%! histfit (randn (100, 1), 8, 'ASDASDASD')
