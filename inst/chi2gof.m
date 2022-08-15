## Copyright (C) 2022 Andreas Bertsatos <abertsatos@biol.uoa.gr>
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
## @deftypefn {Function File} @var{h} = chi2gof (@var{x})
## @deftypefnx {Function File} [@var{h}, @var{p}] = chi2gof (@var{x})
## @deftypefnx {Function File} [@var{p}, @var{h}, @var{stats}] = chi2gof (@var{x})
## @deftypefnx {Function File} [@dots{}] = chi2gof (@var{x}, @var{name1}, @var{value1}, @dots{})
##
## Chi-square goodness-of-fit test.
##
## @code{chi2gof} performs a chi-square goodness-of-fit test for discrete or
## continuous distributions.  The test is performed by grouping the data into
## bins, calculating the observed and expected counts for those bins, and
## computing the chi-square test statistic
## @tex $$ \chi ^ 2 = \sum_{i=1}^N \left (O_i - E_i \right) ^ 2 / E_i $$
## @end tex
## @ifnottex
## SUM((O-E).^2./E),
## @end ifnottex
## where O is the observed counts and E is the expected counts.  This test
## statistic has an approximate chi-square distribution when the counts are
## sufficiently large.
##
## Bins in either tail with an expected count less than 5 are pooled with
## neighboring bins until the count in each extreme bin is at least 5.  If
## bins remain in the interior with counts less than 5, CHI2GOF displays a
## warning.  In that case, you should use fewer bins, or provide bin centers
## or binedges, to increase the expected counts in all bins.
##
## @code{@var{h} = chi2gof (@var{x})} performs a chi-square goodness-of-fit test
## that the data in the vector X are a random sample from a normal distribution
## with mean and variance estimated from @var{x}.  The result is @var{h} = 0 if
## the null hypothesis (that @var{x} is a random sample from a normal
## distribution) cannot be rejected at the 5% significance level, or @var{h} = 1
## if the nullhypothesis can be rejected at the 5% level.  @code{chi2gof} uses
## by default 10 bins ('nbins'), and compares the test statistic to a chi-square
## distribution with 'nbins' - 3 degrees of freedom, to take into account that
## two parameters were estimated.
##
## @code{[@var{h}, @var{p}] = chi2gof (@var{x})} also returns the p-value @var{p},
## which is the probability of observing the given result, or one more extreme,
## by chance if the null hypothesis is true.  If there are not enough degrees of
## freedom to carry out the test, @var{p} is NaN.
##
## @code{[@var{h}, @var{p}, @var{stats}] = chi2gof (@var{x})} also returns a
## @var{stats} structure with the following fields:
##
## @multitable @columnfractions 0.05 0.3 0.65
## @item @tab "chi2stat" @tab Chi-square statistic
## @item @tab "df" @tab Degrees of freedom
## @item @tab "binedges" @tab Vector of bin binedges after pooling
## @item @tab "O" @tab Observed count in each bin
## @item @tab "E" @tab Expected count in each bin
## @end multitable
##
## @code{[@dots{}] = chi2gof (@var{x}, @var{name1}, @var{value1}, @dots{})}
## specifies optional argument name/value pairs chosen from the following list.
##
## @multitable @columnfractions 0.05 0.2 0.75
## @headitem @tab Name @tab Value
## @item @tab "nbins" @tab The number of bins to use.  Default is 10.
## @item @tab "binctrs" @tab A vector of bin centers.
## @item @tab "binedges" @tab A vector of bin binedges.
## @item @tab "cdf" @tab A fully specified cumulative distribution function or a
## a function handle. Alternatively, a cell array whose first element is a
## function handle and all later elements are parameter values, one per cell.
## provided in a cell array whose first element is a function handle, and whose
## later elements are parameter values, one per cell.  The function must take X
## values as its first argument, and other parameters as later arguments.
## @item @tab "expected" @tab A vector with one element per bin specifying the
## expected counts for each bin.
## @item @tab "nparams" @tab The number of estimated parameters; used to adjust
## the degrees of freedom to be 'nbins' - 1 - 'nparams', where 'nbins' is the
## number of bins.
## @item @tab "emin" @tab The minimum allowed expected value for a bin; any bin
## in either tail having an expected value less than this amount is pooled with
## a neighboring bin.  Use the value 0 to prevent pooling.  Default is 5.
## @item @tab "frequency" @tab A vector of the same length as @var{x} containing
## the frequency of the corresponding @var{x} values.
## @item @tab "alpha" @tab An ALPHA value such that the hypothesis is rejected
## if @var{p} < ALPHA.  Default is ALPHA = 0.05.
## @end multitable
##
## You should specify either "cdf" or "expected" parameters, but not both.  If
## your "cdf" input contains extra parameters, these are accounted for
## automatically and there is no need to specify "nparams".  If your "expected"
## input depends on estimated parameters, you should use the "nparams" parameter
## to ensure that the degrees of freedom for the test is correct.
##
## @end deftypefn

function [h, p, stats] = chi2gof (x, varargin)

  ## Check imput arguments
  if (nargin < 1)
    error ("chi2gof: At least one imput argument is required.");
  endif
  if (! isvector(x) || ! isreal(x))
    error ("chi2gof: X must ba a vector of real numbers.");
  endif
  ## Add initial parameters
  nbins = [];
  binctrs = [];
  binedges = [];
  cdf_spec = [];
  expected = [];
  nparams = [];
  emin = 5;
  frequency = [];
  alpha = 0.05;
  ## Parse additional arguments
  numarg = nargin - 1;
  argpos = 1;
  while (numarg)
    argname = varargin{argpos};
    switch (lower (argname))
      case "nbins"
        nbins = varargin{argpos + 1};
      case "ctrs"
        binctrs = varargin{argpos + 1};
      case "edges"
        binedges = varargin{argpos + 1};
      case "cdf"
        cdf_spec = varargin{argpos + 1};
      case "expected"
        expected = varargin{argpos + 1};
      case "nparams"
        nparams = varargin{argpos + 1};
      case "emin"
        emin = varargin{argpos + 1};
      case "frequency"
        frequency = varargin{argpos + 1};
      case "alpha"
        alpha = varargin{argpos + 1};
    endswitch
    numarg -= 2;
    argpos += 2;
  endwhile
  ## Check additional arguments for errors
  if ((! isempty (nbins) + ! isempty (binctrs) + ! isempty (binedges)) > 1)
    error ("chi2gof: Inconsistent Arguments.");
  endif
  if ((! isempty (cdf_spec) + ! isempty (expected)) > 1)
    error ("chi2gof: Conflicted Arguments.");
  endif
  if (! isempty (frequency))
   if (! isvector (frequency) || numel (frequency) != numel (x))
       error ("chi2gof: X and Frequency vectors mismatch.");
   endif
   if (any (frequency < 0))
       error ("chi2gof: Frequency vector contains negative numbers.");
   endif
  endif
  if (! isscalar (emin) || emin < 0 || emin != round (emin) || ! isreal (emin))
    error("chi2gof: 'emin' must be a positive integer.");
  endif
  if (! isempty (nparams))
    if (! isscalar (nparams) || nparams < 0 || nparams != round (nparams) ...
                                           || ! isreal (nparams))
      error ("chi2gof: Wrong number of parameters.");
    endif
  endif
  if (! isscalar (alpha) || ! isreal (alpha) || alpha <= 0 || alpha >= 1)
     error ("chi2gof: Wrong value of alpha.");
  endif
  ## Make X a column vector
  x = x(:);
  ## Parse or create a frequeny vector
  if (isempty (frequency))
     frequency = ones (size (x));
  else
     frequency = frequency(:);
  endif
  ## Remove NaNs if any
  remove_NaNs = isnan (frequency) | isnan (x);
  if (any (remove_NaNs))
    x(remove_NaNs) = [];
    frequency(remove_NaNs) = [];
  endif
  ## Check for bin numbers, centers, or edges and calculate bins accordingly
  if (! isempty (binctrs))
    [Observed, binedges] = calculatebins (x, frequency, "ctrs", binctrs);
  elseif (! isempty (binedges))
    [Observed, binedges] = calculatebins (x, frequency, "edges", binedges);
  else
    if (isempty (nbins))
      if (isempty (expected))
        nbins = 10;                   ## default number of bins
      else
        nbins = length (expected);    ## determined by Expected vector
      endif
    endif
    [Observed, binedges] = calculatebins (x, frequency, "nbins", nbins);
  endif
  Observed = Observed(:);
  nbins = length (Observed);
  ## Calculate expected vector
  cdfargs = {};
  if (! isempty (expected))
    ## Provided as input argument
    if (! isvector (expected) || numel (expected) != nbins)
      error ("chi2gof: Expected counts vector is the wrong size.");
    endif
    if (any (expected < 0))
      error ("chi2gof: Expected counts vector has negative values.");
    endif
    Expected = expected(:);
  else
    ## Calculate from the cdf
    if (isempty (cdf_spec))
      ## Use estimated normal as default
      cdffunc = @normcdf;
      sumfreq = sum (frequency);
      mu = sum (x.*frequency)/sumfreq;
      sigma = sqrt (sum ((x.*frequency - mu) .^ 2) / (sumfreq-1));
      cdfargs = {mu, sigma};
      if (isempty (nparams))
        nparams = 2;
      endif
    elseif (isa (cdf_spec, "function_handle"))
      ## Split function handle to get function name and optional parameters
      cstr = ostrsplit (func2str (cdf_spec), ",");
      ## Simple function handle, no parameters: e.g. @normcdf
      if (isempty (strfind (cstr, "@")) && numel (cstr) == 1)
        cdffunc = str2func (char (strcat ("@", cstr)));
        if (isempty (nparams))
          nparams = numel (cdfargs);
        endif
      ## Complex function handle, no parameters: e.g. @(x) normcdf(x)
      elseif (! isempty (strfind (cstr, "@")) && numel (cstr) == 1)
        ## Remove white spaces
        cstr = char (cstr);
        cstr(strfind (cstr, " ")) = [];
        ## Remove input argument in parentheses
        while (length (strfind (cstr,"(")))
          cstr(index (cstr, "("):index (cstr, ")")) = [];
        endwhile
        cdffunc = str2func (cstr);
        if (isempty (nparams))
          nparams = numel (cdfargs);
        endif
      elseif (! isempty (strfind (cstr, "@")) && numel (cstr) > 1)
        ## Evaluate function name in first cell
        cstr_f = char (cstr(1));
        cstr_f(strfind (cstr_f, " ")) = [];
        cstr_f(index (cstr_f, "("):index (cstr_f, ")")) = [];
        cstr_f(index (cstr_f, "("):end) = [];
        cdffunc = str2func (cstr_f);
        ## Evaluate optional parameters in remaining cells
        cstr_idx = 2;
        while (cstr_idx <= numel (cstr))
          cstr_p = char (cstr(cstr_idx));
          cstr_p(strfind (cstr_p, " ")) = [];
          ## Check for numerical value
          if (isscalar (str2num (cstr_p)))
            cdfargs{cstr_idx - 1} = cstr_p;
          else
            ## Get function handle: e.g. mean
            cstr_p(index (cstr_p, "("):end) = [];
            cdfargs{cstr_idx - 1} = feval (str2func (cstr_p), x .* frequency);
            cstr_idx += 1;
          endif
        endwhile
        if (isempty (nparams))
          nparams = numel (cdfargs);
        endif
      endif
    elseif (iscell (cdf_spec))
      % Get function and args from cell array
      cdffunc = cdf_spec{1};
      cdfargs = cdf_spec(2:end);
      if (isempty (nparams))
        nparams = numel(cdfargs);
      endif
    endif
    if (! is_function_handle (cdffunc))
      error ("chi2gof: Poorly specified cumulative distribution function.");
    else
      cdfname = func2str (cdffunc);
    endif
    ## Calculate only inner bins, since tail probabilitiyis included in the
    ## calculation of expected counts for the first and last bins
    interioredges = binedges(2:end-1);
    ## Compute the cumulative probabilities
    Fcdf = feval (cdffunc, interioredges, cdfargs{:});
    if (! isvector(Fcdf) || numel (Fcdf) != (nbins - 1))
      msg = sprintf("chi2gof: Wrong number of outputs from: %s\n", cdfname);
      error (msg);
    endif
    % Compute the expected values
    Expected = sum(Observed) * diff([0;Fcdf(:);1]);
  endif
  ## Avoid too small expected values
  if (any (Expected < emin))
    [Expected, Observed, binedges] = poolbins (Expected, Observed, binedges, emin);
    nbins = length (Expected);
  end
  ## Compute test statistic
  cstat = sum(((Observed - Expected) .^ 2) ./ Expected);
  ## Calculate degrees of freedom
  if (isempty (nparams))
    nparams = 0;
  endif
  df = nbins - 1 - nparams;
  if (df > 0)
    p = 1 - chi2cdf (cstat, df);
  else
    df = 0;
    p = NaN;
  endif
  h = cast (p <= alpha, "double");
  ## Create 3rd output argument if necessary
  if (nargout > 2)
    stats.chi2stat = cstat;
    stats.df = df;
    stats.edges = binedges;
    stats.O = Observed';
    stats.E = Expected';
  endif
endfunction


function [Expected, Observed, binedges] = poolbins (Expected, ...
                                                    Observed, binedges, emin)
  i = 1;
  j = length(Expected);
  while (i < j - 1 && (Expected(i) < emin || Expected(i + 1) < emin || ...
                       Expected(j) < emin || Expected(j - 1) < emin))
    if (Expected(i) < Expected(j))
      Expected(i+1) = Expected(i+1) + Expected(i);
      Observed(i+1) = Observed(i+1) + Observed(i);
      i = i + 1;
    else
      Expected(j-1) = Expected(j-1) + Expected(j);
      Observed(j-1) = Observed(j-1) + Observed(j);
      j = j - 1;
    endif
  endwhile
  ## Keep only pooled bins
  Expected = Expected(i:j);
  Observed = Observed(i:j);
  binedges(j+1:end-1) = [];
  binedges(2:i) = [];
endfunction

function [Observed, binedges] = calculatebins (x, frequency, binspec, specval)
  lo = double (min (x(:)));
  hi = double (max (x(:)));
  ## Check binspec for bin count, bin centers, or bin edges.
  switch (binspec)
    case "nbins"
      nbins = specval;
      if (isempty (x))
        lo = 0;
        hi = 1;
      endif
      if (lo == hi)
          lo = lo - floor (nbins / 2) - 0.5;
          hi = hi + ceil (nbins / 2) - 0.5;
      endif
      binwidth = (hi - lo) ./ nbins;
      binedges = lo + binwidth * (0:nbins);
      binedges(length (binedges)) = hi;
    case "ctrs"
      binctrs = specval(:)';
      binwidth = diff (binctrs);
      binwidth = [binwidth binwidth(end)];
      binedges = [binctrs(1)-binwidth(1)/2 binctrs+binwidth/2];
    case "edges"
      binedges = specval(:)';
  endswitch
  ## Update bins
  nbins = length (binedges) - 1;
  ## Calculate bin numbers
  if (isempty (x))
    binnum = x;
  elseif (! isequal (binspec, "edges"))
    binedges = binedges + eps(binedges);
    [ignore, binnum] = histc (x, [-Inf binedges(2:end-1) Inf]);
  else
    [ignore, binnum] = histc (x, binedges);
    binnum(binnum == nbins + 1) = nbins;
  end
  ## Remove empty bins
  if (any (binnum == 0))
    frequency(binnum == 0) = [];
    binnum(binnum == 0) = [];
  end
  ## Compute Observed vector
  binnum = binnum(:);
  Observed = accumarray ([ones(size(binnum)), binnum], frequency, [1, nbins]);
endfunction

%!demo
%! x = normrnd (50, 5, 100, 1);
%! [h, p, stats] = chi2gof (x)
%! [h, p, stats] = chi2gof (x, "cdf", @(x)normcdf (x, mean(x), std(x)))
%! [h, p, stats] = chi2gof (x, "cdf", {@normcdf, mean(x), std(x)})
%!demo
%! x = rand (100,1 );
%! n = length (x);
%! binedges = linspace (0, 1, 11);
%! expectedCounts = n * diff (binedges);
%! [h, p, stats] = chi2gof (x, "binedges", binedges, "expected", expectedCounts)
%!demo
%! bins = 0:5;
%! obsCounts = [6 16 10 12 4 2];
%! n = sum(obsCounts);
%! lambdaHat = sum(bins.*obsCounts) / n;
%! expCounts = n * poisspdf(bins,lambdaHat);
%! [h, p, stats] = chi2gof (bins, "binctrs", bins, "frequency", obsCounts, ...
%!                          "expected", expCounts, "nparams",1)

## Test input validation
%!error chi2gof ()
%!error chi2gof ([2,3;3,4])
%!error chi2gof ([1,2,3,4], "nbins", 3, "ctrs", [2,3,4])
%!error chi2gof ([1,2,3,4], "frequency", [2,3,2])
%!error chi2gof ([1,2,3,4], "frequency", [2,3,2,-2])
%!error chi2gof ([1,2,3,4], "frequency", [2,3,2,2], "nparams", i)
%!error chi2gof ([1,2,3,4], "frequency", [2,3,2,2], "alpha", 1.3)
%!error chi2gof ([1,2,3,4], "expected", [-3,2,2])
%!error chi2gof ([1,2,3,4], "expected", [3,2,2], "nbins", 5)
%!error chi2gof ([1,2,3,4], "cdf", @normcdff)
%!test
%! x = [1 2 1 3 2 4 3 2 4 3 2 2];
%! [h, p, stats] = chi2gof (x);
%! assert (h, 0);
%! assert (p, NaN);
%! assert (stats.chi2stat, 0.1205375022748029, 1e-14);
%! assert (stats.df, 0);
%! assert (stats.edges, [1, 2.5, 4], 1e-14);
%! assert (stats.O, [7, 5], 1e-14);
%! assert (stats.E, [6.399995519909668, 5.600004480090332], 1e-14);
