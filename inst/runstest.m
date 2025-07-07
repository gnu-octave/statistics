## Copyright (C) 2013 Nir Krakauer <nkrakauer@ccny.cuny.edu>
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
## @deftypefn  {statistics} {@var{h} =} runstest (@var{x})
## @deftypefnx {statistics} {@var{h} =} runstest (@var{x}, @var{v})
## @deftypefnx {statistics} {@var{h} =} runstest (@var{x}, @qcode{"ud"})
## @deftypefnx {statistics} {@var{h} =} runstest (@dots{}, @var{Name}, @var{Value})
## @deftypefnx {statistics} {[@var{h}, @var{pval}, @var{stats}] =} runstest (@dots{})
##
## Run test for randomness in the vector @var{x}.
##
## @code{@var{h} = runstest (@var{x})} calculates the number of runs of
## consecutive values above or below the mean of @var{x} and tests the null
## hypothesis that the values in the data vector @var{x} come in random order.
## @var{h} is 1 if the test rejects the null hypothesis at the 5% significance
## level, or 0 otherwise.
##
## @code{@var{h} = runstest (@var{x}, @var{v})} tests the null hypothesis based
## on the number of runs of consecutive values above or below the specified
## reference value @var{v}.  Values exactly equal to @var{v} are omitted.
##
## @code{@var{h} = runstest (@var{x}, @qcode{"ud"})} calculates the number of
## runs up or down and tests the null hypothesis that the values in the data
## vector @var{x} follow a trend.  Too few runs indicate a trend, while too
## many runs indicate an oscillation.  Values exactly equal to the preceding
## value are omitted.
##
## @code{@var{h} = runstest (@dots{}, @var{Name}, @var{Value})} specifies
## additional options to the above tests by one or more @var{Name}-@var{Value}
## pair arguments.
##
## @multitable @columnfractions 0.15 0.05 0.8
## @headitem Name @tab @tab Value
## @item @qcode{"alpha"} @tab @tab the significance level. Default is 0.05.
##
## @item @qcode{"method"} @tab @tab a string specifying the method used to
## compute the p-value of the test.  It can be either @qcode{"exact"} to use an
## exact algorithm, or @qcode{"approximate"} to use a normal approximation.  The
## default is @qcode{"exact"} for runs above/below, and for runs up/down when
## the length of x is less than or equal to 50.  When testing for runs up/down
## and the length of @var{x} is greater than 50, then the default is
## @qcode{"approximate"}, and the @qcode{"exact"} method is not available.
##
## @item @qcode{"tail"} @tab @tab a string specifying the alternative hypothesis
## @end multitable
## @multitable @columnfractions 0.2 0.15 0.05 0.5
## @item @tab @qcode{"both"} @tab @tab two-tailed (default)
## @item @tab @qcode{"left"} @tab @tab left-tailed
## @item @tab @qcode{"right"} @tab @tab right-tailed
## @end multitable
##
## @seealso{signrank, signtest}
## @end deftypefn

function [h, pval, stats] = runstest (x, v, varargin)

  ## Check arguments
  if (nargin < 1)
    print_usage;
  endif

  ## Check X being a vector of scalar values
  if (! isvector (x) || ! isnumeric (x))
    error ("runstest: X must be a vector a scalar values.");
  else
    ## Remove missing values (NaNs)
    x(isnan (x)) = [];
  endif

  ## Check second argument being either a scalar reference number or "ud" string
  if (nargin > 1)
    if (isempty (v))
      v = mean (x);
    endif
    if (isnumeric (v) && isscalar (v))
      x = sign (x - v);
      rm = x == 0;
      if (sum (rm) > 0)
        warning ("runstest: %d elements equal to V were omitted.", sum (rm));
      endif
      x(rm) = [];
      N = numel(x);
      UD = false;
    elseif (strcmpi (v, "ud"))
      x = diff (x);
      rm = x == 0;
      if (sum (rm) > 0)
        warning ("runstest: %d repeated elements were omitted.", sum (rm));
      endif
      x(rm) = [];
      N = numel(x) + 1;
      UD = true;
    else
      error ("runstest: V must be either a scalar number or 'ud' char string.");
    endif
    v = v;
  else
    v = mean (x);
    x = sign (x - v);
    rm = x == 0;
    if (sum (rm) > 0)
      warning ("runstest: %d elements equal to 'mean(X)' were omitted.", ...
               sum (rm));
    endif
    x(rm) = [];
    N = numel(x);
    UD = false;
  endif

  ## Get number of runs
  n_up = sum (x==1);
  n_dn = numel (x) - n_up;

  ## Add defaults
  alpha = 0.05;
  if (N < 50 || ! UD)
    method = "exact";
  else
    method = "approximate";
  endif
  tail = "both";

  ## Parse optional arguments and validate parameters
  while (numel (varargin) > 1)
    switch (lower (varargin{1}))
      case "alpha"
        alpha = varargin{2};
        if (! isscalar (alpha) ||
            ! isnumeric (alpha) || alpha <= 0 || alpha >= 1)
          error ("runstest: invalid value for alpha.");
        endif

      case "method"
        method = varargin{2};
        if (! any (strcmpi (method, {"exact", "approximate"})))
          error ("runstest: invalid value for method.");
        endif
        if (strcmpi (method, "exact") && N > 50)
          warning ("runstest: exact method is not available for N > 50.");
          method = "approximate";
        endif

      case "tail"
        tail = varargin{2};
        if (! any (strcmpi (tail, {"both", "left", "right"})))
          error ("runstest: invalid value for tail.");
        endif

      otherwise
        error ("runstest: invalid optional argument.");
    endswitch
    varargin([1:2]) = [];
  endwhile

  ## Do the calculations here
  if (N > 0)
    R_num = sum (x([1:end-1]) != x([2:end])) + 1;
    ##R_num = sum ((x(1:(end-1)) .* x(2:end)) < 0) + 1;   #number of runs

    ## Special case
    if (N == 1)
      z = NaN;
    ## Compute with z statistic
    else
      ## Handle up/down or above/below
      if (UD)
          R_bar = (2 * N - 1) / 3;
          R_std = sqrt ((16 * N - 29) / 90);
      else
          R_bar = 1 + 2 * n_up * n_dn / N;
          R_std = sqrt (2 * n_up * n_dn * (2 * n_up * n_dn - N) / ...
                       (N ^ 2 * (N - 1)));
      end
      ## Handle tail
      if (strcmpi (tail, "both"))
          tc = -0.5 * sign (R_num - R_bar);
      elseif (strcmpi (tail, "left"))
          tc = 0.5;
      else
          tc = -0.5;
      endif
      ## Compute z value
      if (R_std > 0)
        z = (R_num + tc - R_bar) / R_std;
      else
        z = Inf * sign (R_num + tc - R_bar);
      endif
    endif
    ## Exact method
    if (strcmpi (method, "exact"))
      if (UD)
        R_max = N - 1;
        ## Get precalculated results from rundist.mat file
        temp = load ("rundist.mat");
        runD = temp.rundist;
        M = runD{N};
        p = M / sum (M);
        p = p([1:R_max]);
      else
        R_max = 2 * min ([n_up, n_dn]) + 1;
        if (n_up == 0 || n_dn == 0)
          p = 1;
        else
          R_vec = [1:R_max];
          p = zeros (size (R_vec));
          t = mod (R_vec, 2) == 0;
          ## Compute even
          if (any (t))
            k = R_vec(t) / 2;
            p(t) = 2 * exp (logBinoCoeff (n_up - 1, k - 1) + ...
                            logBinoCoeff (n_dn - 1, k - 1) - ...
                            logBinoCoeff (N, n_dn));
          endif
          ## Compute odd
          if (any (! t))
            k = floor (R_vec(! t) / 2);
            logdenom = logBinoCoeff (N, n_dn);
            p(! t) = exp (logBinoCoeff (n_up - 1, k - 1) + ...
                          logBinoCoeff (n_dn - 1, k) - logdenom) + ...
                     exp (logBinoCoeff (n_up - 1, k) + ...
                          logBinoCoeff (n_dn - 1, k - 1) - logdenom);
          endif
        endif
      endif

      if (isempty (p))
        p_ex = 1;
      else
        p_ex = p(R_num);
      end
      p_lo = sum (p([1:R_num-1]));
      p_hi = sum (p([R_num+1:end]));

    else
      ## Compute with z statistic
      p_ex = 0;
      p_lo = normcdf (z);
      p_hi = normcdf (-z);
    end
  ## Assume a constant vector in data
  else
    R_num = NaN;
    p_ex = 1;
    p_lo = 0;
    p_hi = 0;
    z = NaN;
  endif

  ## Compute tail probability
  if (strcmpi (tail, "both"))
    pval = min([1, 2*(p_ex + min ([p_lo, p_hi]))]);
  elseif (strcmpi (tail, "left"))
    pval = p_ex + p_lo;
  else
    pval = p_ex + p_hi;
  endif

  ## Return decision of test
  h = double (pval <= alpha);

  if (nargout > 2)
    stats.nruns = R_num;
    stats.n1 = n_up;
    stats.n0 = n_dn;
    stats.z = z;
  endif

endfunction

## Compute the log of the binomial coefficient
function logBC = logBinoCoeff(N,n)
  logBC = gammaln (N + 1) - gammaln (n + 1) - gammaln (N - n + 1);
endfunction

%!test
%! ## NIST beam deflection data
%! ## http://www.itl.nist.gov/div898/handbook/eda/section4/eda425.htm
%! data = [-213, -564, -35, -15, 141, 115, -420, -360, 203, -338, -431, ...
%!          194, -220, -513, 154, -125, -559, 92, -21, -579, -52, 99, -543, ...
%!         -175, 162, -457, -346, 204, -300, -474, 164, -107, -572, -8, 83, ...
%!         -541, -224, 180, -420, -374, 201, -236, -531, 83, 27, -564, -112, ...
%!          131, -507, -254, 199, -311, -495, 143, -46, -579, -90, 136, ...
%!         -472, -338, 202, -287, -477, 169, -124, -568, 17, 48, -568, -135, ...
%!          162, -430, -422, 172, -74, -577, -13, 92, -534, -243, 194, -355, ...
%!         -465, 156, -81, -578, -64, 139, -449, -384, 193, -198, -538, 110, ...
%!          -44, -577, -6, 66, -552, -164, 161, -460, -344, 205, -281, -504, ...
%!          134, -28, -576, -118, 156, -437, -381, 200, -220, -540, 83, 11, ...
%!         -568, -160, 172, -414, -408, 188, -125, -572, -32, 139, -492, ...
%!         -321, 205, -262, -504, 142, -83, -574, 0, 48, -571, -106, 137, ...
%!         -501, -266, 190, -391, -406, 194, -186, -553, 83, -13, -577, -49, ...
%!          103, -515, -280, 201, 300, -506, 131, -45, -578, -80, 138, -462, ...
%!         -361, 201, -211, -554, 32, 74, -533, -235, 187, -372, -442, 182, ...
%!         -147, -566, 25, 68, -535, -244, 194, -351, -463, 174, -125, -570, ...
%!           15, 72, -550, -190, 172, -424, -385, 198, -218, -536, 96];
%! [h, p, stats] = runstest (data, median (data));
%! expected_h = 1;
%! expected_p = 0.008562;
%! expected_z = 2.6229;
%! assert (h, expected_h);
%! assert (p, expected_p, 1E-6);
%! assert (stats.z, expected_z, 1E-4);

%!shared x
%! x = [45, -60, 1.225, 55.4, -9 27];
%!test
%! [h, p, stats] = runstest (x);
%! assert (h, 0);
%! assert (p, 0.6, 1e-14);
%! assert (stats.nruns, 5);
%! assert (stats.n1, 3);
%! assert (stats.n0, 3);
%! assert (stats.z, 0.456435464587638, 1e-14);
%!test
%! [h, p, stats] = runstest (x, [], "method", "approximate");
%! assert (h, 0);
%! assert (p, 0.6481, 1e-4);
%! assert (stats.z, 0.456435464587638, 1e-14);
%!test
%! [h, p, stats] = runstest (x, [], "tail", "left");
%! assert (h, 0);
%! assert (p, 0.9, 1e-14);
%! assert (stats.z, 1.369306393762915, 1e-14);

%!error<runstest: X must be a vector a scalar values.> runstest (ones (2,20))
%!error<runstest: X must be a vector a scalar values.> runstest (["asdasda"])
%!error<runstest: V must be either a scalar number or> ...
%! runstest ([2 3 4 3 2 3 4], "updown")
%!error<runstest: invalid value for alpha.> ...
%! runstest ([2 3 4 3 2 3 4], [], "alpha", 0)
%!error<runstest: invalid value for alpha.> ...
%! runstest ([2 3 4 3 2 3 4], [], "alpha", [0.02 0.2])
%!error<runstest: invalid value for alpha.> ...
%! runstest ([2 3 4 3 2 3 4], [], "alpha", 1.2)
%!error<runstest: invalid value for alpha.> ...
%! runstest ([2 3 4 3 2 3 4], [], "alpha", -0.05)
%!error<runstest: invalid value for method.> ...
%! runstest ([2 3 4 3 2 3 4], [], "method", "some")
%!error<runstest: invalid value for tail.> ...
%! runstest ([2 3 4 3 2 3 4], [], "tail", "some")
%!error<runstest: invalid optional argument.> ...
%! runstest ([2 3 4 3 2 3 4], [], "option", "some")
