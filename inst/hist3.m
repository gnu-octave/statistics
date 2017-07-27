## Copyright (C) 2015 Carnë Draug <carandraug@octave.org>
##
## This program is free software; you can redistribute it and/or
## modify it under the terms of the GNU General Public License as
## published by the Free Software Foundation; either version 3 of the
## License, or (at your option) any later version.
##
## This program is distributed in the hope that it will be useful, but
## WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
## General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program; if not, see
## <http:##www.gnu.org/licenses/>.

## -*- texinfo -*-
## @deftypefn  {Function File} {} hist3 (@var{X})
## @deftypefnx {Function File} {} hist3 (@var{X}, @var{nbins})
## @deftypefnx {Function File} {} hist3 (@var{X}, @qcode{"Nbins"}, @var{nbins})
## @deftypefnx {Function File} {} hist3 (@var{X}, @var{centers})
## @deftypefnx {Function File} {} hist3 (@var{X}, @qcode{"Ctrs"}, @var{centers})
## @deftypefnx {Function File} {} hist3 (@var{X}, @qcode{"Edges"}, @var{edges})
## @deftypefnx {Function File} {[@var{N}, @var{C}] =} hist3 (@dots{})
## @deftypefnx {Function File} {} hist3 (@dots{}, @var{prop}, @var{val}, @dots{})
## @deftypefnx {Function File} {} hist3 (@var{hax}, @dots{})
## Produce bivariate (2D) histogram counts or plots.
##
## The elements to produce the histogram are taken from the Nx2 matrix
## @var{X}.  Any row with NaN values are ignored.  The actual bins can
## be configured in 3 different: number, centers, or edges of the bins:
##
## @table @asis
## @item Number of bins (default)
## Produces equally spaced bins between the minimum and maximum values
## of @var{X}.  Defined as a 2 element vector, @var{nbins}, one for each
## dimension.  Defaults to @code{[10 10]}.
##
## @item Center of bins
## Defined as a cell array of 2 monotonically increasing vectors,
## @var{centers}.  The width of each bin is determined from the adjacent
## values in the vector with the initial and final bin, extending to Infinity.
##
## @item Edge of bins
## Defined as a cell array of 2 monotonically increasing vectors,
## @var{edges}.  @code{@var{N}(i,j)} contains the number of elements
## in @var{X} for which:
##
## @itemize @w{}
## @item
## @var{edges}@{1@}(i) <= @var{X}(:,1) < @var{edges}@{1@}(i+1)
## @item
## @var{edges}@{2@}(j) <= @var{X}(:,2) < @var{edges}@{2@}(j+1)
## @end itemize
##
## The consequence of this definition is that values outside the initial
## and final edge values are ignored, and that the final bin only contains
## the number of elements exactly equal to the final edge.
##
## @end table
##
## The return values, @var{N} and @var{C}, are the bin counts and centers
## respectively.  These are specially useful to produce intensity maps:
##
## @example
## [counts, centers] = hist3 (data);
## imagesc (centers@{1@}, centers@{2@}, counts)
## @end example
##
## If there is no output argument, or if the axes graphics handle
## @var{hax} is defined, the function will plot a 3 dimensional bar
## graph.  Any extra property/value pairs are passed directly to the
## underlying surface object.
##
## @seealso{hist, histc, lookup, mesh}
## @end deftypefn

function [N, C] = hist3 (X, varargin)
  if (nargin < 1)
    print_usage ();
  endif

  next_argin = 1;
  should_draw = true;
  if (isaxes (X))
    hax = X;
    X = varargin{next_argin++};
  elseif (nargout == 0)
    hax = gca ();
  else
    should_draw = false;
  endif

  if (! ismatrix (X) || columns (X) != 2)
    error ("hist3: X must be a 2 columns matrix");
  endif

  method = "nbins";
  val = [10 10];
  if (numel (varargin) >= next_argin)
    this_arg = varargin{next_argin++};
    if (isnumeric (this_arg))
      method = "nbins";
      val = this_arg;
    elseif (iscell (this_arg))
      method = "ctrs";
      val = this_arg;
    elseif (numel (varargin) >= next_argin
            && any (strcmpi ({"nbins", "ctrs", "edges"}, this_arg)))
      method = tolower (this_arg);
      val = varargin{next_argin++};
    else
      next_argin--;
    endif
  endif

  have_centers = false;
  switch (tolower (method))
    case "nbins"
      [r_edges, c_edges] = edges_from_nbins (X, val);
    case "ctrs"
      have_centers = true;
      centers = val;
      [r_edges, c_edges] = edges_from_centers (val);
    case "centers"
      ## This was supported until 1.2.4 when the Matlab compatible option
      ## 'Ctrs' was added.
      persistent warned = false;
      if (! warned)
        warning ("hist3: option `centers' is deprecated.  Use `ctrs'");
      endif
      have_centers = true;
      centers = val;
      [r_edges, c_edges] = edges_from_centers (val);
    case "edges"
      if (! iscell (val) || numel (val) != 2
          || ! all (cellfun (@isvector, val)))
        error ("hist3: EDGES must be a cell array with 2 vectors");
      endif
      [r_edges] = vec (val{1}, 2);
      [c_edges] = vec (val{2}, 2);
      out_rows = any (X < [r_edges(1) c_edges(1)]
                      | X > [r_edges(end) c_edges(end)], 2);
      X(out_rows,:) = [];
    otherwise
      ## we should never get here...
      error ("hist3: invalid binning method `%s'", method);
  endswitch

  ## We only remove the NaN now, after having computed the bin edges,
  ## because the extremes from each column that define the edges may
  ## be paired with a NaN.  While such values do not appear on the
  ## histogram, they must still be used to compute the histogram
  ## edges.
  X(any (isnan (X), 2), :) = [];

  r_idx = lookup (r_edges, X(:,1), "l");
  c_idx = lookup (c_edges, X(:,2), "l");

  counts_size = [numel(r_edges) numel(c_edges)];
  counts = accumarray ([r_idx, c_idx], 1, counts_size);

  if (should_draw)
    counts = counts.';
    z = zeros ((size (counts) +1) *2);
    z(2:end-1,2:end-1) = kron (counts, ones (2, 2));
    ## Setting the values for the end of the histogram bin like this
    ## seems straight wrong but that's hwo Matlab plots look.
    y = [kron(c_edges, ones (1, 2)) (c_edges(end)*2-c_edges(end-1))([1 1])];
    x = [kron(r_edges, ones (1, 2)) (r_edges(end)*2-r_edges(end-1))([1 1])];
    mesh (hax, x, y, z, "facecolor", [.75 .85 .95], varargin{next_argin:end});
  else
    N = counts;
    if (isargout (2))
      if (! have_centers)
        C = {(r_edges + [diff(r_edges)([1:end end])]/ 2) ...
             (c_edges + [diff(c_edges)([1:end end])]/ 2)};
      else
        C = centers(:)';
        C{1} = vec (C{1}, 2);
        C{2} = vec (C{2}, 2);
      endif
    endif
  endif

endfunction

function [r_edges, c_edges] = edges_from_nbins (X, nbins)
  if (! isnumeric (nbins) || numel (nbins) != 2)
    error ("hist3: NBINS must be a 2 element vector");
  endif
  inits = min (X, [], 1);
  ends  = max (X, [], 1);
  ends -= (ends - inits) ./ vec (nbins, 2);

  ## If any histogram side has an empty range, then still make NBINS
  ## but then place that value at the centre of the centre bin so that
  ## they appear in the centre in the plot.
  single_bins = inits == ends;
  if (any (single_bins))
    inits(single_bins) -= (floor (nbins(single_bins) ./2)) + 0.5;
    ends(single_bins) = inits(single_bins) + nbins(single_bins) -1;
  endif

  r_edges = linspace (inits(1), ends(1), nbins(1));
  c_edges = linspace (inits(2), ends(2), nbins(2));
endfunction

function [r_edges, c_edges] = edges_from_centers (ctrs)
  if (! iscell (ctrs) || numel (ctrs) != 2 || ! all (cellfun (@isvector, ctrs)))
    error ("hist3: CTRS must be a cell array with 2 vectors");
  endif
  r_edges = vec (ctrs{1}, 2);
  c_edges = vec (ctrs{2}, 2);
  r_edges(2:end) -= diff (r_edges) / 2;
  c_edges(2:end) -= diff (c_edges) / 2;
endfunction

%!demo
%! X = [
%!    1    1
%!    1    1
%!    1   10
%!    1   10
%!    5    5
%!    5    5
%!    5    5
%!    5    5
%!    5    5
%!    7    3
%!    7    3
%!    7    3
%!   10   10
%!   10   10];
%! hist3 (X)

%!test
%! N_exp = [ 0  0  0  5 20
%!           0  0 10 15  0
%!           0 15 10  0  0
%!          20  5  0  0  0];
%!
%! n = 100;
%! x = [1:n]';
%! y = [n:-1:1]';
%! D = [x y];
%! N = hist3 (D, [4 5]);
%! assert (N, N_exp);

%!test
%! N_exp = [0  0  0  0  1
%!          0  0  0  0  1
%!          0  0  0  0  1
%!          1  1  1  1 93];
%!
%! n = 100;
%! x = [1:n]';
%! y = [n:-1:1]';
%! D = [x y];
%! C{1} = [1 1.7 3 4];
%! C{2} = [1:5];
%! N = hist3 (D, C);
%! assert (N, N_exp);

## bug 44987
%!test
%! D = [1 1; 3 1; 3 3; 3 1];
%! [c, nn] = hist3 (D, {0:4, 0:4});
%! exp_c = zeros (5);
%! exp_c([7 9 19]) = [1 2 1];
%! assert (c, exp_c);
%! assert (nn, {0:4, 0:4});

%!test
%! for i = 10
%!   assert (size (hist3 (rand (9, 2), "Edges", {[0:.2:1]; [0:.2:1]})), [6 6])
%! endfor

%!test
%! edge_1 = linspace (0, 10, 10);
%! edge_2 = linspace (0, 50, 10);
%! [c, nn] = hist3 ([1:10; 1:5:50]', "Edges", {edge_1, edge_2});
%! exp_c = zeros (10, 10);
%! exp_c([1 12 13 24 35 46 57 68 79 90]) = 1;
%! assert (c, exp_c);
%!
%! assert (nn{1}, edge_1 + edge_1(2)/2, eps*10^4)
%! assert (nn{2}, edge_2 + edge_2(2)/2, eps*10^4)

%!shared X
%! X = [
%!  5  2
%!  5  3
%!  1  4
%!  5  3
%!  4  4
%!  1  2
%!  2  3
%!  3  3
%!  5  4
%!  5  3];

%!test
%! N = zeros (10);
%! N([1 10 53 56 60 91 98 100]) = [1 1 1 1 3 1 1 1];
%! C = {(1.2:0.4:4.8), (2.1:0.2:3.9)};
%! assert (nthargout ([1 2], @hist3, X), {N C}, eps*10^3)

%!test
%! N = zeros (5, 7);
%! N([1 5 17 18 20 31 34 35]) = [1 1 1 1 3 1 1 1];
%! C = {(1.4:0.8:4.6), ((2+(1/7)):(2/7):(4-(1/7)))};
%! assert (nthargout ([1 2], @hist3, X, [5 7]), {N C}, eps*10^3)
%! assert (nthargout ([1 2], @hist3, X, "Nbins", [5 7]), {N C}, eps*10^3)

%!test
%! N = [0 1 0; 0 1 0; 0 0 1; 0 0 0];
%! C = {(2:5), (2.5:1:4.5)};
%! assert (nthargout ([1 2], @hist3, X, "Edges", {(1.5:4.5), (2:4)}), {N C})

%!test
%! N = [0 0 1 0 1 0; 0 0 0 1 0 0; 0 0 1 4 2 0];
%! C = {(1.2:3.2), (0:5)};
%! assert (nthargout ([1 2], @hist3, X, "Ctrs", C), {N C})
%! assert (nthargout ([1 2], @hist3, X, C), {N C})

%!test
%! [~, C] = hist3 (rand (10, 2), "Edges", {[0 .05 .15 .35 .55 .95],
%!                                         [-1 .05 .07 .2 .3 .5 .89 1.2]});
%! C_exp = {[ 0.025  0.1   0.25   0.45  0.75  1.15], ...
%!          [-0.475  0.06  0.135  0.25  0.4   0.695  1.045  1.355]};
%! assert (C, C_exp, eps*10^2)

## Test how handling of out of borders is different whether we are
## defining Centers or Edges.
%!test
%! Xv = repmat ([1:10]', [1 2]);
%!
%! ## Test Centers
%! assert (hist3 (Xv, "Ctrs", {1:10, 1:10}), eye (10))
%!
%! N_exp = eye (6);
%! N_exp([1 end]) = 3;
%! assert (hist3 (Xv, "Ctrs", {3:8, 3:8}), N_exp)
%!
%! N_exp = zeros (8, 6);
%! N_exp([1 2 11 20 29 38 47 48]) = [2 1 1 1 1 1 1 2];
%! assert (hist3 (Xv, "Ctrs", {2:9, 3:8}), N_exp)
%!
%! ## Test Edges
%! assert (hist3 (Xv, "Edges", {1:10, 1:10}), eye (10))
%! assert (hist3 (Xv, "Edges", {3:8, 3:8}), eye (6))
%! assert (hist3 (Xv, "Edges", {2:9, 3:8}), [zeros(1, 6); eye(6); zeros(1, 6)])
%!
%! N_exp = zeros (14);
%! N_exp(3:12, 3:12) = eye (10);
%! assert (hist3 (Xv, "Edges", {-1:12, -1:12}), N_exp)
%!
%! ## Test for Nbins
%! assert (hist3 (Xv), eye (10))
%! assert (hist3 (Xv, [10 10]), eye (10))
%! assert (hist3 (Xv, "nbins", [10 10]), eye (10))
%! assert (hist3 (Xv, [5 5]), eye (5) * 2)
%!
%! N_exp = zeros (7, 5);
%! N_exp([1 9 10 18 26 27 35]) = [2 1 1 2 1 1 2];
%! assert (hist3 (Xv, [7 5]), N_exp)

%!test # bug #51059
%! D = [1 1; NaN 2; 3 1; 3 3; 1 NaN; 3 1];
%! [c, nn] = hist3 (D, {0:4, 0:4});
%! exp_c = zeros (5);
%! exp_c([7 9 19]) = [1 2 1];
%! assert (c, exp_c)
%! assert (nn, {0:4, 0:4})

## Single row of data or cases where all elements have the same value
## on one side of the histogram.
%!test
%! [c, nn] = hist3 ([1 8]);
%! exp_c = zeros (10, 10);
%! exp_c(6, 6) = 1;
%! exp_nn = {-4:5, 3:12};
%! assert (c, exp_c)
%! assert (nn, exp_nn, eps)
%!
%! [c, nn] = hist3 ([1 8], [10 11]);
%! exp_c = zeros (10, 11);
%! exp_c(6, 6) = 1;
%! exp_nn = {-4:5, 3:13};
%! assert (c, exp_c)
%! assert (nn, exp_nn, eps)

## NaNs paired with values defining the histogram edges.
%!test
%! [c, nn] = hist3 ([1 NaN; 2 3; 6 9; 8 NaN]);
%! exp_c = zeros (10, 10);
%! exp_c(2, 1) = 1;
%! exp_c(8, 10) = 1;
%! exp_nn = {linspace(1.35, 7.65, 10) linspace(3.3, 8.7, 10)};
%! assert (c, exp_c)
%! assert (nn, exp_nn, eps*100)

## Columns full of NaNs (recent Matlab versions seem to throw an error
## but this did work like this on R2010b at least).
%!test
%! [c, nn] = hist3 ([1 NaN; 2 NaN; 6 NaN; 8 NaN]);
%! exp_c = zeros (10, 10);
%! exp_nn = {linspace(1.35, 7.65, 10) NaN(1, 10)};
%! assert (c, exp_c)
%! assert (nn, exp_nn, eps*100)

## Behaviour of an empty X after removal of rows with NaN.
%!test
%! [c, nn] = hist3 ([1 NaN; NaN 3; NaN 9; 8 NaN]);
%! exp_c = zeros (10, 10);
%! exp_nn = {linspace(1.35, 7.65, 10) linspace(3.3, 8.7, 10)};
%! assert (c, exp_c)
%! assert (nn, exp_nn, eps*100)
