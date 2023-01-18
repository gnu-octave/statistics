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
## @deftypefn  {statistics} @var{B} = fillmissing (@var{A}, 'constant', @var{v})
## @deftypefnx {statistics} @var{B} = fillmissing (@var{A}, @var{method})
## @deftypefnx {statistics} @var{B} = fillmissing (@var{A}, @var{move_method}, @var{window_size})
## @deftypefnx {statistics} @var{B} = fillmissing (@var{A}, @var{fill_function}, @var{window_size})
## @deftypefnx {statistics} @var{B} = fillmissing (@dots{}, @var{dim})
## @deftypefnx {statistics} @var{B} = fillmissing (@dots{}, @var{PropertyName}, @var{PropertyValue})
## @deftypefnx {statistics} [@var{B}, @var{idx}] = fillmissing (@dots{})
##
## Replace missing entries of array @var{A} either with values in @var{v} or
## as determined by other specified methods.  'missing' values are determined
## by the data type of @var{A} as identified by the function @ref{ismissing},
## curently defined as:
##
## @itemize
## @item
## @code{NaN}: @code{single}, @code{double}
##
## @item
## @code{' '} (white space): @code{char}
##
## @item
## @code{@{''@}} (white space in cell): string cells.
## @end itemize
##
## @var{A} can be a numeric scalar or array, a character vector or array, or
## a cell array of character vectors (a.k.a. string cells).
##
## @var{v} can be a scalar or an array containing values for replacing the
## missing values in @var{A} with a compatible data type for isertion into
## @var{A}. The shape of @var{v} must be a scalar or an array with number
## of elements in @var{v} equal to the number of elements orthoganal to the
## operating dimension. E.g., if @code{size(@var{A})} = [3 5 4], operating
## along @code{dim} = 2 requires @var{v} to contain either 1 or 3x4=12
## elements.
##
## If requested, the optional output @var{idx} will contain a logical array
## the same shape as @var{A} indicating with 1's which locations in @var{A}
## were filled.
##
## Alternate Input Arguments and Values:
## @itemize
## @item @var{method} - replace missing values with:
## @table @code
##
## @item next
## @itemx previous
## @itemx nearest
## next, previous, or nearest non-missing value (nearest defaults to next
## when equidistant as determined by @code{SamplePoints}.)
##
## @item linear
## linear interpolation of neigboring, non-missing values
##
## @item spline
## piecewise cubic spline interpolation of neigboring, non-missing values
##
## @item pchip
## 'shape preserving' piecewise cubic spline interposaliton of neighboring,
## non-missing values
## @end table
##
## @item @var{move_method} - moving window calculated replacement values:
## @table @code
##
## @item movmean
## @itemx movmedian
## moving average or median using a window determined by @var{window_size}.
## @var{window_size} must be either a positive scalar value or a two element
## positive vector of sizes @w{@code{[@var{nb}, @var{na}]}} measured in the
## same units as @code{SamplePoints}.  For scalar values, the window is
## centered on the missing element and includes all data points within a
## distance of half of @var{window_size} on either side of the window center
## point.  Note that for compatability, when using a scalar value, the backward
## window limit is inclusive and the forward limit is exclusive.  If a
## two-element @var{window_size} vector is specified, the window includes all
## points within a distance of @var{nb} backward and @var{na} forward from the
## current element at the window center (both limits inclusive).
## @end table
##
## @item @var{fill_function} - custom method specified as a function handle.
## The supplied fill function must accept three inputs in the following order
## for each missing gap in the data:
## @table @var
## @item A_values -
## elements of @var{A} within the window on either side of the gap as
## determined by @var{window_size}.  (Note these elements can include missing
## values from other nearby gaps.)
## @item A_locs -
## locations of the reference data, @var{A_values}, in terms of the default
## or specified @code{SamplePoints}.
## @item gap_locs -
## location of the gap data points that need to be filled in terms of the
## default or specified @code{SamplePoints}.
## @end table
##
## The supplied function must return a scalar or vector with the same number of
## elements in @var{gap_locs}.  The required @var{window_size} parameter
## follows similar rules as for the moving average and median methods
## described above, with the two exceptions that (1) each gap is processed as a
## single element, rather than gap elements being processed individually, and
## (2) the window extended on either side of the gap has inclusive endpoints
## regardless of how @var{window_size} is specified.
##
## @item @var{dim} - specify a dimension for vector operation (default =
## first non-singeton dimension)
##
## @item @var{PropertyName}-@var{PropertyValue} pairs
## @table @code
## @item SamplePoints
## @var{PropertyValue} is a vector of sample point values representing the
## sorted and unique x-axis values of the data in @var{A}.  If unspecified,
## the default is assumed to be the vector @var{[1 : size (A, dim)]}.  The
## values in @code{SamplePoints} will affect methods and properties that rely
## on the effective distance between data points in @var{A}, such as
## interpolants and moving window functions where the @var{window_size}
## specified for moving window functions is measured relative to the
## @code{SamplePoints}.
##
## @item EndValues
## Apply a separate handling method for missing values at the front or back of
## the array. @var{PropertyValue} can be:
## @itemize
## @item A constant scalar or array with the same shape requirments as @var{v}.
## @item @code{none} - Do not fill end gap values.
## @item @code{extrap} - Use the same procedure as @var{method} to fill the
## end gap values.
## @item Any valid @var{method} listed above except for @code{movmean},
## @code{movmedian}, and @code{fill_function}. Those methods can only be
## applied to end gap values with @code{extrap}.
## @end itemize
##
## @item MissingLocations
## @var{PropertyValue} must be a logical array the same size as @var{A}
## indicating locations of known missing data with a value of @code{true}.
## (cannot be combined with MaxGap)
##
## @item MaxGap
## @var{PropertyValue} is a numeric scalar indicating the maximum gap length
## to fill, and assumes the same distance scale as the sample points. Gap
## length is calculated by the difference in locations of the sample points
## on either side of the gap, and gaps larger than MaxGap are ignored by
## @var{fillmissing}. (cannot be combined with MissingLocations)
## @end table
## @end itemize
##
## Compatibility Notes:
## @itemize
## @item
## Numerical and logical inputs for @var{A} and @var{v} may be specified
## in any combination. The output will be the same class as @var{A}, with the
## @var{v} converted to that data type for filling.  Only @code{single} and
## @code{double} have defined 'missing' values, so except for when the
## @code{missinglocations} option specifies the missing value identification of
## logical and other numeric data types, the output will always be
## @code{@var{B} = @var{A}} with @code{@var{idx} = false(size(@var{A}))}.
## @item
## All interpolation methods can be individually applied to @code{EndValues}.
## @item
## @sc{Matlab}'s @var{fill_function} method currently has several
## inconsistencies with the other methods (tested against version 2022a), and
## Octave's implementation has chosen the following consistent behavior over
## compatibility:  (1) a column full of missing data is considered part of
## @code{EndValues}, (2) such columns are then excluded from
## @var{fill_function} processing because the moving window is always empty.
## (3) operation in dimensions higher than 2 perform identically to operations
## in dims 1 and 2, most notable on vectors.
## @item
## Method "makima" is not yet implemented in @code{interp1}, which is
## used by @code{fillmissing}. Attempting to call this method will produce
## an error until the method is implemented in @code{interp1}.
## @end itemize
##
## @seealso{ismissing, rmmissing, standardizeMissing}
## @end deftypefn

function [A, idx_out] = fillmissing (A, varargin)

  if (nargin < 2)|| (nargin > 12)
     print_usage ();
  endif

  method = varargin{1};

  if (ischar (method))
    method = lower (method);
  elseif (! is_function_handle (method))
    error ("fillmissing: second input must be a string or function handle");
  endif

  sz_A = size (A);
  ndims_A = numel (sz_A);

  dim = [];
  missing_locs = [];
  endgap_method = [];
  endgap_locs = [];
  endgap_val = [];
  fill_vals = [];
  idx_flag = (nargout > 1);
  maxgap = [];
  missinglocations = false;
  reshape_flag = false;
  samplepoints = [];
  standard_samplepoints = true;
  v = [];

  if (idx_flag)
      idx_out = false (sz_A);
  endif

  ## process input arguments
  if (is_function_handle (method))
    ## verify function handle and window
    if ((nargin < 3) || ! isnumeric (varargin{2}) || ...
          ! any( numel (varargin{2})==[1 2]))
      error (["fillmissing: fill function handle must be followed by ", ...
              "a numeric scalar or two-element vector window size"]);
    elseif (nargin (method) < 3)
      error ("fillmissing: fill function must accept at least three inputs.");
    endif
    move_fcn = method;
    method = "movfcn";
    window_size = varargin{2};
    next_varg = 3;

  else
    switch (method)
      case {"previous", "next", "nearest"}
        next_varg = 2;

      case {"linear", "spline", "pchip", "makima"}
        next_varg = 2;
        if (! (isnumeric (A) || islogical (A)))
          error (["fillmissing: interpolation methods only valid for ", ...
                       "numeric input types"]);
        endif

      case "constant"
        if ((nargin < 3))
          error (["fillmissing: 'constant' method must be followed by a ", ...
                      "numeric scalar or array"]);
        endif

        v = varargin{2};
        if ! (isscalar (v) || isempty (v))
          v = v(:);
        endif

        if ((! ischar(v)) && isempty (v))
          error ("fillmissing: a numeric fill value cannot be emtpy");
        endif

        ## type check v against A
        if (iscellstr (A) && ischar (v) && ! iscellstr (v))
          v = {v};
        endif

        if ((! isempty (v)) && ...
             ((isnumeric (A) && ! (isnumeric (v) || islogical (v))) || ...
                 (ischar (A) && ! ischar (v)) || ...
                   (iscellstr (A) && ! (iscellstr (v)))))
          error ("fillmissing: fill value must be the same data type as 'A'");
        endif

        ## v can't be size checked until after processing rest of inputs
        next_varg = 3;

      case {"movmean", "movmedian"}
        if (! (isnumeric (A) || islogical (A)))
          error (["fillmissing: 'movmean' and 'movmedian' methods only ", ...
                             "valid for numeric input types"]);
        endif

        if ((nargin < 3) || ! isnumeric (varargin{2}) || ...
              ! any( numel (varargin{2})==[1 2]))
          error (["fillmissing: moving window method must be followed by ", ...
                   "a numeric scalar or two-element vector"]);
        endif
        window_size = varargin{2};
        next_varg = 3;

      otherwise
        error ("fillmissing: unknown fill method '%s'", method);
    endswitch
  endif

  ## process any more parameters
  if (next_varg < nargin)

    #set dim. if specified, it is the only numeric option that can appear next
    if isnumeric (varargin{next_varg})
      dim = varargin{next_varg};
      if (! (isscalar (dim) && (dim > 0)))
        error ("fillmissing: DIM must be a positive scalar");
      endif
      next_varg++;
    else
      ## default dim is first nonsingleton dimension of A
      if isscalar (A)
        dim = 1;
      else
        dim = find (sz_A > 1, 1, "first");
      endif
    endif
    sz_A_dim = size (A, dim);

    ## process any remaining inputs, must be name-value pairs
    while (next_varg < nargin)

      propname = varargin{next_varg};
      if (next_varg + 1 == nargin)
        ## must be at least one more input with 1st containing value
        error ("fillmissing: properties must be given as name,value pairs");

      else
        propval = varargin{next_varg + 1};
        next_varg = next_varg + 2;

        if (! ischar (propname))
          error ("invalid parameter name specified");
        else
          propname = lower (propname);
        endif

        ## input validation for names and values
        switch (propname)
          case "samplepoints"
            ## val must be sorted, unique, numeric vector the same size
            ## as size(A,dim)
            if (! (isnumeric (propval) && isvector (propval)
                 && (numel (propval) == sz_A_dim) && issorted (propval)
                 && (numel (propval) == numel (unique (propval)))))
              error (["fillmissing: SamplePoints must be a sorted ", ...
                       "non-repeating, numeric vector with %d elements"], ...
                         sz_A_dim);
            endif
            samplepoints = propval(:);
            standard_samplepoints = all (diff (samplepoints, 1, 1) == 1);

          case "endvalues"
            ## for numeric A, val must be numeric scalar, a numeric
            ## array with numel equal to the elements orthogonal to
            ## the dim or certain string methads. For non-numeric A,
            ## "constant" method is not valid.
            if ischar (propval)
              switch (lower (propval))
                case {"extrap", "previous", "next", "nearest", "none", ...
                       "linear", "spline", "pchip", "makima"}
                  endgap_method = propval;

                otherwise
                  error ("fillmissing: invalid EndValues method '%s'", propval);
              endswitch
            elseif (isnumeric (propval))
              if (! (isnumeric (A) || islogical (A)))
                error (["fillmissing: EndValues method 'constant' only ", ...
                          "valid for numeric arrays."]);
              endif
              endgap_method = 'constant';
              endgap_val = propval;

            else
              error (["fillmissing: EndValues must be numeric or a ", ...
                        "valid method name"]);
            endif

          case "missinglocations"

            if !(isnumeric (A) || islogical (A) || isinteger (A) || ...
                   ischar (A) || iscellstr (A))
              error (["fillmissing: MissingLocations option is not ", ...
                        "compatible with data type '%s'"], class (A));
            endif

            if (! isempty (maxgap))
              error (["fillmissing: MissingLocations and MaxGap options ", ...
                        "cannot be used simultaneously"]);
            endif

            ## val must be logical array same size as A
            if (! (islogical (propval) && isequal (sz_A, size (propval))))
              error (["fillmissing: MissingLocations must be a logical ", ...
                       "array the same size as A"]);
            endif

            missinglocations = true;
            missing_locs = propval;

          case "maxgap"
            ## val must be positive numeric scalar
            if (! (isnumeric (propval) && isscalar (propval) && (propval > 0)))
              error ("fillmissing: MaxGap must be a positive numeric scalar");
            endif

            if (! isempty (missing_locs))
              error (["fillmissing: MissingLocations and MaxGap options ", ...
                       "cannot be used simultaneously"]);
            endif

            maxgap = propval;
          case {"replacevalues", "datavariables"}
            error ("fillmissing: the '%s' option has not been implemented", ...
                                                                      propname);

        otherwise
          error ("invalid parameter name '%s'", propname);
        endswitch
      endif
    endwhile

  else
    ## no inputs after method
    ## set default dim
    if isscalar (A)
      dim = 1;
    else
      dim = find (sz_A > 1, 1, "first");
    endif
    sz_A_dim = size (A, dim);
  endif

  ## reduce calls to size and avoid overruns checking sz_A for high dims
  if dim > ndims_A
    sz_A = [sz_A, ones(1, dim - ndims_A)];
    ndims_A = numel (sz_A);
  endif

  ## set defaults for any unspecified parameters
  if (isempty (samplepoints))
    samplepoints = [1 : sz_A_dim]';
  endif
  if isempty (missing_locs)
    missing_locs = ismissing (A);
  endif

  ## endvalues treated separately from interior missing_locs
  if (isempty (endgap_method) || strcmp (endgap_method, "extrap"))
    endgap_method = method;
    if strcmp(endgap_method, "constant")
      endgap_val = v;
    endif
  endif


  ## missingvalues option not compatible with some methods and inputs:
  if (isinteger (A) || islogical (A))
    if (any (ismember (method, ...
             {"linear", "spline", "pchip", "makima", "movmean", "movmedian"})))
      error (["fillmissing: MissingLocations cannot be used ", ...
                "with method '%s' and inputs of type '%s'"], method, class (A));

    elseif (any (ismember (endgap_method, ...
                      {"linear", "spline", "pchip", "makima"})))
      error (["fillmissing: MissingLocations cannot be used with EndValues", ...
              " method '%s' and inputs of type '%s'"], method, class (A));
    endif
  endif


  ## verify size of v and endgap_val for 'constant' methods, resize for A
  orthogonal_size = [sz_A(1:dim-1), 1, sz_A(dim+1:end)]; # orthog. to dim size
  numel_orthogonal = prod (orthogonal_size); # numel perpen. to dim
  if (strcmp (method, "constant") && (! isscalar (v)))
    if (numel (v) != numel_orthogonal)
      error (["fillmissing: fill value 'V' must be a scalar or a %d ", ...
              " element array"], numel_orthogonal);
    else
      v = reshape (v, orthogonal_size);
    endif
  endif

  if (strcmp (endgap_method, "constant") && (! isscalar (endgap_val)))
    if (numel (endgap_val) != numel_orthogonal)
      error (["fillmissing: EndValues must be a scalar or a %d element ", ...
              "array"], numel_orthogonal);
    else
      endgap_val = reshape (endgap_val, orthogonal_size);
    endif
  endif

  ## simplify processing by temporarily permuting A so operation always on dim1
  ## revert permutation at the end
  dim_idx_perm = [1 : ndims_A];
  dim_idx_flip(1 : max(dim, ndims_A)) = {':'};
  dim_idx_flip(1) = [sz_A_dim:-1:1];

  if (dim != 1)
    dim_idx_perm([1, dim]) = [dim, 1];
    A = permute (A, dim_idx_perm);
    sz_A([1, dim]) = sz_A([dim, 1]);
    missing_locs = permute (missing_locs, dim_idx_perm);
    reshape_flag = true;
    orthogonal_size =  [1, sz_A(2:end)];
    if (idx_flag)
      idx_out = false (sz_A);
    endif
    if (! isempty (v) && ! isscalar (v))
      v = permute (v, dim_idx_perm);
    endif
    if (! isempty (endgap_val) && ! isscalar (endgap_val))
      endgap_val = permute (endgap_val, dim_idx_perm);
    endif
  endif

  ## precalculate fill data for several methods
  zero_padding = zeros (orthogonal_size);
  samplepoints_expand = samplepoints(:, ones (1, prod (sz_A(2:end))));

  ##find endgap locations
  if (sz_A_dim < 3)
    ## all missing are endgaps
    endgap_locs = missing_locs;
  else
    ## use cumsums starting from first and last part in dim to find missing
    ## values in and adjacent to end locations.
    endgap_locs = cumprod (missing_locs,1) | ...
                  (cumprod (missing_locs(dim_idx_flip{:}),1))(dim_idx_flip{:});
  endif
  ## remove endgap_locs from missing_locs to avoid double processing
  missing_locs(endgap_locs) = false;

  ## remove elements from missing and end location arrays if maxgap is specified
  if (! isempty (maxgap))
    ## missing_locs: if samplepoints value diff on either side of missing
    ## elements is > maxgap, remove those values.
    ## for endgaps, use diff of inside and missing end samplepoint values
    ## and remove from endgaps

    ## First check gapsize of any interior missings in missing_locs
    if (any (missing_locs(:)))
      ## locations in front of gaps
      loc_before = [diff(missing_locs,1,1); zero_padding] == 1;
      ## locations in back of gaps
      loc_after = diff ([zero_padding; missing_locs],1,1) == -1;

      ## value of samplepoints at front and back locs
      sampvals_before = samplepoints_expand(loc_before);
      sampvals_after = samplepoints_expand(loc_after);

      ## evaluate which gaps are too big to fill
      gaps_to_remove = (sampvals_after - sampvals_before) > maxgap;

      ## convert those gaps into an array element list
      idxs_to_remove = arrayfun ('colon', ...
                              ((find (loc_before))(gaps_to_remove ) + 1), ...
                               ((find (loc_after))(gaps_to_remove ) - 1), ...
                                 "UniformOutput", false);
      ## remove those elements from missing_locs
      missing_locs([idxs_to_remove{:}]) = false;

    endif

    ##then do any endgaps
    if (any (endgap_locs(:)))
      ## if any are all missing, remove for any value of maxgap
      endgap_locs &= ! prod (endgap_locs, 1);

      if ((sz_A_dim < 3) && (abs (samplepoints(2) - samplepoints(1)) > maxgap))
        ## shortcut - all missings are ends and exceed maxgap.
          endgap_locs(:) = false;
      else

        ## check gap size of front endgaps

        ##find loc element after gap
        nextvals =  sum (cumprod (endgap_locs,1)) + 1;
        ## compare diff between values at those points and at base with maxgap
        ends_to_remove = abs (samplepoints(nextvals) - samplepoints(1)) ...
                             > maxgap;
        ## remove any with gap>maxgap
        endgap_locs((cumprod (endgap_locs,1)) & ends_to_remove) = false;

        ## flip, repeat for back endgaps, then unflip and remove.
        nextvals =  sum (cumprod (endgap_locs(dim_idx_flip{:}),1)) + 1;
        ends_to_remove = abs (samplepoints(end:-1:1)(nextvals) ...
                                  - samplepoints(end)) > maxgap;
        endgap_locs((cumprod (...
           endgap_locs(dim_idx_flip{:}), 1)(dim_idx_flip{:})) & ...
                ends_to_remove) = false;
      endif
    endif
  endif

  if (any (strcmp (endgap_method, {"movmean", "movmedian", "movfcn"})))
    ## These methods only called for endgaps with "extrap", so endgaps
    ## are processed together in the missing_locs section.
    missing_locs = missing_locs | endgap_locs;
    endgap_locs(:) = false;
  endif


  ## Actaully fill the missing data

  ## process central missing values (all gaps bound by two valid datapoints)
  ## for each method, calcualte fill_vals, which will be used in assignment
  ## A(missing_locs) = fill_vals, and if idx_flag, populate idx_out
  if (any (missing_locs(:)))
    switch (method)
      case "constant"

        if (isscalar (v))
          fill_vals = v;
        else
          fill_vals = (missing_locs .* v)(missing_locs);
        endif

        if (idx_flag)
          ## if any v are the missing type, those get removed from idx_out
          ## unless using 'missinglocations'
          if (! missinglocations) && any (miss_v = ismissing (v))
            idx_out(missing_locs) = true;
            idx_out(missing_locs & miss_v) = false;

          else
            idx_out(missing_locs) = true;
          endif
        endif

      case {"previous", "next", "nearest", "linear"}
        ## find element locations bounding each gap
        loc_before = [diff(missing_locs, 1, 1); zero_padding] == 1;
        loc_after = diff ([zero_padding; missing_locs], 1, 1) == -1;
        gapsizes = find (loc_after) - find (loc_before) - 1;
        gap_count_idx = [1 : numel(gapsizes); gapsizes'];

        switch (method)
          case "previous"
            fill_vals = repelems (A(loc_before), gap_count_idx)';

          case "next"
            fill_vals = repelems (A(loc_after), gap_count_idx)';

          case {"nearest", "linear"}
            ## determine which missings go with values before or after
            ## gap based on samplevalue distance. (equal dist goes to after)

            ## find sample values before and after gaps
            sampvals_before = samplepoints_expand(loc_before);
            sampvals_after = samplepoints_expand(loc_after);

            ## build cell with linear indices of elements in each gap
            gap_locations = arrayfun ('colon', (find (loc_before)) + 1, ...
                               (find (loc_after)) - 1, "UniformOutput", false);

            ## get sample values at those elements
            [sampvals_in_gaps, ~] = ind2sub (sz_A, [gap_locations{:}]);
            sampvals_in_gaps = samplepoints(sampvals_in_gaps);

            ## expand first and last vectors for each gap point
            Avals_before = repelems (A(loc_before), gap_count_idx)';
            Avals_after = repelems (A(loc_after), gap_count_idx)';

            switch (method)
              case "nearest"
                ## calculate gap mid point for each gap element
                sampvals_midgap = repelems ( ...
                          (sampvals_before + sampvals_after)/2, gap_count_idx)';

                ## generate fill vectors sorting elements into nearest before
                ## or after
                prev_fill = (sampvals_in_gaps < sampvals_midgap);
                next_fill = (sampvals_in_gaps >= sampvals_midgap);
                fill_vals = A(missing_locs);
                fill_vals(prev_fill) = Avals_before(prev_fill);
                fill_vals(next_fill) = Avals_after(next_fill);

              case "linear"
                ## expand samplepoint values for interpolation x-values
                sampvals_before = repelems (sampvals_before, gap_count_idx)';
                sampvals_after = repelems (sampvals_after, gap_count_idx)';

                ## linearly interpolate
                fill_vals = ((Avals_after - Avals_before) ...
                                   ./ (sampvals_after - sampvals_before)) ...
                                   .* (sampvals_in_gaps - sampvals_before) ...
                                   + Avals_before;
            endswitch
        endswitch

        if (idx_flag)
          ## mid gaps will always be filled by above methods.
          idx_out(missing_locs) = true;
        endif

      case {"spline", "pchip", "makima"}
        ## pass more complex interpolations to interp1

        ## TODO: vectorized 'linear' is ~10-100x faster than using interp1.
        ## look to speed these up as well.

        ## identify columns needing interpolation to reduce empty operations
        cols_to_use = any (missing_locs, 1);

        ## missinglocations may send columns with NaN and less than 2
        ## real values resulting in interp1 error. Trim those columns,
        ## prepopulate fill_vals with NaN, mark as filled.

        if (missinglocations)
          fill_vals = NaN (sum (missing_locs(:, cols_to_use)(:)), 1);

          cols_enough_points =  (sum ( ...
                       !isnan(A) & (! missing_locs), 1) > 1) & cols_to_use;

          interp_cols = (cols_enough_points & cols_to_use);

          interp_vals = (missing_locs & cols_enough_points)(missing_locs & ...
                                                                   cols_to_use);

          fill_vals(interp_vals) = other_interpolants (A(:, interp_cols),
                  missing_locs(:, interp_cols), endgap_locs(:, interp_cols), ...
                          method, samplepoints);

        else
          fill_vals = other_interpolants (A(:, cols_to_use),
                  missing_locs(:, cols_to_use), endgap_locs(:, cols_to_use), ...
                          method, samplepoints);
        endif

        if (idx_flag)
          idx_out(missing_locs) = true;
        endif

      case {"movmean","movmedian"}
        ## check window size versus smallest sample gaps. if window smaller,
        ## nothing to do, break out early.
        if ((isscalar (window_size) && ...
                        (window_size/2 >= min (diff (samplepoints)))) || ...
             (isvector (window_size) &&
                       (sum (window_size) >= min (diff (samplepoints)))))

          switch (method)
            case "movmean"
              if sz_A_dim > 1
                allmissing = (missing_locs | endgap_locs)(:,:);

                ## create temporary flattened array for processing,
                A_sum = A(:,:);
                A_sum (allmissing) = 0;

                if (standard_samplepoints && ...
                        all (round (window_size) == window_size))
                  ## window size based on vector elements

                  ## faster codepath for uniform, unit-spacing samplepoints
                  ## and integer valued window sizes.

                  if (isscalar (window_size))
                    window_width = window_size;
                    if (mod (window_size, 2))
                      ## odd window size
                      ## equal number of values on either side of gap
                      window_size = (window_width - 1) .* [0.5, 0.5];
                    else
                      ## even window size
                      ## one extra element on previous side of gap
                      window_size(1) = window_width/2;
                      window_size(2) = window_size(1) - 1;
                    endif
                  else
                    window_width = window_size(1) + window_size(2) + 1;
                  endif

                  ## use columnwise convolution of windowing vector and A for
                  ## vectorized summation.
                  conv_vector = ones (window_width, 1);
                  A_sum = convn (A_sum, conv_vector, ...
                            "full")(1 + window_size(2):end - window_size(1), :);

                  ## get count of values contributing to convolution to account
                  ## for missing elements and to calculate mean.
                  A_sum_count = convn (! allmissing, conv_vector, ...
                            "full")(1 + window_size(2):end - window_size(1), :);

                else
                  ## window size based on sample point distance. Works for non
                  ## integer, non uniform values.

                  ## use A_sum (flattened to 2D), project slice windows in dim3
                  ## automatic broadcasting to get window summations & counts

                  samplepoints_shift = ...
                             samplepoints(:, ones(1, sz_A_dim)) - samplepoints';

                  if (isscalar (window_size))
                    ## [nb, na)
                    window_size = window_size * [-0.5, 0.5];
                    samplepoints_slice_windows = permute (...
                                samplepoints_shift >= window_size(1) & ...
                                  samplepoints_shift < window_size(2), [1,3,2]);
                  else
                    ## [nb, na]
                    window_size(1) = -window_size(1);
                    samplepoints_slice_windows = permute (...
                               samplepoints_shift >= window_size(1) & ...
                                 samplepoints_shift <= window_size(2), [1,3,2]);
                  endif


                  if (missinglocations)
                    ## NaNs left in A_sum will cause all sums to produce NaN

                    ## FIXME: when sum can handle nanflag, the 'else' path
                    ## should be able to be made to handle the vectorized
                    ## summation even with 'missinglocations'
                    A_nan = isnan (A_sum);
                    A_temp = A_sum .* samplepoints_slice_windows;
                    A_temp(!samplepoints_slice_windows & A_nan) = 0;
                    A_sum = permute (sum (A_temp, 1), [3,2,1]);

                  else
                    A_sum = permute (...
                       sum (A_sum .* samplepoints_slice_windows, 1), [3,2,1]);
                  endif

                  A_sum_count = permute (...
                        sum (! allmissing & samplepoints_slice_windows, 1), ...
                          [3,2,1]);
                endif

                ## build fill values
                fill_vals = A(missing_locs); # prefill to include missing vals
                fillable_gaps = missing_locs(:,:) & A_sum_count;
                fill_vals(fillable_gaps(missing_locs(:,:))) = ...
                           A_sum(fillable_gaps) ./ A_sum_count(fillable_gaps);
              endif

            case "movmedian"

              if sz_A_dim > 1

                cols_to_use = any (missing_locs(:,:), 1);

                samplepoints_shift = ...
                             samplepoints(:, ones(1, sz_A_dim)) - samplepoints';
                if (isscalar (window_size))
                  window_size = window_size * [-0.5, 0.5];

                  ## [nb, na)
                  samplepoints_slice_windows = permute (...
                                samplepoints_shift >= window_size(1) & ...
                                  samplepoints_shift < window_size(2), [1,3,2]);
                else
                  window_size(1) = -window_size(1);

                  ## [nb, na]
                  samplepoints_slice_windows = permute (...
                               samplepoints_shift >= window_size(1) & ...
                                 samplepoints_shift <= window_size(2), [1,3,2]);
                endif

                ## use moving window slices to project A and use
                ## custom function for vectorized full array median computation
                A_med = A(:, cols_to_use);
                nan_slice_windows = double (samplepoints_slice_windows);
                nan_slice_windows(! samplepoints_slice_windows) = NaN;
                A_med_slices = A_med .* nan_slice_windows;
                A_med = permute (columnwise_median (A_med_slices), [3 2 1]);

                fillable_gaps = missing_locs(:, cols_to_use);
                fill_vals = A_med(fillable_gaps);
              endif

          endswitch
          if (idx_flag)
            ## Matlab compatiblity - NaNs filled back in by movmean and
            ## movmedian should _not_ show as filled.
            idx_out(fillable_gaps) = true;
            still_nan = missing_locs;
            still_nan(missing_locs) = isnan (fill_vals);
            idx_out(still_nan) = false;
          endif

        endif

      case "movfcn"

        ## for each gap construct:
        ## xval - data values in window on either side of gap, including
        ## other missing values
        ## xloc - sample point values for those xval
        ## gap_loc - sample point values for gap elements
        ## if window has xval fully empty skip processing gap

        ## missing_locs might include endgap_locs
        ## need to build gap locations accounting for both types

        ## missing values can include more than just numeric inputs

        ## windows containing no data points (e.g., endgaps when window
        ## is one sided [3 0] or [0 2], will be dropped from processing,
        ## not being passed to the mov_fcn.

        if (isscalar (window_size))
          window_size = window_size * [-0.5, 0.5];
        else
          window_size(1) = -window_size(1);
        endif

        ## midgap bounds
        loc_before = [diff(missing_locs, 1, 1); zero_padding] == 1;
        loc_after = diff ([zero_padding; missing_locs], 1, 1) == -1;

        ## front/back endgap locations and bounds
        front_gap_locs = logical (cumprod (missing_locs, 1));
        front_next_locs = diff ([zero_padding; front_gap_locs], 1, 1) == -1;

        back_gap_locs = logical ( ...
                  cumprod (missing_locs(dim_idx_flip{:}), 1)(dim_idx_flip{:}));
        back_prev_locs = [diff(back_gap_locs, 1, 1); zero_padding] == 1;

        ## remove gap double counting
        back_gap_locs &= ! front_gap_locs;
        loc_before &= ! back_prev_locs;
        loc_after &= ! front_next_locs;

        ## build gap location array using gap starts and lengths.
        ## simplest to use front / mid / back ordering, track later with sort.
        gap_start_locs = ...
                [find(front_gap_locs & [true; false(sz_A_dim-1,1)])(:); ...
                    find(circshift (loc_before, 1, 1))(:);
                      find(circshift (back_prev_locs, 1, 1))(:)];

        gapsizes = [(sum (front_gap_locs, 1))(any (front_gap_locs, 1))(:);...
                       find(loc_after) - find(loc_before) - 1;...
                          (sum (back_gap_locs, 1))(any (back_gap_locs, 1))(:)];

        ## separate arrayfun/cellfun faster than single fun with
        ## composite anonymous function
        gap_locations = arrayfun ('colon', gap_start_locs, ...
                      gap_start_locs + gapsizes - 1, "UniformOutput", false);

        gap_locations = cellfun('transpose', ...
                                    gap_locations, "UniformOutput", false);

        ## sorting index to bridge front-mid-back and linear index ordering
        [~, gap_full_sort_idx] = sort (vertcat (gap_locations{:}));

        ## remove front or back gaps from gapsizes & gap_locations
        ## if front/back window size = 0, or if full column is missing

        ## index to track empty/removed elements
        removed_element_idx = true (numel (gap_full_sort_idx), 1);
        removed_front_elements = 0;
        removed_back_elements = 0;

        ## simple front/back gap trimming for either window size = 0

        if (! window_size(2))
          ## if no back facing window, ignore front gaps
          removed_front_gap_count = sum (front_gap_locs(1,:));
          removed_front_elements = sum (gapsizes(1 : removed_front_gap_count));
          removed_element_idx(1 : removed_front_elements) = false;
          gap_locations(1 : removed_front_gap_count ) = [];
          gapsizes(1 : removed_front_gap_count ) = [];

        elseif (any (missing_col_gaps = (gapsizes == sz_A_dim)))
          missing_col_elements = ...
            repelems (missing_col_gaps, [1:numel(gapsizes); gapsizes'])';
          removed_element_idx(missing_col_elements) = false;
          gap_locations(missing_col_gaps) = [];
          gapsizes(missing_col_gaps) = [];

        endif

        if (! window_size(1))
          ## if no front facing window, ignore back gaps.
          removed_back_gap_count = sum (back_gap_locs(end,:));
          removed_back_elements = sum (...
                            gapsizes(end - removed_back_gap_count + 1 : end));
          removed_element_idx(end - removed_back_elements + 1 : end) = false;
          gap_locations(end - removed_back_gap_count + 1 : end) = [];
          gapsizes(end - removed_back_gap_count + 1 : end) = [];
        endif

        if (! isempty (gapsizes))
          gap_sample_values = cellfun_subsref (gap_locations, false, ...
                                                     {samplepoints_expand});

          ## build [row,column] locations array for windows around each gap
          window_points_r_c = cell (numel (gapsizes), 2);
          window_points_r_c(:,1) = cellfun (@(x) ...
                ([1:sz_A_dim]')((samplepoints<x(1) & samplepoints >= ...
                    max (x(1) + window_size(1), samplepoints(1))) | ...
                     (samplepoints>x(end) & samplepoints <= ...
                       min (x(end) + window_size(2), samplepoints(end)))), ...
                         gap_sample_values, "UniformOutput", false);
          window_points_r_c(:,2) = cellfun ( ...
               @(x,y) (fix ((x(1)-1)/sz_A_dim)+1)(ones (size(y))), ...
                 gap_locations, window_points_r_c(:,1),"UniformOutput",false);


          ## if any window is emtpy, do not pass that gap to the move_fcn
          empty_gaps = cellfun ('isempty', window_points_r_c(:,1));
          if (any (empty_gaps))
            removed_element_idx(...
                  repelems (empty_gaps, [1:numel(gapsizes); gapsizes'])')...
                      = false;
            window_points_r_c(empty_gaps,:) = [];
            gap_sample_values(empty_gaps) = [];
            gapsizes(empty_gaps) = [];
            gap_locations(empty_gaps) = [];
          endif

          if (! isempty (gapsizes))
            ## Aval = A values at window locations
            ## Aloc = sample values at window locations
            A_window_indexes = cellfun ('sub2ind', {sz_A}, ...
                           window_points_r_c(:,1), window_points_r_c(:,2), ...
                               "UniformOutput", false);
            Aval = cellfun_subsref (A_window_indexes, false, {A});
            Aloc = cellfun_subsref (window_points_r_c(:,1), false, ...
                                                            {samplepoints});

            ##build fill values
            fill_vals_C = cellfun (move_fcn, Aval, Aloc, ...
                        gap_sample_values(:,1), "UniformOutput", false);

            ## check for output of move_fcn having different size than gaps
            if (! all (cellfun ('numel', fill_vals_C) == gapsizes))
              error (["fillmissing: fill function return values must be ", ...
                              "the same size as the gaps"]);
            endif

            [~, gap_trim_sort_idx] = sort (vertcat (gap_locations{:}));
            fill_vals_trim = cell2mat (fill_vals_C);
            if (! isempty (fill_vals_trim))
              fill_vals = A(missing_locs); # prefill to include missing vals
              fill_vals(removed_element_idx(gap_full_sort_idx)) = ...
                                fill_vals_trim(gap_trim_sort_idx);

            if (idx_flag)
              ## for movfcn with A of type having missing:
              ## any outputs still containing class's 'missing' values
              ## are counted as not filled in idx_out, even if the value
              ## was put there by the movfcn. This is true even if
              ## missinglocations is used. If missinglocations changed
              ## a value with no apparent change, it still shows up
              ## as filled.
              ## if A has no missing value (int or logical), then without
              ## missinglocations, idx_out is always empty. With
              ## missinglocations, compatible behavior is undefined as
              ## Matlab 2022a has an apparent bug producing a error message
              ## saying missinglocations with int/logical needs a method that
              ## incldues function handle. Expect behavior should match other
              ## methods, where any processed missing value should be marked
              ## as filled no matter the fill value.
              if ((isnumeric(A) && !isinteger(A)) ||
                      ischar (A) || iscellstr (A))
                idx_out(missing_locs) = ! ismissing (fill_vals);

              elseif (missinglocations)
                ## any missing_loc processed and not skipped must become true
                idx_out(missing_locs) = removed_element_idx(gap_full_sort_idx);
              endif
            endif

            endif
          endif
        endif
    endswitch

    if (! isempty (fill_vals))
      A(missing_locs) = fill_vals;
      fill_vals = [];
    endif

  endif

  ## process endgaps
  if (any (endgap_locs(:)))
    switch (endgap_method)
      case "none"
        endgap_locs(:) = false;

      case "constant"
        if (isscalar (endgap_val))
          fill_vals = endgap_val;
        else
          fill_vals = (endgap_locs .* endgap_val)(endgap_locs);
        endif
        if (idx_flag)
          ## if any v are the missing type, those get removed from idx_out
          ## unless using 'missinglocations'
          idx_out(endgap_locs) = true;
          if (! missinglocations) && any (miss_ev = ismissing (endgap_val))
            idx_out(endgap_locs & miss_ev) = false;
          endif
        endif

      case {"previous", "next", "nearest", "linear", ...
                 "spline", "pchip", "makima"}
        ## all of these methods require sz_A_dim >= 2. shortcut path otherwise
        if (sz_A_dim < 2)
          endgap_locs(:) = false;
        else

          switch (endgap_method)
            case "previous"
              ## remove any gaps at front of array, includes all-missing cols
              endgap_locs (logical (cumprod (endgap_locs,1))) = false;

              if (any (endgap_locs(:)))

                ##find locations of the 'prev' value to use for filling
                subsval_loc = [diff(endgap_locs, 1, 1); zero_padding] == 1;

                ##calculate the number of spots each 'prev' needs to fill
                gapsizes = (sum (endgap_locs, 1))(any (endgap_locs, 1))(:);

                ## construct substitution value vector
                fill_vals = repelems (A(subsval_loc), ...
                                        [1:numel(gapsizes); gapsizes'])';
              endif

            case "next"
              ## remove any gaps at back of array from endgap_locs
              ## includes any all-missing columns
              endgap_locs(logical (cumprod ( ...
                       endgap_locs(dim_idx_flip{:}), 1)(dim_idx_flip{:}))) ...
                          = false;
              if (any (endgap_locs(:)))
                ##find locations of the 'next' value to use for filling
                subsval_loc = diff ([zero_padding; endgap_locs],1,1) == -1;
                ##calculate the number of spots each 'next' needs to fill
                gapsizes = (sum (endgap_locs, 1))(any (endgap_locs, 1))(:);
                ## construct substitution value vector
                fill_vals = repelems (A(subsval_loc), ...
                                        [1:numel(gapsizes); gapsizes'])';
              endif

            case "nearest"

              ## remove any all-missing columns
              endgap_locs &= (! prod (endgap_locs, 1));

              if (any (endgap_locs(:)))
                ## find front end info
                front_gap_locs = logical (cumprod (endgap_locs, 1));
                front_next_loc = diff ( ...
                                    [zero_padding; front_gap_locs], 1, 1) == -1;
                front_gapsizes = (sum (front_gap_locs, 1))(any ...
                                                            (front_gap_locs,1));

                ## find back end info
                back_gap_locs = logical ( ...
                    cumprod (endgap_locs(dim_idx_flip{:}), 1)(dim_idx_flip{:}));
                back_prev_loc = [diff(back_gap_locs, 1, 1); zero_padding] == 1;
                back_gapsizes = (sum (back_gap_locs, 1))(any (back_gap_locs,1));

                ## combine into fill variables
                [~, fb_sort_idx] = sort ...
                                  ([find(front_gap_locs); find(back_gap_locs)]);
                fillval_loc = [find(front_next_loc); find(back_prev_loc)];
                gapsizes = [front_gapsizes; back_gapsizes];

                ## construct substitution value vector with sort order to mix
                ## fronts and backs in column order
                fill_vals = (repelems (A(fillval_loc), ...
                             [1:numel(gapsizes); gapsizes'])')(fb_sort_idx);
              endif


            case "linear"
              ## endgaps not guaranteed to have enough points to interpolate
              cols_to_use = (sum (! (missing_locs | endgap_locs), 1) > 1) ...
                               & any (endgap_locs, 1);
              interp_locs = ! (missing_locs | endgap_locs) & cols_to_use;
              endgap_locs &= cols_to_use;

              if (any (endgap_locs(:)))

                ## process front endgaps
                front_gap_locs = logical (cumprod (endgap_locs, 1));
                fill_vals_front = [];

                if (any (front_gap_locs(:)))
                  front_gapsizes = (sum (front_gap_locs, 1))(any ...
                                                            (front_gap_locs,1));
                  ## collect first data point after gap & expand to gapsize
                  front_interppoint_1 = repelems ( find (...
                        diff ([zero_padding; front_gap_locs], 1, 1) == -1), ...
                                   [1:numel(front_gapsizes); front_gapsizes'])';
                  ## collect second data point after gap & expand to gapsize
                  front_interppoint_2 = repelems ( find ( ...
                           diff ([zero_padding; ((cumsum (interp_locs, 1) .* ...
                               any (front_gap_locs, 1)) == 2)], 1, 1) == 1), ...
                                   [1:numel(front_gapsizes); front_gapsizes'])';
                  front_interp_Avals = A([front_interppoint_1, ...
                                            front_interppoint_2]);
                  front_interp_sampvals = samplepoints_expand( ...
                                    [front_interppoint_1, front_interppoint_2]);
                  front_gap_loc_sampvals = samplepoints_expand(front_gap_locs);

                  ## hack for vector automatic orientation forcing col vector
                  if (isvector (front_interp_Avals))
                    front_interp_Avals = (front_interp_Avals(:)).';
                  endif
                  if (isvector (front_interp_sampvals))
                    front_interp_sampvals = (front_interp_sampvals(:)).';
                  endif

                  ## perform interpolation for every gap point
                  interp_slopes_front = diff (front_interp_Avals, 1, 2) ...
                                          ./ diff (front_interp_sampvals, 1, 2);
                  fill_vals_front = interp_slopes_front .* ...
                                     (front_gap_loc_sampvals - ...
                                              front_interp_sampvals(:,1)) + ...
                                                    front_interp_Avals(:,1);
                endif

                ## process back endgaps
                back_gap_locs = logical ( ...
                    cumprod (endgap_locs(dim_idx_flip{:}), 1)(dim_idx_flip{:}));
                fill_vals_back = [];

                if (any (back_gap_locs(:)))
                  back_gapsizes = (sum ( ...
                                     back_gap_locs, 1))(any (back_gap_locs,1));

                  ## collect last data point before gap & expand to gapsize
                  back_interppoint_2 = repelems ( ...
                      find ([diff(back_gap_locs, 1, 1); zero_padding] == 1), ...
                                [1:numel(back_gapsizes); back_gapsizes'])';

                  ## collect 2nd to last data point before gap & expand to gap
                  back_interppoint_1 = repelems ( ...
                            find ((diff ([zero_padding; ...
                              ((cumsum (interp_locs(dim_idx_flip{:}), 1) .* ...
                                 any (back_gap_locs, 1)) == 2)], ...
                                   1, 1) == 1)(dim_idx_flip{:})), ...
                                     [1:numel(back_gapsizes); back_gapsizes'])';

                  ## build linear interpolant vectors
                  back_interp_Avals = A([back_interppoint_1, ...
                                          back_interppoint_2]);
                  back_interp_sampvals = samplepoints_expand( ...
                                      [back_interppoint_1, back_interppoint_2]);
                  back_gap_loc_sampvals = samplepoints_expand(back_gap_locs);

                  ## hack for vector automatic orientation forcing col vector
                  if (isvector (back_interp_Avals))
                    back_interp_Avals = (back_interp_Avals(:)).';
                  endif
                  if (isvector (back_interp_sampvals))
                    back_interp_sampvals = (back_interp_sampvals(:)).';
                  endif

                  ## perform interpolation for every gap point
                  interp_slopes_back = diff (back_interp_Avals, 1, 2) ...
                                          ./ diff (back_interp_sampvals, 1, 2);
                  fill_vals_back = interp_slopes_back .* ...
                                      (back_gap_loc_sampvals - ...
                                         back_interp_sampvals(:,1)) + ...
                                            back_interp_Avals(:,1);
                endif

                [~, fb_sort_idx] = sort ...
                                 ([find(front_gap_locs); find(back_gap_locs)]);
                fill_vals = [fill_vals_front; fill_vals_back](fb_sort_idx);
              endif

            case {"spline", "pchip", "makima"}
              ## endgap_locs not guaranteed to have 2 points.
              ## need to ignore columns with < 2 values, or with nothing to do
              cols_to_use = (sum (! (endgap_locs | missing_locs), 1) > 1) ...
                               & any (endgap_locs, 1);
              ## trim out unused cols from endgap_locs
              endgap_locs &= cols_to_use;

              if (missinglocations)
                ## missinglocations may send columns with NaN and less than 2
                ## real values resulting in interp1 error. Trim those columns,
                ## prepopulate fill_vals with NaN, mark as filled.

                fill_vals = NaN (sum (endgap_locs(:, cols_to_use)(:)), 1);

                cols_enough_points = (sum ( ...
                             !isnan(A) & (! endgap_locs), 1) > 1) & cols_to_use;

                interp_cols = (cols_enough_points & cols_to_use);

                interp_vals = (endgap_locs & ...
                                 cols_enough_points)(endgap_locs & cols_to_use);

                fill_vals(interp_vals) = other_interpolants ( ...
                              A(:, interp_cols),endgap_locs(:, interp_cols), ...
                                missing_locs(:, interp_cols), endgap_method, ...
                                  samplepoints);

              else

                fill_vals = other_interpolants (
                   A(:, cols_to_use), endgap_locs(:, cols_to_use), ...
                     missing_locs(:, cols_to_use), endgap_method, samplepoints);
              endif
          endswitch
        endif

        if (idx_flag)
          idx_out(endgap_locs) = true;
        endif
    endswitch

    ## some methods remove fill locations, only process if not empty
    if (any (endgap_locs(:)))
      ## replace missings with appropriate fill values
      A(endgap_locs) = fill_vals;
    endif
  endif

  if (reshape_flag)
    ## revert permutation
    A = permute (A, dim_idx_perm);
    if (idx_flag)
      idx_out = permute (idx_out, dim_idx_perm);
    endif
  endif

endfunction


function varargout = cellfun_subsref (x, TF, varargin)
  ## utility fcn for cellfun (@(x) A(x), x, "UniformOutput", true/false)
  ## ~50% faster than anonymous function call
  ## pass A, x, and truefalse.
  ## if nargaut > 1, repeat for C2 with B(x), C3 with C(x), etc.

  x_C = num2cell (struct ("type", "()", "subs", num2cell (x)));
  for (idx = 1 : numel (varargin))
    varargout{idx} = cellfun ('subsref', varargin{idx}, x_C, ...
                                                     "UniformOutput", TF);
  endfor
endfunction



function fill_vals = other_interpolants (data_array, primary_locs,
                                          secondary_locs, method, samplepoints)
  ## use interp1 to perform more complex interpolations. will only be performed
  ## on numerical data.
  ## primary_locs is missing_locs or endgap_locs, whichever the fill_vals are
  ## being returned for.  secondary_locs is the other.
  ##
  ## TODO: splitting out from columnwise cellfun to interp1 would increase
  ## speed a lot, but cannot count on same number of elements being processed
  ## in each column.
  ##
  ## Will error on any columns without at least two non-NaN interpolation
  ## values.
  ##
  ## Matlab incompatibility - using interp1, if 'missinglocations' sends through
  ## data_array values with NaN, interp1 will ignore those and interpolate with
  ## other data points.  Matlab will instead return NaN.

  ## find logical data and empty location indices for each column to be
  ## used in columnwise call to interp1 (cast as columnwise cell arrays)

  interp_data_locs = num2cell ((! (primary_locs | secondary_locs)), 1);
  interp_locs_tofill = num2cell (primary_locs, 1);

  ## build cell arrays with sample and interp values to pass to interp1
  [A_interpvalues, interp_samplevals] = cellfun_subsref (interp_data_locs, ...
                              false, num2cell (data_array, 1), {samplepoints});

  interp_empty_samplevals = cellfun_subsref (interp_locs_tofill, false, ...
                                                               {samplepoints});

  ## generate fill_vals using interp1 for missing locations.
  fill_vals = vertcat (cellfun ('interp1', interp_samplevals, ...
               A_interpvalues, interp_empty_samplevals, {method}, ...
                {'extrap'}, 'UniformOutput', false){:});
endfunction



function med = columnwise_median (x)
  ## takes a column of values, ignores any NaN values, retuns the median
  ## of what's left.  returns NaN if no values.
  ## uses only built-in fns to avoid 3x 'median' slowdown

  sz_x = size (x);

  if (isempty (x))
    med = NaN ([1, sz_x(2:end)]);

  elseif (isvector (x))
    x = x(! isnan (x));
    x = sort (x);
    n = numel (x);
    if (mod (n, 2))
      ## odd
      med = x((n+1)/2);
    elseif (n > 0)
      ## even
      med = (x(n/2) + x((n/2)+1))/2;
    else
      ## only called for types with missing = NaN
      med = NaN;
    endif

  else
    x = sort (x, 1); # NaNs sent to bottom
    n = sum (! isnan (x), 1);
    odd = logical (mod (n, 2)); # 0 even or zero, 1 odd
    even = !odd & (n != 0);

    if (ismatrix (x))
      med = NaN ([1, sz_x(2)]);
      med_idx_odd = sub2ind (sz_x, (n(odd)+1)/2, (1 : sz_x(2))(odd));
      med_idx_even = sub2ind (sz_x, [n(even)/2; n(even)/2 + 1], ...
                                  (1 : sz_x(2))(even)([1 1], :));

    else #nD arrays
      sz_x_flat = [sz_x(1), prod(sz_x(2:end))];
      med = NaN ([1, sz_x(2:end)]);
      med_idx_odd = sub2ind (sz_x_flat, ...
                      ((n(odd)+1)/2)(:)', (1 : sz_x_flat(2))(odd)(:)');
      med_idx_even = sub2ind (sz_x_flat, ...
                                [(n(even)/2)(:)'; (n(even)/2 + 1)(:)'], ...
                                   (1 : sz_x_flat(2))(even)([1 1], :));
    endif

    med(odd) = x(med_idx_odd);
    med(even) = sum (x(med_idx_even), 1) / 2;
  endif

endfunction


%!assert (fillmissing ([1 2 3], "constant", 99), [1 2 3])
%!assert (fillmissing ([1 2 NaN], "constant", 99), [1 2 99])
%!assert (fillmissing ([NaN 2 NaN], "constant", 99), [99 2 99])
%!assert (fillmissing ([1 2 3]', "constant", 99), [1 2 3]')
%!assert (fillmissing ([1 2 NaN]', "constant", 99), [1 2 99]')

%!assert (fillmissing ([1 2 3; 4 5 6], "constant", 99), [1 2 3; 4 5 6])
%!assert (fillmissing ([1 2 NaN; 4 NaN 6], "constant", 99), [1 2 99; 4 99 6])
%!assert (fillmissing ([NaN 2 NaN; 4 NaN 6], "constant", [97, 98, 99]), [97 2 99; 4 98 6])

%!test
%! x = cat (3, [1, 2, NaN; 4, NaN, 6], [NaN, 2, 3; 4, 5, NaN]);
%! y = cat (3, [1, 2, 99; 4, 99, 6], [99, 2, 3; 4, 5, 99]);
%! assert (fillmissing (x, "constant", 99), y);
%! y = cat (3, [1, 2, 96; 4, 95, 6], [97, 2, 3; 4, 5, 99]);
%! assert (fillmissing (x, "constant", [94:99]), y);
%! assert (fillmissing (x, "constant", [94:99]'), y);
%! assert (fillmissing (x, "constant", permute ([94:99], [1 3 2])), y);
%! assert (fillmissing (x, "constant", [94, 96, 98; 95, 97, 99]), y);
%! assert (fillmissing (x, "constant", [94:99], 1), y);
%! y = cat (3, [1, 2, 96; 4, 97, 6], [98, 2, 3; 4, 5, 99]);
%! assert (fillmissing (x, "constant", [96:99], 2), y);
%! y = cat (3, [1, 2, 98; 4, 97, 6], [94, 2, 3; 4, 5, 99]);
%! assert (fillmissing (x, "constant", [94:99], 3), y);
%! y = cat (3, [1, 2, 92; 4, 91, 6], [94, 2, 3; 4, 5, 99]);
%! assert (fillmissing (x, "constant", [88:99], 99), y);

%!test
%! x = reshape([1:24],4,3,2);
%! x([1, 6, 7, 9, 12, 14, 16, 19, 22, 23]) = NaN;
%! y = x;
%! y([1,6,7,9 12, 14, 16, 19, 22, 23]) = [94 95 95 96 96 97 97 98 99 99];
%! assert (fillmissing (x, "constant", [94:99], 1), y);
%! y([1,6,7,9 12, 14, 16, 19, 22, 23]) = [92 93 94 92 95 97 99 98 97 98];
%! assert (fillmissing (x, "constant", [92:99], 2), y);
%! y([1,6,7,9 12, 14, 16, 19, 22, 23]) = [88 93 94 96 99 89 91 94 97 98];
%! assert (fillmissing (x, "constant", [88:99], 3), y);
%! y([1,6,7,9 12, 14, 16, 19, 22, 23]) = [76 81 82 84 87 89 91 94 97 98];
%! assert (fillmissing (x, "constant", [76:99], 99), y);

## tests with different endvalues behavior
%!assert (fillmissing ([1 2 3], "constant", 99, "endvalues", 88), [1 2 3])
%!assert (fillmissing ([1 NaN 3], "constant", 99, "endvalues", 88), [1 99 3])
%!assert (fillmissing ([1 2 NaN], "constant", 99, "endvalues", 88), [1 2 88])
%!assert (fillmissing ([NaN 2 3], "constant", 99, "endvalues", 88), [88 2 3])
%!assert (fillmissing ([NaN NaN 3], "constant", 99, "endvalues", 88), [88 88 3])
%!assert (fillmissing ([1 NaN NaN], "constant", 99, "endvalues", 88), [1 88 88])
%!assert (fillmissing ([NaN 2 NaN], "constant", 99, "endvalues", 88), [88 2 88])
%!assert (fillmissing ([NaN 2 NaN]', "constant", 99, "endvalues", 88), [88 2 88]')
%!assert (fillmissing ([1 NaN 3 NaN 5], "constant", 99, "endvalues", 88), [1 99 3 99 5])
%!assert (fillmissing ([1 NaN NaN NaN 5], "constant", 99, "endvalues", 88), [1 99 99 99 5])
%!assert (fillmissing ([NaN NaN NaN NaN 5], "constant", 99, "endvalues", 88), [88 88 88 88 5])
%!assert (fillmissing ([1 NaN 3 4 NaN], "constant", 99, "endvalues", 88), [1 99 3 4 88])
%!assert (fillmissing ([1 NaN 3 4 NaN], "constant", 99, 1, "endvalues", 88), [1 88 3 4 88])
%!assert (fillmissing ([1 NaN 3 4 NaN], "constant", 99, 1, "endvalues", "extrap"), [1 99 3 4 99])

%!test
%! x = reshape ([1:24], 3, 4, 2);
%! y = x;
%! x([1,2,5,6,8,10,13,16,18,19,20,21,22]) = NaN;
%! y([1,2,5,6,10,13,16,18,19,20,21,22])= 88; y([8])=99;
%! assert (fillmissing (x, "constant", 99, "endvalues", 88), y);
%! assert (fillmissing (x, "constant", 99, 1, "endvalues", 88), y);
%! y = x; y([1,2,5,8,10,13,16,19,22])= 88; y([6,18,20,21])=99;
%! assert (fillmissing (x, "constant", 99, 2, "endvalues", 88), y);
%! y(y==99) = 88;
%! assert (fillmissing (x, "constant", 99, 3, "endvalues", 88), y);
%! assert (fillmissing (x, "constant", 99, 4, "endvalues", 88), y);
%! assert (fillmissing (x, "constant", 99, 99, "endvalues", 88), y);
%! y([8]) = 94;
%! assert (fillmissing (x, "constant", [92:99], 1, "endvalues", 88), y);
%! y([6,8,18,20,21]) = [96,88,99,98,99];
%! assert (fillmissing (x, "constant", [94:99], 2, "endvalues", 88), y);
%! y = x; y(isnan(y)) = 88;
%! assert (fillmissing (x, "constant", [88:99], 3, "endvalues", 88), y);
%! y = x; y(isnan(y)) = [82,82,83,83,94,85,86,87,87,88,88,88,89];
%! assert (fillmissing (x, "constant", [92:99], 1, "endvalues", [82:89]), y);
%! y = x; y(isnan(y)) = [84,85,85,96,85,84,87,87,99,87,98,99,87];
%! assert (fillmissing (x, "constant", [94:99], 2, "endvalues", [84:89]), y);
%! y = x; y(isnan(y)) = [68,69,72,73,75,77,68,71,73,74,75,76,77];
%! assert (fillmissing (x, "constant", [88:99], 3, "endvalues", [68:79]), y);
%! assert (fillmissing (x, "constant", [88:93;94:99]', 3, "endvalues", [68:73;74:79]'), y)

%!test
%! x = reshape([1:24],4,3,2);
%! x([1, 6, 7, 9, 12, 14, 16, 19, 22, 23]) = NaN;
%! y = x;
%! y([1,6,7,9 12, 14, 16, 19, 22, 23]) = [94 95 95 96 96 97 97 98 99 99];
%! assert (fillmissing (x, "constant", [94:99], 1), y);
%! y([1,6,7,9 12, 14, 16, 19, 22, 23]) = [92 93 94 92 95 97 99 98 97 98];
%! assert (fillmissing (x, "constant", [92:99], 2), y);
%! y([1,6,7,9 12, 14, 16, 19, 22, 23]) = [88 93 94 96 99 89 91 94 97 98];
%! assert (fillmissing (x, "constant", [88:99], 3), y);
%! y([1,6,7,9 12, 14, 16, 19, 22, 23]) = [76 81 82 84 87 89 91 94 97 98];
%! assert (fillmissing (x, "constant", [76:99], 99), y);

## next/previous tests
%!assert (fillmissing ([1 2 3], "previous"), [1 2 3])
%!assert (fillmissing ([1 2 3], "next"), [1 2 3])
%!assert (fillmissing ([1 2 3]', "previous"), [1 2 3]')
%!assert (fillmissing ([1 2 3]', "next"), [1 2 3]')
%!assert (fillmissing ([1 2 NaN], "previous"), [1 2 2])
%!assert (fillmissing ([1 2 NaN], "next"), [1 2 NaN])
%!assert (fillmissing ([NaN 2 NaN], "previous"), [NaN 2 2])
%!assert (fillmissing ([NaN 2 NaN], "next"), [2 2 NaN])
%!assert (fillmissing ([1 NaN 3], "previous"), [1 1 3])
%!assert (fillmissing ([1 NaN 3], "next"), [1 3 3])
%!assert (fillmissing ([1 2 NaN; 4 NaN 6], "previous", 1), [1 2 NaN; 4 2 6])
%!assert (fillmissing ([1 2 NaN; 4 NaN 6], "previous", 2), [1 2 2; 4 4 6])
%!assert (fillmissing ([1 2 NaN; 4 NaN 6], "previous", 3), [1 2 NaN; 4 NaN 6])
%!assert (fillmissing ([1 2 NaN; 4 NaN 6], "next", 1), [1 2 6; 4 NaN 6])
%!assert (fillmissing ([1 2 NaN; 4 NaN 6], "next", 2), [1 2 NaN; 4 6 6])
%!assert (fillmissing ([1 2 NaN; 4 NaN 6], "next", 3), [1 2 NaN; 4 NaN 6])

%!test
%! x = reshape([1:24],4,3,2);
%! x([1, 6, 7, 9, 12, 14, 16, 19, 22, 23]) = NaN;
%! y = x;
%! y([1, 6, 7, 9, 14, 19, 22, 23]) = [2 8 8 10 15 20 24 24];
%! assert (fillmissing (x, "next", 1), y);
%! y = x;
%! y([1, 6, 7, 14, 16]) = [5, 10, 11, 18, 20];
%! assert (fillmissing (x, "next", 2), y);
%! y = x;
%! y([1, 6, 9, 12]) = [13 18 21 24];
%! assert (fillmissing (x, "next", 3), y);
%! assert (fillmissing (x, "next", 99), x);
%! y = x;
%! y([6, 7, 12, 14, 16, 19, 22, 23]) = [5 5 11 13 15 18 21 21];
%! assert (fillmissing (x, "previous", 1), y);
%! y = x;
%! y([6, 7, 9, 12, 19, 22, 23]) = [2 3 5 8 15 18 15];
%! assert (fillmissing (x, "previous", 2), y);
%! y = x;
%! y([14, 16, 22, 23]) = [2 4 10 11];
%! assert (fillmissing (x, "previous", 3), y);
%! assert (fillmissing (x, "previous", 99), x);

## next/previous tests with different endvalue behavior
%!assert (fillmissing ([1 2 3], "constant", 0, "endvalues", "previous"), [1 2 3])
%!assert (fillmissing ([1 2 3], "constant", 0, "endvalues", "next"), [1 2 3])
%!assert (fillmissing ([1 NaN 3], "constant", 0, "endvalues", "previous"), [1 0 3])
%!assert (fillmissing ([1 NaN 3], "constant", 0, "endvalues", "next"), [1 0 3])
%!assert (fillmissing ([1 2 NaN], "constant", 0, "endvalues", "previous"), [1 2 2])
%!assert (fillmissing ([1 2 NaN], "constant", 0, "endvalues", "next"), [1 2 NaN])
%!assert (fillmissing ([1 NaN NaN], "constant", 0, "endvalues", "previous"), [1 1 1])
%!assert (fillmissing ([1 NaN NaN], "constant", 0, "endvalues", "next"), [1 NaN NaN])
%!assert (fillmissing ([NaN 2 3], "constant", 0, "endvalues", "previous"), [NaN 2 3])
%!assert (fillmissing ([NaN 2 3], "constant", 0, "endvalues", "next"), [2 2 3])
%!assert (fillmissing ([NaN NaN 3], "constant", 0, "endvalues", "previous"), [NaN NaN 3])
%!assert (fillmissing ([NaN NaN 3], "constant", 0, "endvalues", "next"), [3 3 3])
%!assert (fillmissing ([NaN NaN NaN], "constant", 0, "endvalues", "previous"), [NaN NaN NaN])
%!assert (fillmissing ([NaN NaN NaN], "constant", 0, "endvalues", "next"), [NaN NaN NaN])
%!assert (fillmissing ([NaN 2 NaN 4 NaN], "constant", 0, "endvalues", "previous"), [NaN 2 0 4 4])
%!assert (fillmissing ([NaN 2 NaN 4 NaN], "constant", 0, "endvalues", "next"), [2 2 0 4 NaN])
%!assert (fillmissing ([NaN 2 NaN 4 NaN], "constant", 0, 1, "endvalues", "previous"), [NaN 2 NaN 4 NaN])
%!assert (fillmissing ([NaN 2 NaN 4 NaN], "constant", 0, 1, "endvalues", "next"), [NaN 2 NaN 4 NaN])
%!assert (fillmissing ([NaN 2 NaN 4 NaN], "constant", 0, 2, "endvalues", "previous"), [NaN 2 0 4 4])
%!assert (fillmissing ([NaN 2 NaN 4 NaN], "constant", 0, 2, "endvalues", "next"), [2 2 0 4 NaN])
%!assert (fillmissing ([NaN 2 NaN 4 NaN], "constant", 0, 3, "endvalues", "previous"), [NaN 2 NaN 4 NaN])
%!assert (fillmissing ([NaN 2 NaN 4 NaN], "constant", 0, 3, "endvalues", "next"), [NaN 2 NaN 4 NaN])

%!test
%! x = reshape ([1:24], 3, 4, 2);
%! x([1,2,5,6,8,10,13,16,18,19,20,21,22]) = NaN;
%! y = x;
%! y([5,6,8,18])=[4,4,0,17];
%! assert (fillmissing (x, "constant", 0, "endvalues", "previous"), y);
%! assert (fillmissing (x, "constant", 0, 1, "endvalues", "previous"), y);
%! y = x;
%! y([6,10,18,20,21])=[0,7,0,0,0];
%! assert (fillmissing (x, "constant", 0, 2, "endvalues", "previous"), y);
%! y = x;
%! y([16,19,21])=[4,7,9];
%! assert (fillmissing (x, "constant", 0, 3, "endvalues", "previous"), y);
%! assert (fillmissing (x, "constant", 0, 4, "endvalues", "previous"), x);
%! assert (fillmissing (x, "constant", 0, 99, "endvalues", "previous"), x);
%! y = x;
%! y([1,2,8,10,13,16,22])=[3,3,0,11,14,17,23];
%! assert (fillmissing (x, "constant", 0, "endvalues", "next"), y);
%! assert (fillmissing (x, "constant", 0, 1, "endvalues", "next"), y);
%! y = x;
%! y([1,2,5,6,8,18,20,21])=[4,11,11,0,11,0,0,0];
%! assert (fillmissing (x, "constant", 0, 2, "endvalues", "next"), y);
%! y = x;
%! y([2,5])=[14,17];
%! assert (fillmissing (x, "constant", 0, 3, "endvalues", "next"), y);
%! assert (fillmissing (x, "constant", 0, 4, "endvalues", "next"), x);
%! assert (fillmissing (x, "constant", 0, 99, "endvalues", "next"), x);

##tests for nearest
%!assert (fillmissing ([1 2 3], "nearest"), [1 2 3])
%!assert (fillmissing ([1 2 3]', "nearest"), [1 2 3]')
%!assert (fillmissing ([1 2 NaN], "nearest"), [1 2 2])
%!assert (fillmissing ([NaN 2 NaN], "nearest"), [2 2 2])
%!assert (fillmissing ([1 NaN 3], "nearest"), [1 3 3])
%!assert (fillmissing ([1 2 NaN; 4 NaN 6], "nearest", 1), [1 2 6; 4 2 6])
%!assert (fillmissing ([1 2 NaN; 4 NaN 6], "nearest", 2), [1 2 2; 4 6 6])
%!assert (fillmissing ([1 2 NaN; 4 NaN 6], "nearest", 3), [1 2 NaN; 4 NaN 6])
%!assert (fillmissing ([1 NaN 3 NaN 5], "nearest"), [1 3 3 5 5])
%!assert (fillmissing ([1 NaN 3 NaN 5], "nearest", "samplepoints", [0 1 2 3 4]), [1 3 3 5 5])
%!assert (fillmissing ([1 NaN 3 NaN 5], "nearest", "samplepoints", [0.5 1 2 3 5]), [1 1 3 3 5])

%!test
%! x = reshape([1:24],4,3,2);
%! x([1, 6, 7, 9, 12, 14, 16, 19, 22, 23]) = NaN;
%! y = x;
%! y([1, 6, 7, 9, 12, 14, 16, 19, 22, 23]) = [2 5 8 10 11 15 15 20 21 24];
%! assert (fillmissing (x, "nearest", 1), y);
%! y = x;
%! y([1, 6, 7, 9, 12, 14, 16, 19, 22, 23]) = [5 10 11 5 8 18 20 15 18 15];
%! assert (fillmissing (x, "nearest", 2), y);
%! y = x;
%! y([1, 6, 9, 12, 14, 16, 22, 23]) = [13 18 21 24 2 4 10 11];
%! assert (fillmissing (x, "nearest", 3), y);
%! assert (fillmissing (x, "nearest", 99), x);

##tests for nearest with diff endvalue behavior
%!assert (fillmissing ([1 2 3], "constant", 0, "endvalues", "nearest"), [1 2 3])
%!assert (fillmissing ([1 NaN 3], "constant", 0, "endvalues", "nearest"), [1 0 3])
%!assert (fillmissing ([1 2 NaN], "constant", 0, "endvalues", "nearest"), [1 2 2])
%!assert (fillmissing ([1 NaN NaN], "constant", 0, "endvalues", "nearest"), [1 1 1])
%!assert (fillmissing ([NaN 2 3], "constant", 0, "endvalues", "nearest"), [2 2 3])
%!assert (fillmissing ([NaN NaN 3], "constant", 0, "endvalues", "nearest"), [3 3 3])
%!assert (fillmissing ([NaN NaN NaN], "constant", 0, "endvalues", "nearest"), [NaN NaN NaN])
%!assert (fillmissing ([NaN 2 NaN 4 NaN], "constant", 0, "endvalues", "nearest"), [2 2 0 4 4])
%!assert (fillmissing ([NaN 2 NaN 4 NaN], "constant", 0, 1, "endvalues", "nearest"), [NaN 2 NaN 4 NaN])
%!assert (fillmissing ([NaN 2 NaN 4 NaN], "constant", 0, 2, "endvalues", "nearest"), [2 2 0 4 4])
%!assert (fillmissing ([NaN 2 NaN 4 NaN], "constant", 0, 3, "endvalues", "nearest"), [NaN 2 NaN 4 NaN])

%!test
%! x = reshape ([1:24], 3, 4, 2);
%! x([1,2,5,6,8,10,13,16,18,19,20,21,22]) = NaN;
%! y = x;
%! y([1,2,5,6,8,10,13,16,18,22])=[3 3 4 4 0 11 14 17 17 23];
%! assert (fillmissing (x, "constant", 0, "endvalues", "nearest"), y);
%! assert (fillmissing (x, "constant", 0, 1, "endvalues", "nearest"), y);
%! y = x;
%! y([1,2,5,6,8,10,18,20,21])=[4 11 11 0 11 7 0 0 0];
%! assert (fillmissing (x, "constant", 0, 2, "endvalues", "nearest"), y);
%! y = x;
%! y([2,5,16,19,21])=[14 17 4 7 9];
%! assert (fillmissing (x, "constant", 0, 3, "endvalues", "nearest"), y);
%! assert (fillmissing (x, "constant", 0, 99, "endvalues", "nearest"), x);

##tests for linear
%!assert (fillmissing ([1 2 3], "linear"), [1 2 3])
%!assert (fillmissing ([1 2 3]', "linear"), [1 2 3]')
%!assert (fillmissing ([1 2 NaN], "linear"), [1 2 3])
%!assert (fillmissing ([NaN 2 NaN], "linear"), [NaN 2 NaN])
%!assert (fillmissing ([1 NaN 3], "linear"), [1 2 3])
%!assert (fillmissing ([1 2 NaN; 4 NaN 6], "linear", 1), [1 2 NaN; 4 NaN 6])
%!assert (fillmissing ([1 2 NaN; 4 NaN 6], "linear", 2), [1 2 3; 4 5 6])
%!assert (fillmissing ([1 2 NaN; 4 NaN 6], "linear", 3), [1 2 NaN; 4 NaN 6])
%!assert (fillmissing ([1 NaN 3 NaN 5], "linear"), [1 2 3 4 5])
%!assert (fillmissing ([1 NaN 3 NaN 5], "linear", "samplepoints", [0 1 2 3 4]), [1 2 3 4 5])
%!assert (fillmissing ([1 NaN 3 NaN 5], "linear", "samplepoints", [0 1.5 2 5 14]), [1 2.5 3 3.5 5], eps)

%!test
%! x = reshape([1:24],4,3,2);
%! x([1, 6, 7, 9, 12, 14, 16, 19, 22, 23]) = NaN;
%! assert (fillmissing (x, "linear", 1), reshape([1:24],4,3,2));
%! y = reshape([1:24],4,3,2);
%! y([1 9 14 19 22 23]) = NaN;
%! assert (fillmissing (x, "linear", 2), y);
%! y = reshape([1:24],4,3,2);
%! y([1, 6, 7, 9, 12, 14, 16, 19, 22, 23]) = NaN;
%! assert (fillmissing (x, "linear", 3), y);
%! assert (fillmissing (x, "linear", 99), x);

##tests for linear with diff endvalue behavior
%!assert (fillmissing ([1 2 3], "linear", "endvalues", 0), [1 2 3])
%!assert (fillmissing ([1 NaN 3], "linear", "endvalues", 0), [1 2 3])
%!assert (fillmissing ([1 2 NaN], "linear", "endvalues", 0), [1 2 0])
%!assert (fillmissing ([1 NaN NaN], "linear", "endvalues", 0), [1 0 0])
%!assert (fillmissing ([NaN 2 3], "linear", "endvalues", 0), [0 2 3])
%!assert (fillmissing ([NaN NaN 3], "linear", "endvalues", 0), [0 0 3])
%!assert (fillmissing ([NaN NaN NaN], "linear", "endvalues", 0), [0 0 0])
%!assert (fillmissing ([NaN 2 NaN 4 NaN], "linear", "endvalues", 0), [0 2 3 4 0])
%!assert (fillmissing ([NaN 2 NaN 4 NaN], "linear", 1, "endvalues", 0), [0 2 0 4 0])
%!assert (fillmissing ([NaN 2 NaN 4 NaN], "linear", 2, "endvalues", 0), [0 2 3 4 0])
%!assert (fillmissing ([NaN 2 NaN 4 NaN], "linear", 3, "endvalues", 0), [0 2 0 4 0])

%!test
%! x = reshape ([1:24], 3, 4, 2);
%! x([1,2,5,6,8,10,13,16,18,19,20,21,22]) = NaN;
%! y = x;
%! y([1,2,5,6,10,13,16,18,19,20,21,22])=0; y(8)=8;
%! assert (fillmissing (x, "linear", "endvalues", 0), y);
%! assert (fillmissing (x, "linear", 1, "endvalues", 0), y);
%! y = x;
%! y([1,2,5,8,10,13,16,19,22])=0; y([6,18,20,21])=[6,18,20,21];
%! assert (fillmissing (x, "linear", 2, "endvalues", 0), y);
%! y = x;
%! y(isnan(y))=0;
%! assert (fillmissing (x, "linear", 3, "endvalues", 0), y);
%! assert (fillmissing (x, "linear", 99, "endvalues", 0), y);

##tests with linear only on endvalues
%!assert (fillmissing ([1 2 3], "constant", 99, "endvalues", "linear"), [1 2 3])
%!assert (fillmissing ([1 NaN 3], "constant", 99, "endvalues", "linear"), [1 99 3])
%!assert (fillmissing ([1 NaN 3 NaN], "constant", 99, "endvalues", "linear"), [1 99 3 4])
%!assert (fillmissing ([NaN 2 NaN 4 NaN], "constant", 99, "endvalues", "linear"), [1 2 99 4 5])
%!assert (fillmissing ([NaN 2 NaN NaN], "constant", 99, "endvalues", "linear"), [NaN 2 NaN NaN])
%!assert (fillmissing ([NaN 2 NaN 4 NaN], "constant", 99, "endvalues", "linear", "samplepoints", [1 2 3 4 5]), [1 2 99 4 5])
%!assert (fillmissing ([NaN 2 NaN 4 NaN], "constant", 99, "endvalues", "linear", "samplepoints", [0 2 3 4 10]), [0 2 99 4 10])

##test other interpolants
%! x = reshape ([1:24], 3, 4, 2);
%! x([1,2,5,6,8,10,13,16,18,19,20,21,22]) = NaN;
%! y = x;
%! y([1,6,10,18,20,21]) = [2.5, 5, 8.5, 17.25, 21, 21.75];
%! assert (fillmissing (x, "linear", 2, "samplepoints", [2 4 8 10]), y, eps);
%! y([1,6,10,18,20,21]) = [2.5, 4.5, 8.5, 17.25, 21.5, 21.75];
%! assert (fillmissing (x, "spline", 2, "samplepoints", [2 4 8 10]), y, eps);
%! y([1,6,10,18,20,21]) = [2.5, 4.559386973180077, 8.5, 17.25, 21.440613026819925, 21.75];
%! assert (fillmissing (x, "pchip", 2, "samplepoints", [2 4 8 10]), y, 10*eps);

## known fail: makima method not yet implemented in interp1
%!test <60965>
%! x = reshape ([1:24], 3, 4, 2);
%! x([1,2,5,6,8,10,13,16,18,19,20,21,22]) = NaN;
%! y = x;
%! y([1,6,10,18,20,21]) = [2.5, 4.609523809523809, 8.5, 17.25, 21.390476190476186, 21.75];
%! assert (fillmissing (x, "makima", 2, "samplepoints", [2 4 8 10]), y, 10*eps);

##test other interpolants code path on endvalues
%!assert (fillmissing ([1 2 3], "constant", 99, "endvalues", "spline"), [1 2 3])
%!assert (fillmissing ([1 NaN 3], "constant", 99, "endvalues", "spline"), [1 99 3])
%!assert (fillmissing ([1 NaN 3 NaN], "constant", 99, "endvalues", "spline"), [1 99 3 4])
%!assert (fillmissing ([NaN 2 NaN 4 NaN], "constant", 99, "endvalues", "spline"), [1 2 99 4 5])
%!assert (fillmissing ([NaN 2 NaN NaN], "constant", 99, "endvalues", "spline"), [NaN 2 NaN NaN])
%!assert (fillmissing ([NaN 2 NaN 4 NaN], "constant", 99, "endvalues", "spline", "samplepoints", [1 2 3 4 5]), [1 2 99 4 5])
%!assert (fillmissing ([NaN 2 NaN 4 NaN], "constant", 99, "endvalues", "spline", "samplepoints", [0 2 3 4 10]), [0 2 99 4 10])


## test movmean
%!assert (fillmissing ([1 2 3], "movmean", 1), [1 2 3])
%!assert (fillmissing ([1 2 NaN], "movmean", 1), [1 2 NaN])
%!assert (fillmissing ([1 2 3], "movmean", 2), [1 2 3])
%!assert (fillmissing ([1 2 3], "movmean", [1 0]), [1 2 3])
%!assert (fillmissing ([1 2 3]', "movmean", 2), [1 2 3]')
%!assert (fillmissing ([1 2 NaN], "movmean", 2), [1 2 2])
%!assert (fillmissing ([1 2 NaN], "movmean", [1 0]), [1 2 2])
%!assert (fillmissing ([1 2 NaN], "movmean", [1 0]'), [1 2 2])
%!assert (fillmissing ([NaN 2 NaN], "movmean", 2), [NaN 2 2])
%!assert (fillmissing ([NaN 2 NaN], "movmean", [1 0]), [NaN 2 2])
%!assert (fillmissing ([NaN 2 NaN], "movmean", [0 1]), [2 2 NaN])
%!assert (fillmissing ([NaN 2 NaN], "movmean", [0 1.1]), [2 2 NaN])
%!assert (fillmissing ([1 NaN 3 NaN 5], "movmean", [3 0]), [1 1 3 2 5])
%!assert (fillmissing ([1 2 NaN; 4 NaN 6], "movmean", 3, 1), [1 2 6; 4 2 6])
%!assert (fillmissing ([1 2 NaN; 4 NaN 6], "movmean", 3, 2), [1 2 2; 4 5 6])
%!assert (fillmissing ([1 2 NaN; 4 NaN 6], "movmean", 3, 3), [1 2 NaN; 4 NaN 6])
%!assert (fillmissing ([1 NaN 3 NaN 5], "movmean", 99), [1 3 3 3 5])
%!assert (fillmissing ([1 NaN 3 NaN 5], "movmean", 99, 1), [1 NaN 3 NaN 5])
%!assert (fillmissing ([1 NaN 3 NaN 5]', "movmean", 99, 1), [1 3 3 3 5]')
%!assert (fillmissing ([1 NaN 3 NaN 5], "movmean", 99, 2), [1 3 3 3 5])
%!assert (fillmissing ([1 NaN 3 NaN 5]', "movmean", 99, 2), [1 NaN 3 NaN 5]')

%!assert (fillmissing ([1 NaN NaN NaN 5], "movmean", 3, "samplepoints", [1 2 3 4 5]), [1 1 NaN 5 5])
%!assert (fillmissing ([1 NaN NaN NaN 5], "movmean", [1 1], "samplepoints", [1 2 3 4 5]), [1 1 NaN 5 5])
%!assert (fillmissing ([1 NaN NaN NaN 5], "movmean", [1.5 1.5], "samplepoints", [1 2 3 4 5]), [1 1 NaN 5 5])
%!assert (fillmissing ([1 NaN NaN NaN 5], "movmean", 4, "samplepoints", [1 2 3 4 5]), [1 1 1 5 5])
%!assert (fillmissing ([1 NaN NaN NaN 5], "movmean", [2 2], "samplepoints", [1 2 3 4 5]), [1 1 3 5 5])
%!assert (fillmissing ([1 NaN NaN NaN 5], "movmean", 4.0001, "samplepoints", [1 2 3 4 5]), [1 1 3 5 5])
%!assert (fillmissing ([1 NaN NaN NaN 5], "movmean", 3, "samplepoints", [1.5 2 3 4 5]), [1 1 1 5 5])
%!assert (fillmissing ([1 NaN NaN NaN 5], "movmean", 3, "samplepoints", [1 2 3 4 4.5]), [1 1 NaN 5 5])
%!assert (fillmissing ([1 NaN NaN NaN 5], "movmean", 3, "samplepoints", [1.5 2 3 4 4.5]), [1 1 1 5 5])
%!assert (fillmissing ([1 NaN NaN NaN 5], "movmean", [1.5 1.5], "samplepoints", [1.5 2 3 4 5]), [1 1 1 5 5])
%!assert (fillmissing ([1 NaN NaN NaN 5], "movmean", [1.5 1.5], "samplepoints", [1 2 3 4 4.5]), [1 1 5 5 5])
%!assert (fillmissing ([1 NaN NaN NaN 5], "movmean", [1.5 1.5], "samplepoints", [1.5 2 3 4 4.5]), [1 1 3 5 5])

%!test
%! x = reshape ([1:24], 3, 4, 2);
%! x([1,2,5,6,8,10,13,16,18,19,20,21,22]) = NaN;
%! y = x;
%! y([2,5,8,10,13,16,18,22]) = [3,4,8,11,14,17,17,23];
%! assert (fillmissing (x, "movmean", 3), y);
%! assert (fillmissing (x, "movmean", [1 1]), y);
%! assert (fillmissing (x, "movmean", 3, "endvalues", "extrap"), y);
%! assert (fillmissing (x, "movmean", 3, "samplepoints", [1 2 3]), y);
%! y = x;
%! y([1,6,8,10,18,20,21]) = [4,6,11,7,15,20,24];
%! assert (fillmissing (x, "movmean", 3, 2), y);
%! assert (fillmissing (x, "movmean", [1 1], 2), y);
%! assert (fillmissing (x, "movmean", 3, 2, "endvalues", "extrap"), y);
%! assert (fillmissing (x, "movmean", 3, 2, "samplepoints", [1 2 3 4]), y);
%! y([1,18]) = NaN; y(6) = 9;
%! assert (fillmissing (x, "movmean", 3, 2, "samplepoints", [0 2 3 4]), y);
%! y = x;
%! y([1,2,5,6,10,13,16,18,19,20,21,22]) = 99; y(8) = 8;
%! assert (fillmissing (x, "movmean", 3, "endvalues", 99), y);
%! y = x;
%! y([1,2,5,8,10,13,16,19,22]) = 99; y([6,18,20,21]) = [6,15,20,24];
%! assert (fillmissing (x, "movmean", 3, 2, "endvalues", 99), y);


## test movmedian
%!assert (fillmissing ([1 2 3], "movmedian", 1), [1 2 3])
%!assert (fillmissing ([1 2 NaN], "movmedian", 1), [1 2 NaN])
%!assert (fillmissing ([1 2 3], "movmedian", 2), [1 2 3])
%!assert (fillmissing ([1 2 3], "movmedian", [1 0]), [1 2 3])
%!assert (fillmissing ([1 2 3]', "movmedian", 2), [1 2 3]')
%!assert (fillmissing ([1 2 NaN], "movmedian", 2), [1 2 2])
%!assert (fillmissing ([1 2 NaN], "movmedian", [1 0]), [1 2 2])
%!assert (fillmissing ([1 2 NaN], "movmedian", [1 0]'), [1 2 2])
%!assert (fillmissing ([NaN 2 NaN], "movmedian", 2), [NaN 2 2])
%!assert (fillmissing ([NaN 2 NaN], "movmedian", [1 0]), [NaN 2 2])
%!assert (fillmissing ([NaN 2 NaN], "movmedian", [0 1]), [2 2 NaN])
%!assert (fillmissing ([NaN 2 NaN], "movmedian", [0 1.1]), [2 2 NaN])
%!assert (fillmissing ([1 NaN 3 NaN 5], "movmedian", [3 0]), [1 1 3 2 5])
%!assert (fillmissing ([1 2 NaN; 4 NaN 6], "movmedian", 3, 1), [1 2 6; 4 2 6])
%!assert (fillmissing ([1 2 NaN; 4 NaN 6], "movmedian", 3, 2), [1 2 2; 4 5 6])
%!assert (fillmissing ([1 2 NaN; 4 NaN 6], "movmedian", 3, 3), [1 2 NaN; 4 NaN 6])
%!assert (fillmissing ([1 NaN 3 NaN 5], "movmedian", 99), [1 3 3 3 5])
%!assert (fillmissing ([1 NaN 3 NaN 5], "movmedian", 99, 1), [1 NaN 3 NaN 5])
%!assert (fillmissing ([1 NaN 3 NaN 5]', "movmedian", 99, 1), [1 3 3 3 5]')
%!assert (fillmissing ([1 NaN 3 NaN 5], "movmedian", 99, 2), [1 3 3 3 5])
%!assert (fillmissing ([1 NaN 3 NaN 5]', "movmedian", 99, 2), [1 NaN 3 NaN 5]')

%!assert (fillmissing ([1 NaN NaN NaN 5], "movmedian", 3, "samplepoints", [1 2 3 4 5]), [1 1 NaN 5 5])
%!assert (fillmissing ([1 NaN NaN NaN 5], "movmedian", [1 1], "samplepoints", [1 2 3 4 5]), [1 1 NaN 5 5])
%!assert (fillmissing ([1 NaN NaN NaN 5], "movmedian", [1.5 1.5], "samplepoints", [1 2 3 4 5]), [1 1 NaN 5 5])
%!assert (fillmissing ([1 NaN NaN NaN 5], "movmedian", 4, "samplepoints", [1 2 3 4 5]), [1 1 1 5 5])
%!assert (fillmissing ([1 NaN NaN NaN 5], "movmedian", [2 2], "samplepoints", [1 2 3 4 5]), [1 1 3 5 5])
%!assert (fillmissing ([1 NaN NaN NaN 5], "movmedian", 4.0001, "samplepoints", [1 2 3 4 5]), [1 1 3 5 5])
%!assert (fillmissing ([1 NaN NaN NaN 5], "movmedian", 3, "samplepoints", [1.5 2 3 4 5]), [1 1 1 5 5])
%!assert (fillmissing ([1 NaN NaN NaN 5], "movmedian", 3, "samplepoints", [1 2 3 4 4.5]), [1 1 NaN 5 5])
%!assert (fillmissing ([1 NaN NaN NaN 5], "movmedian", 3, "samplepoints", [1.5 2 3 4 4.5]), [1 1 1 5 5])
%!assert (fillmissing ([1 NaN NaN NaN 5], "movmedian", [1.5 1.5], "samplepoints", [1.5 2 3 4 5]), [1 1 1 5 5])
%!assert (fillmissing ([1 NaN NaN NaN 5], "movmedian", [1.5 1.5], "samplepoints", [1 2 3 4 4.5]), [1 1 5 5 5])
%!assert (fillmissing ([1 NaN NaN NaN 5], "movmedian", [1.5 1.5], "samplepoints", [1.5 2 3 4 4.5]), [1 1 3 5 5])

%!test
%! x = reshape ([1:24], 3, 4, 2);
%! x([1,2,5,6,8,10,13,16,18,19,20,21,22]) = NaN;
%! y = x;
%! y([2 5 8 10 13 16 18 22]) = [3 4 8 11 14 17 17 23];
%! assert (fillmissing (x, "movmedian", 3), y);
%! assert (fillmissing (x, "movmedian", [1 1]), y);
%! assert (fillmissing (x, "movmedian", 3, "endvalues", "extrap"), y);
%! assert (fillmissing (x, "movmedian", 3, "samplepoints", [1 2 3]), y);
%! y = x;
%! y([1 6 8 10 18 20 21]) = [4 6 11 7 15 20 24];
%! assert (fillmissing (x, "movmedian", 3, 2), y);
%! assert (fillmissing (x, "movmedian", [1 1], 2), y);
%! assert (fillmissing (x, "movmedian", 3, 2, "endvalues", "extrap"), y);
%! assert (fillmissing (x, "movmedian", 3, 2, "samplepoints", [1 2 3 4]), y);
%! y([1,18]) = NaN; y(6) = 9;
%! assert (fillmissing (x, "movmedian", 3, 2, "samplepoints", [0 2 3 4]), y);
%! y = x;
%! y([1,2,5,6,10,13,16,18,19,20,21,22]) = 99; y(8) = 8;
%! assert (fillmissing (x, "movmedian", 3, "endvalues", 99), y);
%! y = x;
%! y([1,2,5,8,10,13,16,19,22]) = 99; y([6,18,20,21]) = [6,15,20,24];
%! assert (fillmissing (x, "movmedian", 3, 2, "endvalues", 99), y);

## test movfcn
%!assert (fillmissing ([1 2 3], @(x,y,z) x+y+z, 2), [1 2 3])
%!assert (fillmissing ([1 2 NaN], @(x,y,z) x+y+z, 1), [1 2 NaN])
%!assert (fillmissing ([1 2 3], @(x,y,z) x+y+z, 2), [1 2 3])
%!assert (fillmissing ([1 2 3], @(x,y,z) x+y+z, [1 0]), [1 2 3])
%!assert (fillmissing ([1 2 3]', @(x,y,z) x+y+z, 2), [1 2 3]')
%!assert (fillmissing ([1 2 NaN], @(x,y,z) x+y+z, 2), [1 2 7])
%!assert (fillmissing ([1 2 NaN], @(x,y,z) x+y+z, [1 0]), [1 2 7])
%!assert (fillmissing ([1 2 NaN], @(x,y,z) x+y+z, [1 0]'), [1 2 7])
%!assert (fillmissing ([NaN 2 NaN], @(x,y,z) x+y+z, 2), [5 2 7])
%!assert (fillmissing ([NaN 2 NaN], @(x,y,z) x+y+z, [1 0]), [NaN 2 7])
%!assert (fillmissing ([NaN 2 NaN], @(x,y,z) x+y+z, [0 1]), [5 2 NaN])
%!assert (fillmissing ([NaN 2 NaN], @(x,y,z) x+y+z, [0 1.1]), [5 2 NaN])
%!assert (fillmissing ([1 2 NaN NaN 3 4], @(x,y,z) x+y+z, 2),[1 2 7 12 3 4])
%!assert (fillmissing ([1 2 NaN NaN 3 4], @(x,y,z) x+y+z, 0.5),[1 2 NaN NaN 3 4])

%!function A = testfcn (x,y,z)
%!  if isempty (y)
%!    A = z;
%!  elseif (numel (y) == 1)
%!    A = repelem (x(1), numel(z));
%!  else
%!    A = interp1 (y, x, z, "linear","extrap");
%!  endif
%!endfunction
%!assert (fillmissing ([1 NaN 3 NaN 5], @testfcn, [3 0]), [1 1 3 NaN 5])
%!assert (fillmissing ([1 2 NaN; 4 NaN 6], @testfcn, 3, 1), [1 2 6; 4 2 6])
%!assert (fillmissing ([1 2 NaN; 4 NaN 6], @testfcn, 3, 2), [1 2 2; 4 5 6])
%!assert (fillmissing ([1 2 NaN; 4 NaN 6], @testfcn, 3, 3), [1 2 NaN; 4 NaN 6])
%!assert (fillmissing ([1 NaN 3 NaN 5], @testfcn, 99), [1 2 3 4 5])
%!assert (fillmissing ([1 NaN 3 NaN 5], @testfcn, 99, 1), [1 NaN 3 NaN 5]) ##known not-compatible. matlab bug ML2022a: [1 1 3 1 5]
%!assert (fillmissing ([1 NaN 3 NaN 5]', @testfcn, 99, 1), [1 2 3 4 5]')
%!assert (fillmissing ([1 NaN 3 NaN 5], @testfcn, 99, 2), [1 2 3 4 5])
%!assert (fillmissing ([1 NaN 3 NaN 5]', @testfcn, 99, 2), [1 NaN 3 NaN 5]') ##known not-compatible. matlab bug ML2022a: [1 1 3 1 5]'
%!assert (fillmissing ([1 NaN 3 NaN 5], @testfcn, 99, 3), [1 NaN 3 NaN 5])
%!assert (fillmissing ([1 NaN 3 NaN 5]', @testfcn, 99, 3), [1 NaN 3 NaN 5]')
%!assert (fillmissing ([1 NaN NaN NaN 5], @testfcn, 3, "samplepoints", [1 2 3 4 5]), [1 2 3 4 5])
%!assert (fillmissing ([1 NaN NaN NaN 5], @testfcn, [1 1], "samplepoints", [1 2 3 4 5]), [1 2 3 4 5])
%!assert (fillmissing ([1 NaN NaN NaN 5], @testfcn, [1.5 1.5], "samplepoints", [1 2 3 4 5]), [1 2 3 4 5])
%!assert (fillmissing ([1 NaN NaN NaN 5], @testfcn, 4, "samplepoints", [1 2 3 4 5]), [1 2 3 4 5])
%!assert (fillmissing ([1 NaN NaN NaN 5], @testfcn, [2 2], "samplepoints", [1 2 3 4 5]), [1 2 3 4 5])
%!assert (fillmissing ([1 NaN NaN NaN 5], @testfcn, 3, "samplepoints", [1 2 2.5 3 3.5]), [1 2.6 3.4 4.2 5], 10*eps)
%!assert (fillmissing ([NaN NaN 3 NaN 5], @testfcn, 99, 1), [NaN NaN 3 NaN 5]) ##known not-compatible. matlab bug ML2022a: [1 1 3 1 5]

## known noncompatible. for move_fcn method, ML2021b (1) ignores windowsize
## for full missing column and processes it anyway, (2) doesn't consider it
## part of endvalues unlike all other methods, (3) ignores samplepoint values
##  when calcuating move_fcn results. should review against future versions.
%!test
%!function A = testfcn (x,y,z)
%!  if isempty (y)
%!    A = z;
%!  elseif (numel (y) == 1)
%!    A = repelem (x(1), numel(z));
%!  else
%!    A = interp1 (y, x, z, "linear","extrap");
%!  endif
%!endfunction
%! x = reshape ([1:24], 3, 4, 2);
%! x([1,2,5,6,8,10,13,16,18,19,20,21,22]) = NaN;
%! y = x;
%! y([1,2,5,6,8,10,13,16,18,22]) = [3,3,4,4,8,11,14,17,17,23];
%! assert (fillmissing (x, @testfcn, 3), y);
%! assert (fillmissing (x, @testfcn, [1 1]), y);
%! assert (fillmissing (x, @testfcn, 3, "endvalues", "extrap"), y);
%! assert (fillmissing (x, @testfcn, 3, "samplepoints", [1 2 3]), y);
%! y= x;
%! y(isnan(x)) = 99; y(8) = 8;
%! assert (fillmissing (x, @testfcn, 3, "endvalues", 99), y)
%! y = x;
%! y([1,2,5,6,8,10,18,20,21]) = [4,11,11,6,11,7,18,20,21];
%! assert (fillmissing (x, @testfcn, 3, 2), y);
%! assert (fillmissing (x, @testfcn, [1 1], 2), y);
%! assert (fillmissing (x, @testfcn, 3, 2, "endvalues", "extrap"), y);
%! assert (fillmissing (x, @testfcn, 3, 2, "samplepoints", [1 2 3 4]), y);
%! y(1) = NaN; y([6,18,21]) = [9,24,24];
%! assert (fillmissing (x, @testfcn, 3, 2, "samplepoints", [0 2 3 4]), y);
%! y = x;
%! y([1,2,5,6,10,13,16,18,19,20,21,22]) = 99; y(8) = [8];
%! assert (fillmissing (x, @testfcn, 3, "endvalues", 99), y);
%! y([6,18,20,21]) = [6,18,20,21]; y(8)=99;
%! assert (fillmissing (x, @testfcn, 3, 2, "endvalues", 99), y);
%! y([6,18,20,21]) = 99;
%! assert (fillmissing (x, @testfcn, 3, 3, "endvalues", 99), y);

##test maxgap for mid and end points
%!assert (fillmissing ([1 2 3], "constant", 0, "maxgap", 1), [1 2 3])
%!assert (fillmissing ([1 2 3], "constant", 0, "maxgap", 99), [1 2 3])
%!assert (fillmissing ([1 NaN 3], "constant", 0, "maxgap", 1), [1 NaN 3])
%!assert (fillmissing ([1 NaN 3], "constant", 0, "maxgap", 1.999), [1 NaN 3])
%!assert (fillmissing ([1 NaN 3], "constant", 0, "maxgap", 2), [1 0 3])
%!assert (fillmissing ([1 NaN NaN 4], "constant", 0, "maxgap", 2), [1 NaN NaN 4])
%!assert (fillmissing ([1 NaN NaN 4], "constant", 0, "maxgap", 3), [1 0 0 4])
%!assert (fillmissing ([1 NaN 3 NaN 5], "constant", 0, "maxgap", 2), [1 0 3 0 5])
%!assert (fillmissing ([NaN 2 NaN], "constant", 0, "maxgap", 0.999), [NaN 2 NaN])
%!assert (fillmissing ([NaN 2 NaN], "constant", 0, "maxgap", 1), [0 2 0])
%!assert (fillmissing ([NaN 2 NaN NaN], "constant", 0, "maxgap", 1), [0 2 NaN NaN])
%!assert (fillmissing ([NaN 2 NaN NaN], "constant", 0, "maxgap", 2), [0 2 0 0])
%!assert (fillmissing ([NaN NaN NaN], "constant", 0, "maxgap", 1), [NaN NaN NaN])
%!assert (fillmissing ([NaN NaN NaN], "constant", 0, "maxgap", 3), [NaN NaN NaN])
%!assert (fillmissing ([NaN NaN NaN], "constant", 0, "maxgap", 999), [NaN NaN NaN])
%!assert (fillmissing ([1 NaN 3 NaN 5], "constant", 0, "maxgap", 2, "samplepoints", [0 1 2 3 5]), [1 0 3 NaN 5])
%!assert (fillmissing ([1 NaN 3 NaN 5]', "constant", 0, "maxgap", 2, "samplepoints", [0 1 2 3 5]), [1 0 3 NaN 5]')
%!assert (fillmissing ([1 NaN 3 NaN 5], "constant", 0, "maxgap", 2, "samplepoints", [0 2 3 4 5]), [1 NaN 3 0 5])
%!assert (fillmissing ([1 NaN 3 NaN 5; 1 NaN 3 NaN 5], "constant", 0, 2, "maxgap", 2, "samplepoints", [0 2 3 4 5]), [1 NaN 3 0 5; 1 NaN 3 0 5])

%!test
%! x = cat (3, [1, 2, NaN; 4, NaN, NaN], [NaN, 2, 3; 4, 5, NaN]);
%! assert (fillmissing (x, "constant", 0, "maxgap", 0.1), x);
%! y = x;
%! y([4,7,12]) = 0;
%! assert (fillmissing (x, "constant", 0, "maxgap", 1), y);
%! assert (fillmissing (x, "constant", 0, 1, "maxgap", 1), y);
%! y = x;
%! y([5,7,12]) = 0;
%! assert (fillmissing (x, "constant", 0, 2, "maxgap", 1), y);
%! y = x;
%! y([4,5,7]) = 0;
%! assert (fillmissing (x, "constant", 0, 3, "maxgap", 1), y);

## 2nd output
## verify consistent with dim
%!test
%! x = cat (3, [1, 2, NaN; 4, NaN, NaN], [NaN, 2, 3; 4, 5, NaN]);
%! [~, idx] = fillmissing (x, "constant", 0, "maxgap", 1);
%! assert (idx, logical (cat (3, [0 0 0; 0 1 0], [1 0 0; 0 0 1])));
%! [~, idx] = fillmissing (x, "constant", 0, 1, "maxgap", 1);
%! assert (idx, logical (cat (3, [0 0 0; 0 1 0], [1 0 0; 0 0 1])));
%! [~, idx] = fillmissing (x, "constant", 0, 2, "maxgap", 1);
%! assert (idx, logical (cat (3, [0 0 1; 0 0 0], [1 0 0; 0 0 1])));
%! [~, idx] = fillmissing (x, "constant", 0, 3, "maxgap", 1);
%! assert (idx, logical (cat (3, [0 0 1; 0 1 0], [1 0 0; 0 0 0])));

## verify idx matches when methods leave gaps unfilled, or when fill looks
## the same
%!test
%! x = [NaN, 2, 3];
%! [~,idx] = fillmissing (x, "previous");
%! assert (idx, logical ([0 0 0]));
%! [~,idx] = fillmissing (x, "movmean", 1);
%! assert (idx, logical ([0 0 0]));
%! x = [1:3;4:6;7:9];
%! x([2,4,7,9]) = NaN;
%! [~,idx] = fillmissing (x, "linear");
%! assert (idx, logical ([0 1 0;1 0 0;0 0 0]));
%! [~,idx] = fillmissing (x, "movmean", 2);
%! assert (idx, logical ([0 0 0;1 0 0;0 0 1]));
%! [A, idx] = fillmissing ([1 2 3 NaN NaN], 'movmean',2);
%! assert (A, [1 2 3 3 NaN]);
%! assert (idx, logical([0 0 0 1 0]));
%! [A, idx] = fillmissing ([1 2 3 NaN NaN], 'movmean',3);
%! assert (A, [1 2 3 3 NaN]);
%! assert (idx, logical([0 0 0 1 0]));
%! [A, idx] = fillmissing ([1 2 NaN NaN NaN], 'movmedian', 2);
%! assert (A, [1 2 2 NaN NaN]);
%! assert (idx, logical([0 0 1 0 0]));
%! [A, idx] = fillmissing ([1 2 3 NaN NaN], 'movmedian', 3);
%! assert (A, [1 2 3 3 NaN]);
%! assert (idx, logical([0 0 0 1 0]));
%! [A, idx] = fillmissing ([1 NaN 1 NaN 1],  @(x,y,z) z, 3);
%! assert (A, [1 2 1 4 1]);
%! assert (idx, logical([0 1 0 1 0]));
%! [A, idx] = fillmissing ([1 NaN 1 NaN 1],  @(x,y,z) NaN (size (z)), 3);
%! assert (A, [1 NaN 1 NaN 1]);
%! assert (idx, logical([0 0 0 0 0]));

#test missinglocations
%!assert (fillmissing ([1 2 3], "constant", 99, "missinglocations", logical([0 0 0])), [1 2 3])
%!assert (fillmissing ([1 2 3], "constant", 99, "missinglocations", logical([1 1 1])), [99 99 99])
%!assert (fillmissing ([1 NaN 2 3 NaN], "constant", 99, "missinglocations", logical([1 0 1 0 1])), [99 NaN 99 3 99])
%!assert (fillmissing ([1 NaN 3 NaN 5], "constant", NaN, "missinglocations", logical([0 1 1 1 0])), [1 NaN NaN NaN 5])
%!assert (fillmissing (["foo ";" bar"], "constant", 'X', "missinglocations", logical([0 0 0 0; 0 0 0 0])), ["foo ";" bar"])
%!assert (fillmissing (["foo ";" bar"], "constant", 'X', "missinglocations", logical([1 0 1 0; 0 1 1 0])), ["XoX ";" XXr"])
%!assert (fillmissing ({"foo","", "bar"}, "constant", 'X', "missinglocations", logical([0 0 0])), {"foo","", "bar"})
%!assert (fillmissing ({"foo","", "bar"}, "constant", 'X', "missinglocations", logical([1 1 0])), {"X","X","bar"})
%!test
%! [~,idx] = fillmissing ([1 NaN 3 NaN 5], "constant", NaN);
%! assert (idx, logical([0 0 0 0 0]));
%! [~,idx] = fillmissing ([1 NaN 3 NaN 5], "constant", NaN, "missinglocations", logical([0 1 1 1 0]));
%! assert (idx, logical([0 1 1 1 0]));
%! [A, idx] = fillmissing ([1 2 NaN 1 NaN], 'movmean', 3.1, 'missinglocations', logical([0 0 1 1 0]));
%! assert (A, [1 2 2 NaN NaN]);
%! assert (idx, logical([0 0 1 0 0]));
%! [A, idx] = fillmissing ([1 2 NaN NaN NaN], 'movmean', 2, 'missinglocations', logical([0 0 1 1 0]));
%! assert (A, [1 2 2 NaN NaN]);
%! assert (idx, logical([0 0 1 0 0]));
%! [A, idx] = fillmissing ([1 2 NaN 1 NaN], 'movmean', 3, 'missinglocations', logical([0 0 1 1 0]));
%! assert (A, [1 2 2 NaN NaN]);
%! assert (idx, logical([0 0 1 0 0]));
%! [A, idx] = fillmissing ([1 2 NaN NaN NaN], 'movmean', 3, 'missinglocations', logical([0 0 1 1 0]));
%! assert (A, [1 2 2 NaN NaN]);
%! assert (idx, logical([0 0 1 0 0]));
%! [A, idx] = fillmissing ([1 2 NaN NaN NaN], 'movmedian', 2, 'missinglocations', logical([0 0 1 1 0]));
%! assert (A, [1 2 2 NaN NaN]);
%! assert (idx, logical([0 0 1 0 0]));
%! [A, idx] = fillmissing ([1 2 NaN NaN NaN], 'movmedian', 3, 'missinglocations', logical([0 0 1 1 0]));
%! assert (A, [1 2 2 NaN NaN]);
%! assert (idx, logical([0 0 1 0 0]));
%! [A, idx] = fillmissing ([1 2 NaN NaN NaN], 'movmedian', 3.1, 'missinglocations', logical([0 0 1 1 0]));
%! assert (A, [1 2 2 NaN NaN]);
%! assert (idx, logical([0 0 1 0 0]));
%! [A, idx] = fillmissing ([1 NaN 1 NaN 1],  @(x,y,z) ones (size (z)), 3, "missinglocations", logical([0 1 0 1 1]));
%! assert (A, [1 1 1 1 1]);
%! assert (idx, logical([0 1 0 1 1]));
%! [A, idx] = fillmissing ([1 NaN 1 NaN 1],  @(x,y,z) NaN (size (z)), 3, "missinglocations", logical([0 1 0 1 1]));
%! assert (A, [1 NaN 1 NaN NaN]);
%! assert (idx, logical([0 0 0 0 0]));


##Test char and cellstr
%!assert (fillmissing (' foo bar ', "constant", 'X'), 'XfooXbarX')
%!assert (fillmissing ([' foo';'bar '], "constant", 'X'), ['Xfoo';'barX'])
%!assert (fillmissing ([' foo';'bar '], "next"), ['bfoo';'bar '])
%!assert (fillmissing ([' foo';'bar '], "next", 1), ['bfoo';'bar '])
%!assert (fillmissing ([' foo';'bar '], "previous"), [' foo';'baro'])
%!assert (fillmissing ([' foo';'bar '], "previous", 1), [' foo';'baro'])
%!assert (fillmissing ([' foo';'bar '], "nearest"), ['bfoo';'baro'])
%!assert (fillmissing ([' foo';'bar '], "nearest", 1), ['bfoo';'baro'])
%!assert (fillmissing ([' foo';'bar '], "next", 2), ['ffoo';'bar '])
%!assert (fillmissing ([' foo';'bar '], "previous", 2), [' foo';'barr'])
%!assert (fillmissing ([' foo';'bar '], "nearest", 2), ['ffoo';'barr'])
%!assert (fillmissing ([' foo';'bar '], "next", 3), [' foo';'bar '])
%!assert (fillmissing ([' foo';'bar '], "previous", 3), [' foo';'bar '])
%!assert (fillmissing ([' foo';'bar '], "nearest", 3), [' foo';'bar '])
%!assert (fillmissing ({'foo','bar'}, "constant", 'a'), {'foo','bar'})
%!assert (fillmissing ({'foo','bar'}, "constant", {'a'}), {'foo','bar'})
%!assert (fillmissing ({'foo', '', 'bar'}, "constant", 'a'), {'foo', 'a', 'bar'})
%!assert (fillmissing ({'foo', '', 'bar'}, "constant", {'a'}), {'foo', 'a', 'bar'})
%!assert (fillmissing ({'foo', '', 'bar'}, "previous"), {'foo', 'foo', 'bar'})
%!assert (fillmissing ({'foo', '', 'bar'}, "next"), {'foo', 'bar', 'bar'})
%!assert (fillmissing ({'foo', '', 'bar'}, "nearest"), {'foo', 'bar', 'bar'})
%!assert (fillmissing ({'foo', '', 'bar'}, "previous", 2), {'foo', 'foo', 'bar'})
%!assert (fillmissing ({'foo', '', 'bar'}, "next", 2), {'foo', 'bar', 'bar'})
%!assert (fillmissing ({'foo', '', 'bar'}, "nearest", 2), {'foo', 'bar', 'bar'})
%!assert (fillmissing ({'foo', '', 'bar'}, "previous", 1), {'foo', '', 'bar'})
%!assert (fillmissing ({'foo', '', 'bar'}, "previous", 1), {'foo', '', 'bar'})
%!assert (fillmissing ({'foo', '', 'bar'}, "next", 1), {'foo', '', 'bar'})
%!assert (fillmissing ({'foo', '', 'bar'}, "nearest", 1), {'foo', '', 'bar'})
%!assert (fillmissing ("abc ", @(x,y,z) x+y+z, 2), "abcj")
%!assert (fillmissing ({'foo', '', 'bar'}, @(x,y,z) x(1), 3), {'foo','foo','bar'})

%!test
%! [A, idx] = fillmissing (" a b c", "constant", " ");
%! assert (A, " a b c");
%! assert (idx, logical([0 0 0 0 0 0]));
%! [A, idx] = fillmissing ({"foo", "", "bar", ""}, "constant", "");
%! assert (A, {"foo", "", "bar", ""});
%! assert (idx, logical([0 0 0 0]));
%! [A, idx] = fillmissing ({"foo", "", "bar", ""}, "constant", {""});
%! assert (A, {"foo", "", "bar", ""});
%! assert (idx, logical([0 0 0 0]));
%! [A,idx] = fillmissing (' f o o ', @(x,y,z) repelem ("a", numel (z)), 3);
%! assert (A, "afaoaoa");
%! assert (idx, logical([1 0 1 0 1 0 1]));
%! [A,idx] = fillmissing (' f o o ', @(x,y,z) repelem (" ", numel (z)), 3);
%! assert (A, " f o o ");
%! assert (idx, logical([0 0 0 0 0 0 0]));
%! [A,idx] = fillmissing ({'','foo',''}, @(x,y,z) repelem ({'a'}, numel (z)), 3);
%! assert (A, {'a','foo','a'});
%! assert (idx, logical([1 0 1]));
%! [A,idx] = fillmissing ({'','foo',''}, @(x,y,z) repelem ({''}, numel (z)), 3);
%! assert (A, {'','foo',''});
%! assert (idx, logical([0 0 0]));


##types without a defined 'missing' (currently logical, int) that can be filled
%!assert (fillmissing (logical ([1 0 1 0 1]), "constant", true), logical ([1 0 1 0 1]))
%!assert (fillmissing (logical ([1 0 1 0 1]), "constant", false, 'missinglocations', logical([1 0 1 0 1])), logical ([0 0 0 0 0]))
%!assert (fillmissing (logical ([1 0 1 0 1]), "previous",  'missinglocations', logical([1 0 1 0 1])), logical ([1 0 0 0 0]))
%!assert (fillmissing (logical ([1 0 1 0 1]), "next",  'missinglocations', logical([1 0 1 0 1])), logical ([0 0 0 0 1]))
%!assert (fillmissing (logical ([1 0 1 0 1]), "nearest", 'missinglocations', logical([1 0 1 0 1])), logical ([0 0 0 0 0]))
%!assert (fillmissing (logical ([1 0 1 0 1]),  @(x,y,z) false(size(z)), 3), logical ([1 0 1 0 1]))
%!assert (fillmissing (logical ([1 0 1 0 1]),  @(x,y,z) false(size(z)), 3, 'missinglocations', logical([1 0 1 0 1])), logical ([0 0 0 0 0]))
%!assert (fillmissing (logical ([1 0 1 0 1]),  @(x,y,z) false(size(z)), [2 0], 'missinglocations', logical([1 0 1 0 1])), logical ([1 0 0 0 0]))
%!test
%! x = logical ([1 0 1 0 1]);
%! [~,idx] = fillmissing (x, "constant", true);
%! assert (idx, logical([0 0 0 0 0]));
%! [~,idx] = fillmissing (x, "constant", false, 'missinglocations', logical([1 0 1 0 1]));
%! assert (idx, logical([1 0 1 0 1]));
%! [~,idx] = fillmissing (x, "constant", true, 'missinglocations', logical([1 0 1 0 1]));
%! assert (idx, logical([1 0 1 0 1]));
%! [~,idx] = fillmissing (x, "previous", 'missinglocations', logical([1 0 1 0 1]));
%! assert (idx, logical([0 0 1 0 1]));
%! [~,idx] = fillmissing (x, "next",  'missinglocations', logical([1 0 1 0 1]));
%! assert (idx, logical([1 0 1 0 0]));
%! [~,idx] = fillmissing (x, "nearest", 'missinglocations', logical([1 0 1 0 1]));
%! assert (idx, logical([1 0 1 0 1]));
%! [~,idx] = fillmissing (x, @(x,y,z) false(size(z)), 3);
%! assert (idx, logical ([0 0 0 0 0]))
%! [~,idx] = fillmissing (x, @(x,y,z) false(size(z)), 3, 'missinglocations', logical([1 0 1 0 1]));
%! assert (idx, logical ([1 0 1 0 1]))
%! [~,idx] = fillmissing (x, @(x,y,z) false(size(z)), [2 0], 'missinglocations', logical([1 0 1 0 1]));
%! assert (idx, logical ([0 0 1 0 1]))

%!assert (fillmissing (int32 ([1 2 3 4 5]), "constant", 0), int32 ([1 2 3 4 5]))
%!assert (fillmissing (int32 ([1 2 3 4 5]), "constant", 0,'missinglocations', logical([1 0 1 0 1])), int32 ([0 2 0 4 0]))
%!assert (fillmissing (int32 ([1 2 3 4 5]), "previous", 'missinglocations', logical([1 0 1 0 1])), int32 ([1 2 2 4 4]))
%!assert (fillmissing (int32 ([1 2 3 4 5]), "next", 'missinglocations', logical([1 0 1 0 1])), int32 ([2 2 4 4 5]))
%!assert (fillmissing (int32 ([1 2 3 4 5]), "nearest", 'missinglocations', logical([1 0 1 0 1])), int32 ([2 2 4 4 4]))
%!assert (fillmissing (int32 ([1 2 3 4 5]), @(x,y,z) z+10, 3), int32 ([1 2 3 4 5]))
%!assert (fillmissing (int32 ([1 2 3 4 5]), @(x,y,z) z+10, 3, 'missinglocations', logical([1 0 1 0 1])), int32 ([11 2 13 4 15]))
%!assert (fillmissing (int32 ([1 2 3 4 5]), @(x,y,z) z+10, [2 0], 'missinglocations', logical([1 0 1 0 1])), int32 ([1 2 13 4 15]))
%!test
%! x = int32 ([1 2 3 4 5]);
%! [~,idx] = fillmissing (x, "constant", 0);
%! assert (idx, logical([0 0 0 0 0]));
%! [~,idx] = fillmissing (x, "constant", 0, 'missinglocations', logical([1 0 1 0 1]));
%! assert (idx, logical([1 0 1 0 1]));
%! [~,idx] = fillmissing (x, "constant", 3, 'missinglocations', logical([0 0 1 0 0]));
%! assert (idx, logical([0 0 1 0 0]));
%! [~,idx] = fillmissing (x, "previous", 'missinglocations', logical([1 0 1 0 1]));
%! assert (idx, logical([0 0 1 0 1]));
%! [~,idx] = fillmissing (x, "next", 'missinglocations', logical([1 0 1 0 1]));
%! assert (idx, logical([1 0 1 0 0]));
%! [~,idx] = fillmissing (x, "nearest", 'missinglocations', logical([1 0 1 0 1]));
%! assert (idx, logical([1 0 1 0 1]));
%! [~,idx] = fillmissing (x, @(x,y,z) z+10, 3);
%! assert (idx, logical([0 0 0 0 0]));
%! [~,idx] = fillmissing (x, @(x,y,z) z+10, 3, 'missinglocations', logical([1 0 1 0 1]));
%! assert (idx, logical([1 0 1 0 1]));
%! [~,idx] = fillmissing (x, @(x,y,z) z+10, [2 0], 'missinglocations', logical([1 0 1 0 1]));
%! assert (idx, logical([0 0 1 0 1]));

## other data type passthrough
%!test
%! [A, idx] = fillmissing ([struct struct], "constant", 1);
%! assert (A, [struct struct])
%! assert (idx, [false false])

## Test input validation and error messages
%!error <Invalid call> fillmissing ()
%!error <Invalid call> fillmissing (1)
%!error <Invalid call> fillmissing (1,2,3,4,5,6,7,8,9,10,11,12,13)
%!error <second input must be a> fillmissing (1, 2)
%!error <unknown fill method 'foo'> fillmissing (1, "foo")
%!error <fill function must accept at least> fillmissing (1, @(x) x, 1)
%!error <fill function must accept at least> fillmissing (1, @(x,y) x+y, 1)
%!error <interpolation methods only valid for numeric> fillmissing ("a b c", "linear")
%!error <interpolation methods only valid for numeric> fillmissing ({'a','b'}, "linear")
%!error <'movmean' and 'movmedian' methods only valid for numeric> fillmissing ("a b c", "movmean", 2)
%!error <'movmean' and 'movmedian' methods only valid for numeric> fillmissing ({'a','b'}, "movmean", 2)
%!error <'constant' method must be followed by> fillmissing (1, "constant")
%!error <a numeric fill value cannot be emtpy> fillmissing (1, "constant", [])
%!error <fill value must be the same data type> fillmissing (1, "constant", "a")
%!error <fill value must be the same data type> fillmissing ("a", "constant", 1)
%!error <fill value must be the same data type> fillmissing ("a", "constant", {"foo"})
%!error <fill value must be the same data type> fillmissing ({"foo"}, "constant", 1)
%!error <moving window method must be followed by> fillmissing (1, "movmean")
%!error <moving window method must be followed by> fillmissing (1, "movmedian")
%!error <DIM must be a positive scalar> fillmissing (1, "constant", 1, 0)
%!error <DIM must be a positive scalar> fillmissing (1, "constant", 1, -1)
%!error <DIM must be a positive scalar> fillmissing (1, "constant", 1, [1 2])
%!error <properties must be given as> fillmissing (1, "constant", 1, "samplepoints")
%!error <properties must be given as> fillmissing (1, "constant", 1, "foo")
%!error <properties must be given as> fillmissing (1, "constant", 1, 1, "foo")
%!error <invalid parameter name specified> fillmissing (1, "constant", 1, 2, {1}, 4)
%!error <SamplePoints must be a> fillmissing ([1 2 3], "constant", 1, 2, "samplepoints", [1 2])
%!error <SamplePoints must be a> fillmissing ([1 2 3], "constant", 1, 2, "samplepoints", [3 1 2])
%!error <SamplePoints must be a> fillmissing ([1 2 3], "constant", 1, 2, "samplepoints", [1 1 2])
%!error <SamplePoints must be a> fillmissing ([1 2 3], "constant", 1, 2, "samplepoints", "abc")
%!error <SamplePoints must be a> fillmissing ([1 2 3], "constant", 1, 2, "samplepoints", logical([1 1 1]))
%!error <SamplePoints must be a> fillmissing ([1 2 3], "constant", 1, 1, "samplepoints", [1 2 3])
%!error <EndValues method 'constant' only valid> fillmissing ('foo', "next", "endvalues", 1)
%!error <invalid EndValues method 'foo'> fillmissing (1, "constant", 1, 1, "endvalues", "foo")
%!error <EndValues must be a scalar or a 1 element array> fillmissing ([1 2 3], "constant", 1, 2, "endvalues", [1 2 3])
%!error <EndValues must be a scalar or a 3 element array> fillmissing ([1 2 3], "constant", 1, 1, "endvalues", [1 2])
%!error <EndValues must be a scalar or a 12 element array> fillmissing (randi(5,4,3,2), "constant", 1, 3, "endvalues", [1 2])
%!error <EndValues must be numeric or a> fillmissing (1, "constant", 1, 1, "endvalues", {1})
%!error <invalid parameter name 'foo'> fillmissing (1, "constant", 1, 2, "foo", 4)
%!error <MissingLocations option is not compatible with> fillmissing (struct, "constant", 1, "missinglocations", false)
%!error <MissingLocations and MaxGap options> fillmissing (1, "constant", 1, 2, "maxgap", 1, "missinglocations", false)
%!error <MissingLocations and MaxGap options> fillmissing (1, "constant", 1, 2, "missinglocations", false, "maxgap", 1)
%!error <the 'replacevalues' option has not> fillmissing (1, "constant", 1, "replacevalues", true)
%!error <the 'datavariables' option has not> fillmissing (1, "constant", 1, "datavariables", 'Varname')
%!error <MissingLocations must be a> fillmissing (1, "constant", 1, 2, "missinglocations", 1)
%!error <MissingLocations must be a> fillmissing (1, "constant", 1, 2, "missinglocations", 'a')
%!error <MissingLocations must be a> fillmissing (1, "constant", 1, 2, "missinglocations", [true false])
%!error <MissingLocations cannot be used with method> fillmissing (true, "linear", "missinglocations", true)
%!error <MissingLocations cannot be used with method> fillmissing (int8(1), "linear", "missinglocations", true)
%!error <MissingLocations cannot be used with EndValues method> fillmissing (true, "next", "missinglocations", true, "EndValues", "linear")
%!error <MissingLocations cannot be used with EndValues method> fillmissing (true, "next", "EndValues", "linear", "missinglocations", true)
%!error <MissingLocations cannot be used with EndValues method> fillmissing (int8(1), "next", "missinglocations", true, "EndValues", "linear")
%!error <MissingLocations cannot be used with EndValues method> fillmissing (int8(1), "next", "EndValues", "linear", "missinglocations", true)
%!error <MaxGap must be a positive numeric scalar> fillmissing (1, "constant", 1, 2, "maxgap", true)
%!error <MaxGap must be a positive numeric scalar> fillmissing (1, "constant", 1, 2, "maxgap", 'a')
%!error <MaxGap must be a positive numeric scalar> fillmissing (1, "constant", 1, 2, "maxgap", [1 2])
%!error <MaxGap must be a positive numeric scalar> fillmissing (1, "constant", 1, 2, "maxgap", 0)
%!error <MaxGap must be a positive numeric scalar> fillmissing (1, "constant", 1, 2, "maxgap", -1)
%!error <fill value 'V' must be a scalar or a 1> fillmissing ([1 2 3], "constant", [1 2 3])
%!error <fill value 'V' must be a scalar or a 1> fillmissing ([1 2 3]', "constant", [1 2 3])
%!error <fill value 'V' must be a scalar or a 1> fillmissing ([1 2 3]', "constant", [1 2 3], 1)
%!error <fill value 'V' must be a scalar or a 1> fillmissing ([1 2 3], "constant", [1 2 3], 2)
%!error <fill value 'V' must be a scalar or a 6> fillmissing (randi(5,4,3,2), "constant", [1 2], 1)
%!error <fill value 'V' must be a scalar or a 8> fillmissing (randi(5,4,3,2), "constant", [1 2], 2)
%!error <fill value 'V' must be a scalar or a 12> fillmissing (randi(5,4,3,2), "constant", [1 2], 3)
%!error <fill function handle must be followed by> fillmissing (1, @(x,y,z) x+y+z)
%!error <fill function return values must be the same size> fillmissing ([1 NaN 2], @(x,y,z) [1 2], 2)
