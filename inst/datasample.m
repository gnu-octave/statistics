## Copyright (C) 2021 Stefano Guidoni <ilguido@users.sf.net>
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
## @deftypefn  {statistics} @var{y} = datasample (@var{data}, @var{k})
## @deftypefnx {statistics} @var{y} = datasample (@var{data}, @var{k}, @var{dim})
## @deftypefnx {statistics} @var{y} = datasample (@dots{}, @var{Name}, @var{Value})
## @deftypefnx {statistics} [@var{y} @var{idcs}] = datasample (@dots{})
##
## Randomly sample data.
##
## Return @var{k} observations randomly sampled from @var{data}.  @var{data} can
## be a vector or a matrix of any data.  When @var{data} is a matrix or a
## n-dimensional array, the samples are the subarrays of size n - 1, taken along
## the dimension @var{dim}. The default value for @var{dim} is 1, that is the
## row vectors when sampling a matrix.
##
## Output @var{y} is the returned sampled data. Optional output @var{idcs} is
## the vector of the indices to build @var{y} from @var{data}.
##
## Additional options are set through pairs of parameter name and value.
## Available parameters are:
##
## @table @code
## @item @qcode{Replace}
## a logical value that can be @code{true} (default) or @code{false}: when set
## to @code{true}, @code{datasample} returns data sampled with replacement.
##
## @item @qcode{Weigths}
## a vector of positive numbers that sets the probability of each element.  It
## must have the same size as @var{data} along dimension @var{dim}.
##
## @end table
##
##
## @end deftypefn
##
## @seealso{rand, randi, randperm, randsample}

function [y, idcs] = datasample (data, k, varargin)

  ## check input
  if ( nargin < 2 )
    print_usage ();
  endif

  ## data: some data, any type, any format but cell
  ## MATLAB compatibility: there are no "table" or "dataset array" types in
  ## Octave
  if (iscell (data))
    error ("datasample: data must be a vector or matrix");
  endif

  ## k, a positive integer
  if ((! isnumeric (k) || ! isscalar (k)) || (! (floor (k) == k)) || (k <= 0))
    error ("datasample: k must be a positive integer scalar");
  endif

  dim = 1;
  replace = true;
  weights = [];
  if ( nargin > 2 )
    pair_index = 1;

    if (! ischar (varargin{1}))
      ## it must be dim
      dim = varargin{1};

      ## the (Name, Value) pairs start further
      pair_index += 1;

      ## dim, another positive integer
      if ((! isscalar (dim)) || (! (floor (dim) == dim)) || (dim <= 0))
        error ("datasample: DIM must be a positive integer scalar");
      endif
    endif

    ## (Name, Value) pairs
    while (pair_index < (nargin - 2))
      switch (lower (varargin{pair_index}))
        case "replace"
          if (! islogical (varargin{pair_index + 1}))
            error ("datasample: expected a logical value for 'Replace'");
          endif
          replace = varargin{pair_index + 1};
        case "weights"
          if ((! isnumeric (varargin{pair_index + 1})) ||
              (! isvector (varargin{pair_index + 1})) ||
              (any (varargin{pair_index + 1} < 0)))
            error (["datasample: the sampling weights must be defined as a " ...
              "vector of positive values"]);
          endif
          weights = varargin{pair_index + 1};
        otherwise
          error ("datasample: unknown property %s", varargin{pair_index});
      endswitch
      pair_index += 2;
    endwhile
  endif

  ## get the size of the population to sample
  if (isvector (data))
    imax = length (data);
  else
    imax = size (data, dim);
  endif

  if (isempty (weights))
    ## all elements have the same probability of being chosen
    ## this is easy

    ## with or without replacement
    if (replace)
      idcs = randi (imax, k, 1);
    else
      idcs = randperm (imax, k);
    endif
  else
    ## first check if the weights vector is right
    if (imax != length (weights))
      error (["datasample: the size of the vector of sampling weights must"...
        " be equal to the size of the sampled data"]);
    endif

    if (replace)
      ## easy case:
      ## normalize the weights,
      weights_n = cumsum (weights ./ sum (weights));
      weights_n(end) = 1; # just to be sure
      ## then choose k numbers uniformly between 0 and 1
      samples = rand (k, 1);

      ## we have subdivided the space between 0 and 1 accordingly to the
      ## weights vector: we have just to map back the random numbers to the
      ## indices of the orginal dataset
      for iter = 1 : k
        idcs(iter) = find (weights_n >= samples(iter), 1);
      endfor
    else
      ## complex case
      ## choose k numbers uniformly between 0 and 1
      samples = rand (k, 1);

      for iter = 1 : k
        ## normalize the weights
        weights_n = cumsum (weights ./ sum (weights));
        weights_n(end) = 1; # just to be sure

        idcs(iter) = find (weights_n >= samples(iter), 1);

        ## remove the element from the set, i. e. set its probability to zero
        weights(idcs(iter)) = 0;
      endfor
    endif
  endif

  ## let's get the actual data from the original set
  if (isvector (data))
    ## data is a vector
    y = data(idcs);
  else
    vS = size (data);

    if (length (vS) == 2)
      ## data is a 2-dimensional matrix
      if (dim == 1)
        y = data(idcs, :);
      else
        y = data(:, idcs);
      endif
    else
      ## data is an n-dimensional matrix
      s = "y = data(";
      for iter = 1 : length (vS)
        if (iter == dim)
          s = [s "idcs,"];
        else
          s = [s ":,"];
        endif
      endfor
      s = [s ":);"];
      eval (s);
    endif
  endif

endfunction

## some tests
%!error datasample();
%!error datasample(1);
%!error <data must be a vector or matrix> datasample({1, 2, 3}, 1);
%!error <k must be a positive integer scalar> datasample([1 2], -1);
%!error <k must be a positive integer scalar> datasample([1 2], 1.5);
%!error <k must be a positive integer scalar> datasample([1 2], [1 1]);
%!error <k must be a positive integer scalar> datasample([1 2], 'g', [1 1]);
%!error <DIM must be a positive integer scalar> datasample([1 2], 1, -1);
%!error <DIM must be a positive integer scalar> datasample([1 2], 1, 1.5);
%!error <DIM must be a positive integer scalar> datasample([1 2], 1, [1 1]);
%!error <Replace> datasample([1 2], 1, 1, "Replace", -2);
%!error <weights must be defined> datasample([1 2], 1, 1, "Weights", "abc");
%!error <weights must be defined> datasample([1 2], 1, 1, "Weights", [1 -2 3]);
%!error <weights must be defined> datasample([1 2], 1, 1, "Weights", ones (2));
%!error <weights must be equal> datasample([1 2], 1, 1, "Weights", [1 2 3]);

%!test
%! dat = randn (10, 4);
%! assert (size (datasample (dat, 3, 1)), [3 4]);

%!test
%! dat = randn (10, 4);
%! assert (size (datasample (dat, 3, 2)), [10 3]);

