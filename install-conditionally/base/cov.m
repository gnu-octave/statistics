## Copyright (C) 1995-2017 Kurt Hornik
## Copyright (C) 2017 Nicholas R. Jankowski
##
## This program is free software: you can redistribute it and/or
## modify it under the terms of the GNU General Public License as
## published by the Free Software Foundation, either version 3 of the
## License, or (at your option) any later version.
##
## This program is distributed in the hope that it will be useful, but
## WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
## General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program; see the file COPYING.  If not, see
## <http://www.gnu.org/licenses/>.

## -*- texinfo -*-
## @deftypefn  {} {} cov (@var{x})
## @deftypefnx {} {} cov (@var{x}, @var{opt})
## @deftypefnx {} {} cov (@var{x}, @var{y})
## @deftypefnx {} {} cov (@var{x}, @var{y}, @var{opt})
## @deftypefnx {} {} cov (@var{x}, @var{y}, @var{opt}, @var{NaN-option})
## Compute the covariance matrix.
##
## If each row of @var{x} and @var{y} is an observation, and each column is
## a variable, then the @w{(@var{i}, @var{j})-th} entry of
## @code{cov (@var{x}, @var{y})} is the covariance between the @var{i}-th
## variable in @var{x} and the @var{j}-th variable in @var{y}.
## @tex
## $$
## \sigma_{ij} = {1 \over N-1} \sum_{i=1}^N (x_i - \bar{x})(y_i - \bar{y})
## $$
## where $\bar{x}$ and $\bar{y}$ are the mean values of @var{x} and @var{y}.
## @end tex
## @ifnottex
##
## @example
## cov (@var{x}) = 1/(N-1) * SUM_i (@var{x}(i) - mean(@var{x})) * (@var{y}(i) - mean(@var{y}))
## @end example
##
## where @math{N} is the length of the @var{x} and @var{y} vectors.
##
## @end ifnottex
##
## If called with one argument, compute @code{cov (@var{x}, @var{x})}, the
## covariance between the columns of @var{x}.
##
## If called with two vector arguments, compute 
## @code{cov (@var{x}, @var{y})}, the covariance between two random variables
## @var{x} and @var{y}. The output will be the 2 by 2 covariance matrix. 
##
## If called with two matrix arguments, the matrices are treated as vectors and 
## covariance is computed as @code{cov (@var{x}(:), @var{y}(:))}. The output 
## will be the 2 by 2 covariance matrix. 
##
## The optional argument @var{opt} determines the type of normalization to use.
## Valid values are
##
## @table @asis
## @item 0:
##   normalize with @math{N-1}, provides the best unbiased estimator of the
## covariance [default]
##
## @item 1:
##   normalize with @math{N}, this provides the second moment around the mean
## @end table
##
## The optional argument @var{NaN-option} controls how @code{cov} deals with NaN
## values in the data.  The three valid values are
##
## @table @asis
## @item includenan:
##   leave NaN values in @var{x} and @var{y}. Output will follow the normal 
##   rules for handling NaN values in arithemtic operations [default]
##
## @item omitnans:
##   rows containing NaN values are trimmed from both @var{x} and @var{y} prior
##   to calculating the covariance. (A NaN in one variable will that row from
##   both @var{x} and @var{y}.)
##
## @item partialnans:
##   rows containing NaN values are ignored from both @var{x} and @var{y} 
##   independently when for each @var{i}-th and @var{j}-th covariance 
##   calculation. This may result in a different number of observations,
##   @math{N}, being used to calculated each element of the covariance matrix.
## @end table
##
## Compatibility Note: This version of @code{cov} attempts to maintain full 
##   compatibility with @sc{matlab}'s cov function by treating @var{x} and 
##   @var{y} as two univariate distributions regardless of shape, resulting in
##   a 2x2 output matrix.  Previous versions of cov in Octave treated rows
##   of @var{x} and @var{y} as multivariate random variables. Code relying on 
##   Octave's previous definition will need to be changed when running this newer
##   version of @code{cov}.
## @seealso{corr}
## @end deftypefn

## Author: KH <Kurt.Hornik@wu-wien.ac.at>
## Author: Nicholas Jankowski <jankowskin@asme.org>
## Description: Compute covariances

function c = cov (x, varargin)

  %%input sorting

  switch nargin
    case 1
      [y, opt, handlenan] = deal ({"no_y", 0, "includenan"}{:});
      
    case 4
      [y, opt, handlenan] = deal (varargin{:});      
       
    case {2,3}
      [y, opt, handlenan] = deal ({"no_y", 0, "includenan"}{:});

      for vararg_idx = 1 : (nargin-1)
        v = varargin{vararg_idx};
        if ischar (v)
          if (vararg_idx == 1 && nargin == 3)
            error ('cov: NaN handling string must be the last input');
          else
            handlenan = v;
         
          endif
          
        else
          if (isscalar(v) && (v == 0 || v == 1))
            opt = v;
          
          elseif (vararg_idx ~= 2)
            y = v;
          
          else
            print_usage();
            
          endif
        endif
        
      endfor
    
    otherwise
      print_usage();
  
  endswitch

  ## check sorted X
  if ~((isnumeric (x) || islogical (x)) && (ndims (x) == 2))
    error ("cov: X must be a numeric 2-D matrix or vector");
  endif

  ##vector x needs to be column for calulations, flip before any nantrimming
  if (isrow (x))
    x = x';
  endif

  if ~(strcmp (y, "no_y"))
  ## check sorted Y assuming one is given
    if ~((isnumeric (y) || islogical (y)) && (ndims (y) == 2))
      error ("cov: Y must be a numeric 2-D matrix or vector");
    endif

    if (numel (x) ~= numel (y))
      error ("cov: X and Y must have the same number of elements");
    endif
  endif

  ## check sorted opt
  if (opt ~= 0 && opt ~= 1)
    error ("cov: normalization factor OPT must be either 0 or 1");
  endif
    
  ## check sorted NaN handling switch, correct for spelling, adjust x and y
  switch handlenan
    case {"includenan"}
      ## okay, do nothing
    case {"omitrows", "omitrow"}
      handlenan = "omitrows";
      if (strcmp (y, "no_y"))
        #trim out rows with nans from x
        x = x(~any (isnan (x), 2), :);
      else
        nan_locs = any (isnan ([x(:), y(:)]), 2);
        x = x(~nan_locs);
        y = y(~nan_locs);
               
      endif
      
    case {"partialrows", "partialrow"}
      handlenan = "partialrows";
      if ~(strcmp (y, "no_y"))
        ##no need to handle anything differently for single input
        x_nan_locs = any (isnan (x(:)), 2);
        y_nan_locs = any (isnan (y(:)), 2);
        both_nan_locs = any (isnan ([x(:), y(:)]), 2);
        
        x_xytrim  = x(~both_nan_locs);
        y_xytrim  = y(~both_nan_locs);
        
        x = x(~x_nan_locs);
        y = y(~y_nan_locs);
        
      endif
        
    otherwise
      error (["cov: unknown NaN handling parameter, '", handlenan, "'"]);
  endswitch   

  ## opt being single shouldn't affect output
  if (isa (opt, "single"))
    opt = double (opt);
  endif
    
  ## end input sorting/checking
 
  ## Primary handling difference is whether there are one or two inputs.
  if (strcmp (y, "no_y"))
    
    ## Special case, scalar has zero covariance
    if (isscalar (x)) 
      if isnan (x)
        c = NaN;
      else
        c = 0;
      endif
      
      if (isa (x, "single"))
        c = single (c);
      endif

      return;
    
    elseif (isempty (x)) %not scalar x, check if empty
      sx = size (x);
     
      if all (sx == 0)
        c = NaN;
      elseif (sx(1) > 0)
        c = [];
      else
        c = NaN (sx(2));
      endif
      
      if (isa (x, "single"))
        c = single (c);
      endif
      
      return;
      
    else %not scalar x, not empty, no y, generate covariance matrix
   
      if strcmp (handlenan, 'partialrows')
        ## if 'partialrows', need to calc each element separately with 'omitrows'
        ##  c(i,j) is cov(x(:,i),x(:,j),'omitrows)
        ## TODO:  find more efficient method. maybe can flatten recursion
        
        szx = size (x);

        for rw = 1:szx(1)
          for cl = 1:szx(2)
            c(rw,cl) = (cov (x(:,rw), x(:,cl), 'omitrows'))(2);
          endfor
        endfor
        return
      
      else
        ## if some elements are NaN, they propagate through calc as needed
        
        n = rows (x);
        x = center (x, 1);
        
        if n == 1
          c = x' * x; %% to handle case of omitrows trimming to 1 row
        else
          c = x' * x / (n - 1 + opt); %will preserve single type
        endif
        
        return;   
      endif
    
    endif

  else %there is a y
    
    if (isscalar (x))
      if (isnan (x) || isnan (y))
        c = NaN (2, 2);
        if ~isnan (x)
          c(1,1) = 0;
        elseif ~isnan (y)
          c(2,2) = 0;
        endif
        
        if (isa (x, "single") || isa (y, "single") )
          c = single (c);
        endif
        
        return;
      
      else  %scalar but neither a nan... both should be numbers... 
      
        if (isa (x, "single") || isa (y, "single") )
          c = single ([0 0; 0 0]);
        else
          c = [0 0; 0 0];
        endif

        return;
                
      endif  
    
    else % matrix or vector handled the same way, generate 2x2 covariance matrix

      if (isempty (x) || isempty (y))
         if (isa (x, "single") || isa (y, "single") )
          c = single (NaN (2, 2));
        else
          c = NaN (2, 2);
        endif
        
        return;
      endif   

      if (~strcmp (handlenan, 'partialrows'))
      
        denom = numel(x) - 1 + opt;
        x = center (x(:), 1);
        y = center (y(:), 1);
        
        c1 = x' * x;
        c23 = x' * y;
        c4 = y' * y;
        
        c = [c1, c23; c23, c4] ./ denom;
        
        return;
        
      else
        ## 'partialrows': handle each element separatley
        
        denom_xy = numel (x_xytrim) - 1 + opt;
        x_xytrim = center (x_xytrim(:), 1);
        y_xytrim = center (y_xytrim(:), 1);
        c23 = (x_xytrim' * y_xytrim) ./ denom_xy;
        
        denom_x = numel(x) - 1 + opt;
        x = center (x(:), 1);
        c1 = (x' * x) ./ denom_x;

        denom_y = numel(y) - 1 + opt;
        y = center (y(:), 1);
        c4 = (y' * y) ./ denom_y;
         
        c = [c1, c23; c23, c4];
        
        return;
  
      endif
    
    endif

  endif

endfunction


%!test
%! x = rand (10);
%! cx1 = cov (x);
%! cx2 = cov (x, x);
%! assert (size (cx1) == [10, 10] && size (cx2) == [2, 2]);
%! assert (cx2 - cx2(1), zeros (2, 2), eps);

%!test
%! x = [1:3]';
%! y = [3:-1:1]';
%! assert (cov (x, y), [1 -1; -1 1]);
%! assert (cov (x, flipud (y)), ones (2, 2));
%! assert (cov ([x, y]), [1 -1; -1 1]);

%!test
%! x = single ([1:3]');
%! y = single ([3:-1:1]');
%! assert (cov (x, y), single ([1 -1; -1 1]));
%! assert (cov (x, flipud (y)), single (ones (2, 2)));
%! assert (cov ([x, y]), single ([1 -1; -1 1]));

%!test
%! x = [0 2 4];
%! y = [3 2 1];
%! z = [4 -2; -2 1];
%! assert (cov (x, y), z);
%! assert (cov (single (x), y), single (z));
%! assert (cov (x, single (y)), single (z));
%! assert (cov (single (x), single (y)), single (z));
 
%!test
%! x = [1:5];
%! c = cov (x);
%! assert (c, 2.5);

%!test
%! x = [1:5];
%! c = cov (x, 0);
%! assert (c, 2.5);
%! c = cov (x, 1);
%! assert (c, 2);
%! c = cov (x, single (1));
%! assert (c, double (2));

%!test
%! x = [5 0 3 7; 1 -5 7 3; 4 9 8 10];
%! b = [13/3 53/6 -3 17/3; 53/6 151/3 6.5 145/6; -3 6.5 7 1; 17/3 145/6 1 37/3];
%! assert (cov (x), b, 50*eps);

%!test
%! x = [3 6 4];
%! y = [7 12 -9];
%! assert (cov (x, y), [7, 20.5; 20.5, 361]./3, 50*eps);

%!test
%! x = [2 0 -9; 3 4 1];
%! y = [5 2 6; -4 4 9];
%! assert (cov (x, y), [66.5, -20.8; -20.8, 58.4]./3, 50*eps);

%!test
%! x = [1 3 -7; 3 9 2; -5 4 6];
%! assert (cov (x, 1), [104 46 -92; 46 62 47; -92 47 266]./9, 50*eps);

%!test
%! x = [1 0; 1 0];
%! y = [1 2; 1 1];
%! assert (cov (x, y), [1/3 -1/6; -1/6 0.25], 50*eps);
%! assert (cov (x, y(:)), [1/3 -1/6; -1/6 0.25], 50*eps);
%! assert (cov (x, y(:)'), [1/3 -1/6; -1/6 0.25], 50*eps);
%! assert (cov (x', y(:)), [1/3 1/6; 1/6 0.25], 50*eps);
%! assert (cov (x(:), y), [1/3 -1/6; -1/6 0.25], 50*eps);
%! assert (cov (x(:)', y), [1/3 -1/6; -1/6 0.25], 50*eps);

%!assert (cov (5), 0)
%!assert (cov (single (5)), single (0))
%!assert (cov (1, 3), zeros (2, 2))
%!assert (cov (5, 0), 0)
%!assert (cov (5, 1), 0)
%!assert (cov (5, 2), zeros (2, 2))
%!assert (cov (5, 99), zeros (2, 2))
%!assert (cov (logical(0), logical(0)), double(0))
%!assert (cov (0, logical(0)), double(0))
%!assert (cov (logical(0), 0), double(0))
%!assert (cov (logical([0 1; 1 0]), logical([0 1; 1 0])), double ([1 1;1 1]./3))

## Test empty and NaN handling (bug #48690)
## TODO: verify compatibily for matlab > 2016b
!assert (cov ([]), NaN)
%!assert (cov (single ([])), single (NaN))
%!assert (cov ([], []), NaN (2, 2))
%!assert (cov (single ([]), single([])),  single (NaN (2, 2)))
%!assert (cov ([], single ([])), single (NaN (2, 2)))
%!assert (cov (single ([]), []), single (NaN (2, 2)))
%!assert (cov (ones(2, 0)), [])
%!assert (cov (ones(0, 2)), NaN (2, 2)) 
%!assert (cov (ones(0, 6)), NaN (6, 6)) 
%!assert (cov (ones(2, 0), []), NaN (2, 2))
%!assert (cov (NaN), NaN)
%!assert (cov (NaN, NaN), NaN (2, 2))
%!assert (cov (5, NaN), [0, NaN; NaN, NaN])
%!assert (cov (NaN, 5), [NaN, NaN; NaN, 0])
%!assert (cov (single (NaN)), single (NaN))
%!assert (cov (NaN (2, 2)), NaN (2, 2))
%!assert (cov (single (NaN (2, 2))), single (NaN (2, 2)))
%!assert (cov (NaN(2, 9)), NaN(9, 9))
%!assert (cov (NaN(9, 2)), NaN(2, 2))
%!assert (cov ([NaN, 1, 2, NaN]), NaN)
%!assert (cov ([1, NaN, 2, NaN]), NaN)

## Test nan handling parameter, 1 input
%!test
%! x = [1 3 -7; NaN 9 NaN; -5 4 6];
%! y1 = [NaN NaN NaN;NaN 31/3 NaN;NaN NaN NaN];
%! y2 = [28 NaN -15;NaN NaN NaN;-15 NaN 103/3];
%! y3 = [18 -3 -39; -3 0.5 6.5; -39 6.5 84.5];
%! assert (cov (x), y1, 50*eps);
%! assert (cov (x'), y2, 50*eps);
%! assert (cov (x, 'includenan'), y1, 50*eps);
%! assert (cov (x', 'includenan'), y2, 50*eps);
%! assert (cov (x, 'omitrows'), y3, 50*eps);
%! assert (cov (x', 'omitrows'), zeros(3, 3), 50*eps);
%! y3(2,2) = 31/3;
%! assert (cov (x, 'partialrows'), y3, 50*eps);
%! y2(isnan (y2)) = 0;
%! assert (cov (x', 'partialrows'), y2, 50*eps);

## Test nan handling parameter, 2 inputs
%!test
%! x = magic (3);
%! x(1) = NaN;
%! y = magic (3)';
%! assert (cov (x, y), [NaN, NaN; NaN, 7.5]);
%! assert (cov (x', y), [NaN, NaN; NaN, 7.5]);
%! assert (cov (x, y'), [NaN, NaN; NaN, 7.5]);
%! assert (cov (x', y'), [NaN, NaN; NaN, 7.5]);
%! assert (cov (x, y, 'omitrows'), [57/8 303/56; 303/56 57/8]);
%! assert (cov (x', y, 'omitrows'), [57/8 57/8; 57/8 57/8]);
%! assert (cov (x, y', 'omitrows'), [57/8 57/8; 57/8 57/8]);
%! assert (cov (x', y', 'omitrows'), [57/8 303/56; 303/56 57/8]);
%! assert (cov (x, y, 'partialrows'), [57/8 303/56; 303/56 7.5]);
%! assert (cov (x', y, 'partialrows'), [57/8 57/8; 57/8 7.5]);
%! assert (cov (x, y', 'partialrows'), [57/8 57/8; 57/8 7.5]);
%! assert (cov (x', y', 'partialrows'), [57/8 303/56; 303/56 7.5]);
%! assert (cov (y, x), [7.5, NaN; NaN, NaN]);
%! assert (cov (y', x), [7.5, NaN; NaN, NaN]);
%! assert (cov (y, x'), [7.5, NaN; NaN, NaN]);
%! assert (cov (y', x'), [7.5, NaN; NaN, NaN]);
%! assert (cov (y, x, 'omitrows'), [57/8 303/56; 303/56 57/8]);
%! assert (cov (y', x, 'omitrows'), [57/8 57/8; 57/8 57/8]);
%! assert (cov (y, x', 'omitrows'), [57/8 57/8; 57/8 57/8]);
%! assert (cov (y', x', 'omitrows'), [57/8 303/56; 303/56 57/8]);
%! assert (cov (y, x, 'partialrows'), [7.5 303/56; 303/56 57/8]);
%! assert (cov (y', x, 'partialrows'), [7.5 57/8; 57/8 57/8]);
%! assert (cov (y, x', 'partialrows'), [7.5 57/8; 57/8 57/8]);
%! assert (cov (y', x', 'partialrows'), [7.5 303/56; 303/56 57/8]);

## Test nan handling parameter, 2 inputs, vectors
%!test
%! x = [1:5];
%! y = [10:-2:2];
%! assert (cov (x, y), [2.5 -5; -5 10]);
%! assert (cov ([x NaN], [y 1]), [NaN NaN; NaN, 73/6], 50*eps);
%! assert (cov ([x NaN], [y 1], 'omitrows'), [2.5 -5; -5 10],50*eps);
%! assert (cov ([x NaN], [y 1], 'partialrows'), [2.5 -5; -5 73/6],50*eps);
%! assert (cov ([x 1], [y NaN]), [8/3 NaN; NaN, NaN],50*eps);
%! assert (cov ([x 1], [y NaN], 'omitrows'), [2.5 -5; -5 10],50*eps);
%! assert (cov ([x 1], [y NaN], 'partialrows'), [8/3 -5; -5 10],50*eps);
%! assert (cov ([NaN x], [y NaN], 'omitrows'), [5/3 -10/3; -10/3 20/3],50*eps);
%! assert (cov ([NaN x], [y NaN], 'partialrows'), [2.5 -10/3; -10/3 10],50*eps);

## Test nan handling parameter, one matrix trimmed to vector
%!test
%! x = magic(3);
%! y = magic (3) - 2;
%! assert (cov (x, y), [7.5 7.5; 7.5 7.5]);
%! x(3:4) = NaN;
%! assert (cov (x), [NaN, NaN, NaN; NaN, NaN, NaN; NaN, NaN, 7]);
%! assert (cov (x, 'omitrows'), zeros (3, 3));
%! assert (cov (x, 'partialrows'), [12.5 0 -2.5; 0 8 -10; -2.5 -10 7]);
%! assert (cov (x, y), [NaN, NaN; NaN, 7.5]);
%! assert (cov (x, y, 'omitrows'), [46/7, 46/7; 46/7, 46/7], 50*eps);
%! assert (cov (x, y, 'partialrows'), [46/7, 46/7; 46/7, 7.5], 50*eps);
%! assert (cov (y, x), [7.5, NaN; NaN, NaN]);
%! assert (cov (y, x, 'omitrows'), [46/7, 46/7; 46/7, 46/7], 50*eps);
%! assert (cov (y, x, 'partialrows'), [7.5, 46/7; 46/7, 46/7], 50*eps);
 
## Test input validation
%!error cov ()
%!error cov (1, 2, 3, 4)
%!error cov (5,[1 2])
%!error cov ([1; 2], ["A", "B"])
%!error cov (ones (2, 2, 2))
%!error cov (ones (2, 2), ones (2, 2, 2))
%!error cov (ones (2, 2), ones (3, 2))
%!error cov (1, [])
%!error cov (ones (1, 0, 2))
%!error cov ([1, 2],ones(1, 0, 2))
%!error cov (1, 2, [])
%!error cov (1, 2, 1, [])
