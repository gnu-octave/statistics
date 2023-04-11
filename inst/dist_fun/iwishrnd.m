## Copyright (C) 2013 Nir Krakauer <nkrakauer@ccny.cuny.edu>
##
## This file is part of the statistics package for GNU Octave.
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
## @deftypefn  {statistics} {[@var{W}, @var{DI}] =} iwishrnd (@var{Tau}, @var{df}, @var{DI}, @var{n}=1)
##
## Return a random matrix sampled from the inverse Wishart distribution with
## given parameters.
##
## Inputs: the @math{p x p} positive definite matrix @var{Tau} and scalar
## degrees of freedom parameter @var{df} (and optionally the transposed Cholesky
## factor @var{DI} of @var{Sigma} = @code{inv(Tau)}).
##
## @var{df} can be non-integer as long as @math{@var{df} > d}
##
## Output: a random @math{p x p}  matrix @var{W} from the inverse
## Wishart(@var{Tau}, @var{df}) distribution. (@code{inv(W)} is from the
## Wishart(@code{inv(Tau)}, @var{df}) distribution.) If @var{n} > 1,
## then @var{W} is @var{p} x @var{p} x @var{n} and holds @var{n} such random
## matrices. (Optionally, the transposed Cholesky factor @var{DI} of @var{Sigma}
## is also returned.)
##
## Averaged across many samples, the mean of @var{W} should approach
## @var{Tau} / (@var{df} - @var{p} - 1).
##
## @subheading References
##
## @enumerate
## @item
## Yu-Cheng Ku and Peter Bloomfield (2010), Generating Random Wishart Matrices
## with Fractional Degrees of Freedom in OX,
## http://www.gwu.edu/~forcpgm/YuChengKu-030510final-WishartYu-ChengKu.pdf
## @end enumerate
##
## @seealso{iwishpdf, wishpdf, wishrnd}
## @end deftypefn

function [W, DI] = iwishrnd (Tau, df, DI, n = 1)

  if (nargin < 2)
    print_usage ();
  endif

  if (nargin < 3 || isempty (DI))
    try
      D = chol (inv (Tau));
    catch
      error (strcat (["iwishrnd: Cholesky decomposition failed;"], ...
                     [" TAU probably not positive definite."]));
    end_try_catch
    DI = D';
  else
    D = DI';
  endif

  w = wishrnd ([], df, D, n);

  if (n > 1)
    p = size (D, 1);
    W = nan (p, p, n);
  endif

  for i = 1:n
    W(:, :, i) = inv (w(:, :, i));
  endfor

endfunction



%!assert(size (iwishrnd (1,2,1)), [1, 1]);
%!assert(size (iwishrnd ([],2,1)), [1, 1]);
%!assert(size (iwishrnd ([3 1; 1 3], 2.00001, [], 1)), [2, 2]);
%!assert(size (iwishrnd (eye(2), 2, [], 3)), [2, 2, 3]);

%% Test input validation
%!error iwishrnd ()
%!error iwishrnd (1)
%!error iwishrnd ([-3 1; 1 3],1)
%!error iwishrnd ([1; 1],1)
