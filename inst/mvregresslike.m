## Copyright (C) 2026 Andreas Bertsatos <abertsatos@biol.uoa.gr>
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
## @deftypefn  {statistics} {@var{nlogL} =} mvregresslike (@var{X}, @var{Y}, @var{beta}, @var{Sigma}, @var{alg})
## @deftypefnx {statistics} {[@var{nlogL}, @var{COVB}] =} mvregresslike (@dots{})
##
## Negative log-likelihood for a multivariate regression model.
##
## @code{mvregresslike (@var{X}, @var{Y}, @var{beta}, @var{Sigma}, @var{alg})}
## returns the negative log-likelihood @var{nlogL} of the multivariate normal
## regression model with responses @var{Y} (an @var{n}-by-@var{d} matrix, one
## row per observation), coefficients @var{beta}, and residual covariance
## @var{Sigma} (@var{d}-by-@var{d}).
##
## @var{X} specifies the design.  It is either a numeric @var{n}-by-@var{p}
## matrix, in which case the same @var{p} predictors apply to every response and
## @var{beta} is @var{p}-by-@var{d}; or a cell array of @var{n} design matrices,
## each @var{d}-by-@var{K}, in which case @var{beta} is @var{K}-by-1.
##
## @var{alg} selects how missing responses (@code{NaN} entries of @var{Y}) are
## handled: @qcode{"ecm"} (the default) and @qcode{"cwls"} use every observed
## response through the marginal likelihood of the observed components, while
## @qcode{"mvn"} discards any observation that has a missing response.  With no
## missing data all three agree.
##
## The optional second output @var{COVB} is the covariance matrix of the
## coefficient estimates, computed as the inverse of the observed Fisher
## information at @var{beta} and @var{Sigma}.  With missing data and the
## @qcode{"ecm"}/@qcode{"cwls"} algorithms this is the standard observed-data
## covariance and can differ from @sc{matlab}'s value (which uses a different
## information convention) at the @code{1e-3} level; @var{nlogL} agrees exactly.
##
## @seealso{mvregress}
## @end deftypefn

function [nlogL, COVB] = mvregresslike (X, Y, beta, Sigma, alg)

  if (nargin < 4)
    print_usage ();
  endif
  if (nargin < 5 || isempty (alg))
    alg = "ecm";
  endif
  alg = lower (char (alg));
  if (! any (strcmp (alg, {"ecm", "cwls", "mvn"})))
    error ("mvregresslike: ALG must be 'ecm', 'cwls', or 'mvn'.");
  endif

  [Xcell, bvec, n, d, K] = mvr_design (X, Y, beta);
  if (! isequal (size (Sigma), [d, d]))
    error ("mvregresslike: Sigma must be a %d-by-%d matrix.", d, d);
  endif

  ## 'mvn' discards observations with any missing response.
  listwise = strcmp (alg, "mvn");

  nlogL = 0;
  Info = zeros (K, K);
  for i = 1:n
    o = ! isnan (Y(i, :));
    if (! any (o) || (listwise && ! all (o)))
      continue;
    endif
    Xio = Xcell{i}(o, :);
    r = Y(i, o)' - Xio * bvec;
    So = Sigma(o, o);
    nlogL += 0.5 * (sum (o) * log (2*pi) + logdet (So) + r' * (So \ r));
    Info += Xio' * (So \ Xio);
  endfor

  if (nargout > 1)
    COVB = inv (Info);
  endif

endfunction

## Normalise the design to a cell array of per-observation d-by-K matrices and
## the coefficient vector to K-by-1.
function [Xcell, bvec, n, d, K] = mvr_design (X, Y, beta)
  if (! (isnumeric (Y) && ismatrix (Y)))
    error ("mvregresslike: Y must be a numeric matrix.");
  endif
  [n, d] = size (Y);
  if (iscell (X))
    if (numel (X) != n)
      error ("mvregresslike: X must have one design matrix per observation.");
    endif
    K = columns (X{1});
    Xcell = X(:);
    bvec = beta(:);
  else
    if (rows (X) != n)
      error ("mvregresslike: X must have as many rows as Y.");
    endif
    p = columns (X);
    K = p * d;
    Xcell = cell (n, 1);
    for i = 1:n
      Xcell{i} = kron (eye (d), X(i, :));
    endfor
    bvec = beta(:);
  endif
  if (numel (bvec) != K)
    error ("mvregresslike: beta has the wrong number of elements.");
  endif
endfunction

## log(det(A)) via the Cholesky factor (A is a positive-definite covariance).
function ld = logdet (A)
  ld = 2 * sum (log (diag (chol (A))));
endfunction

%!demo
%! ## Negative log-likelihood of a two-response regression at the true params.
%! X = [ones(20,1), (1:20)'/20];
%! B = [1 -2; 0.5 3];
%! Y = X * B + 0.3 * randn (20, 2);
%! nll = mvregresslike (X, Y, B, cov (Y - X*B))

%!test  # complete data: nll agrees across algorithms
%! X = [ones(15,1), linspace(-1,1,15)'];
%! B = [2 -1 0.5; 1 0.3 -0.4];
%! Y = X * B + 0.1 * cos ((1:15)' * [1 2 3]);
%! S = cov (Y - X*B);
%! n1 = mvregresslike (X, Y, B, S, "ecm");
%! n2 = mvregresslike (X, Y, B, S, "cwls");
%! n3 = mvregresslike (X, Y, B, S, "mvn");
%! assert (n1, n2, 1e-12);
%! assert (n1, n3, 1e-12);

%!test  # COVB of complete common design equals kron (Sigma, inv (X'X))
%! X = [ones(15,1), linspace(-1,1,15)'];
%! B = [2 -1 0.5; 1 0.3 -0.4];
%! Y = X * B + 0.1 * cos ((1:15)' * [1 2 3]);
%! S = cov (Y - X*B);
%! [~, COVB] = mvregresslike (X, Y, B, S, "ecm");
%! assert (COVB, kron (S, inv (X'*X)), 1e-10);

%!test  # cell and numeric designs give the same nll
%! X = [ones(10,1), (1:10)'];
%! B = [1 2; -1 0.5];
%! Y = X * B + 0.2 * sin ((1:10)' * [1 2]);
%! S = cov (Y - X*B);
%! Xc = cell (10, 1);
%! for i = 1:10, Xc{i} = kron (eye (2), X(i,:)); end
%! assert (mvregresslike (Xc, Y, B(:), S), mvregresslike (X, Y, B, S), 1e-12);

%!error <Invalid call> mvregresslike (1, 2, 3)
%!error <ALG must be> mvregresslike ([1;2], [1;2], 1, 1, "xxx")
%!error <Sigma must be a 2-by-2> mvregresslike (ones(3,2), ones(3,2), ones(2,2), 1)
