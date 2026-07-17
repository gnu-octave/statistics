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


## MATLAB-verified parity (mvregresslike R2026a): 25 observations, 3 responses,
## a shared 2-column design, at four anchor points and complete/missing data.
%!shared Xc, Yc, Ym, anch
%! Xc = [1 -0.5382438937; ...
%!  1 0.8672321576; ...
%!  1 0.9759864635; ...
%!  1 0.3373902524; ...
%!  1 -0.9960940966; ...
%!  1 -0.5232140317; ...
%!  1 -1.297447477; ...
%!  1 0.9173885891; ...
%!  1 0.1766016286; ...
%!  1 0.7551799357; ...
%!  1 -0.5914999598; ...
%!  1 1.844389637; ...
%!  1 1.816922249; ...
%!  1 -0.1238333503; ...
%!  1 -1.110601355; ...
%!  1 -0.6809058803; ...
%!  1 0.0141693264; ...
%!  1 -0.05955061046; ...
%!  1 -0.6610110145; ...
%!  1 0.3059509151; ...
%!  1 -0.4090578458; ...
%!  1 -1.281854823; ...
%!  1 -0.2849028295; ...
%!  1 -0.06478685589; ...
%!  1 1.000383189];
%! Yc = [0.03719753357 -1.298719829 2.592641544; ...
%!  -0.8732979121 0.2860080857 0.07019509407; ...
%!  0.9664790872 2.93781363 2.07467313; ...
%!  1.928227972 0.000735060807 1.023278991; ...
%!  1.099525737 -1.790260694 1.621821039; ...
%!  0.8212389245 -1.973919957 3.698786362; ...
%!  0.1466274959 -1.562219847 2.319552949; ...
%!  1.316125907 0.7610511082 0.9228417259; ...
%!  3.19879516 1.119205548 3.035118373; ...
%!  3.324100363 1.40533639 2.521210797; ...
%!  -0.5492763312 -2.535610011 -0.009321693727; ...
%!  0.7506317092 1.042220218 -0.1580843947; ...
%!  2.65978063 1.200673839 0.1312316546; ...
%!  1.060827996 -1.429695245 2.854978118; ...
%!  0.3193218045 -0.6292560647 3.238543733; ...
%!  0.08284787504 -0.7632940769 3.223147972; ...
%!  1.943818586 -2.359810723 2.132345279; ...
%!  -1.007337725 0.4325819195 1.790333779; ...
%!  0.5483627737 -2.463255725 2.103406495; ...
%!  0.4463418718 0.4193706852 1.806900862; ...
%!  2.151700093 -0.6796218054 0.2488978108; ...
%!  1.601013282 -0.7661239912 3.004258316; ...
%!  1.51658466 -0.6677671869 2.37476419; ...
%!  0.8655905179 -1.911105684 1.609121594; ...
%!  -1.036725518 0.3884447934 0.2753353942];
%! Ym = [0.03719753357 -1.298719829 2.592641544; ...
%!  -0.8732979121 NaN 0.07019509407; ...
%!  0.9664790872 2.93781363 2.07467313; ...
%!  1.928227972 0.000735060807 1.023278991; ...
%!  1.099525737 -1.790260694 NaN; ...
%!  0.8212389245 -1.973919957 3.698786362; ...
%!  0.1466274959 NaN 2.319552949; ...
%!  1.316125907 0.7610511082 0.9228417259; ...
%!  3.19879516 1.119205548 3.035118373; ...
%!  3.324100363 1.40533639 2.521210797; ...
%!  -0.5492763312 -2.535610011 -0.009321693727; ...
%!  0.7506317092 1.042220218 -0.1580843947; ...
%!  2.65978063 1.200673839 0.1312316546; ...
%!  1.060827996 NaN 2.854978118; ...
%!  0.3193218045 -0.6292560647 3.238543733; ...
%!  0.08284787504 -0.7632940769 3.223147972; ...
%!  1.943818586 -2.359810723 2.132345279; ...
%!  -1.007337725 0.4325819195 1.790333779; ...
%!  0.5483627737 -2.463255725 NaN; ...
%!  0.4463418718 0.4193706852 1.806900862; ...
%!  2.151700093 -0.6796218054 0.2488978108; ...
%!  1.601013282 -0.7661239912 3.004258316; ...
%!  1.51658466 -0.6677671869 2.37476419; ...
%!  0.8655905179 -1.911105684 1.609121594; ...
%!  -1.036725518 0.3884447934 0.2753353942];
%! Bt = [1 -0.5 2; ...
%!  0.3 1.2 -0.8];
%! S0 = cov (Yc);
%! anch = {Bt(:), S0; Bt(:)*0, eye(3); ...
%!         Bt(:)+0.25, 1.5*S0; Bt(:)-0.4, S0+0.2*eye(3)};

%!test  # nll vs MATLAB, complete data (all three algorithms agree)
%! ref = [111.8681339; 179.5391646; 120.2543928; 116.9416069];
%! for k = 1:4
%!   for alg = {"ecm", "cwls", "mvn"}
%!     assert (mvregresslike (Xc, Yc, anch{k,1}, anch{k,2}, alg{1}), ref(k), 1e-6);
%!   endfor
%! endfor

%!test  # nll vs MATLAB, missing data: ecm/cwls use all observed, mvn is listwise
%! ref_ecm = [105.2323495; 169.1339808; 112.8053693; 110.1131617];
%! ref_mvn = [91.35955956; 146.830606; 97.1466466; 96.69358201];
%! for k = 1:4
%!   assert (mvregresslike (Xc, Ym, anch{k,1}, anch{k,2}, "ecm"), ref_ecm(k), 1e-6);
%!   assert (mvregresslike (Xc, Ym, anch{k,1}, anch{k,2}, "cwls"), ref_ecm(k), 1e-6);
%!   assert (mvregresslike (Xc, Ym, anch{k,1}, anch{k,2}, "mvn"), ref_mvn(k), 1e-6);
%! endfor

%!test  # COVB vs MATLAB (complete data) = kron (Sigma, inv (X'X))
%! for k = 1:4
%!   [~, COVB] = mvregresslike (Xc, Yc, anch{k,1}, anch{k,2}, "ecm");
%!   assert (COVB, kron (anch{k,2}, inv (Xc'*Xc)), 1e-8);
%! endfor

%!error <Invalid call> mvregresslike (1, 2, 3)
%!error <ALG must be> mvregresslike ([1;2], [1;2], 1, 1, "xxx")
%!error <Sigma must be a 2-by-2> mvregresslike (ones(3,2), ones(3,2), ones(2,2), 1)
