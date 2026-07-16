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
## FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
## details.
##
## You should have received a copy of the GNU General Public License along with
## this program; if not, see <http://www.gnu.org/licenses/>.

## -*- texinfo -*-
## @deftypefn  {statistics} {@var{B} =} rotatefactors (@var{A})
## @deftypefnx {statistics} {@var{B} =} rotatefactors (@var{A}, @var{Name}, @var{Value}, @dots{})
## @deftypefnx {statistics} {[@var{B}, @var{T}] =} rotatefactors (@dots{})
##
## Rotate a factor-loading matrix.
##
## @code{@var{B} = rotatefactors (@var{A})} rotates the @math{D * M} factor
## loadings matrix @var{A} (@math{D} observed variables, @math{M} factors) to
## the @qcode{'varimax'} criterion and returns the rotated loadings @var{B}, the
## same size as @var{A}.
##
## @code{[@var{B}, @var{T}] = rotatefactors (@dots{})} also returns the
## @math{M * M} rotation matrix @var{T}, so that @code{@var{B} = @var{A} *
## @var{T}}.  For the orthogonal methods @var{T} is orthonormal
## (@code{@var{T}' * @var{T}} is the identity); for the oblique methods
## (@qcode{'promax'} and oblique @qcode{'procrustes'}) it is a general invertible
## matrix.
##
## The rotation is controlled by @var{Name}/@var{Value} pairs:
##
## @table @asis
## @item @qcode{'Method'}
## The rotation criterion, one of:
##
## @table @asis
## @item @qcode{'varimax'} (default)
## Orthomax with a criterion coefficient of 1; maximizes the variance of the
## squared loadings within each factor.
##
## @item @qcode{'quartimax'}
## Orthomax with a coefficient of 0; simplifies the description of each variable.
##
## @item @qcode{'equamax'}
## Orthomax with a coefficient of @math{M / 2}.
##
## @item @qcode{'parsimax'}
## Orthomax with a coefficient of @math{D (M - 1) / (D + M - 2)}.
##
## @item @qcode{'orthomax'}
## General orthomax with the coefficient given by @qcode{'Coeff'}.
##
## @item @qcode{'promax'}
## Oblique rotation obtained by fitting an oblique transformation to a target
## built from a @qcode{'varimax'} solution raised to the power @qcode{'Power'}.
##
## @item @qcode{'procrustes'}
## Rotation towards the @qcode{'Target'} matrix, either orthogonal or oblique
## according to @qcode{'Type'}.
## @end table
##
## @item @qcode{'Normalize'}
## @qcode{'on'} (default) applies Kaiser normalization (each row of @var{A} is
## scaled to unit length before the orthomax rotation and unscaled afterwards);
## @qcode{'off'} disables it.  Ignored by @qcode{'procrustes'}.
##
## @item @qcode{'Reltol'}
## Relative convergence tolerance for the iterative orthomax rotation.  The
## default is @code{sqrt (eps)}.
##
## @item @qcode{'Maxit'}
## Maximum number of iterations for the iterative orthomax rotation.  The default
## is 250.
##
## @item @qcode{'Coeff'}
## The orthomax coefficient used when @qcode{'Method'} is @qcode{'orthomax'}.
## The default is 1 (equivalent to @qcode{'varimax'}).
##
## @item @qcode{'Power'}
## The power used to build the @qcode{'promax'} target, a scalar greater than or
## equal to 1.  The default is 4.
##
## @item @qcode{'Target'}
## The target loadings matrix for @qcode{'procrustes'}, the same size as
## @var{A}.  Required for that method.
##
## @item @qcode{'Type'}
## @qcode{'orthogonal'} (default) or @qcode{'oblique'}, selecting the kind of
## @qcode{'procrustes'} rotation.
## @end table
##
## @strong{Note on the orthomax family:} for coefficients up to 1
## (@qcode{'varimax'}, @qcode{'quartimax'}, and small @qcode{'orthomax'}) the
## rotation follows the same successive-SVD iteration as MATLAB and stops at the
## same relative tolerance.  For larger coefficients (@qcode{'equamax'},
## @qcode{'parsimax'}) that iteration does not converge, so a monotonically
## convergent pairwise algorithm is used instead; it reaches the same optimum as
## MATLAB to that solution's own convergence precision.
##
## @seealso{factoran, pca, pcacov, procrustes}
## @end deftypefn

function [B, T] = rotatefactors (A, varargin)

  ## Input validation
  if (nargin < 1)
    print_usage ();
  endif
  if (! (isnumeric (A) && ismatrix (A) && ndims (A) == 2))
    error ("rotatefactors: A must be a numeric matrix.");
  endif
  if (! isreal (A))
    error ("rotatefactors: A must be real.");
  endif

  [d, m] = size (A);

  ## Parse Name/Value options
  optNames = {'Method', 'Normalize', 'Reltol', 'Maxit', 'Coeff', ...
              'Power', 'Target', 'Type'};
  dfValues = {'varimax', 'on', sqrt(eps), 250, 1, 4, [], 'orthogonal'};
  [Method, Normalize, Reltol, Maxit, Coeff, Power, Target, Type, rem] = ...
    parsePairedArguments (optNames, dfValues, varargin(:));
  if (! isempty (rem))
    error ("rotatefactors: unknown or unpaired optional argument.");
  endif

  if (! (ischar (Method) && isrow (Method)))
    error ("rotatefactors: 'Method' must be a character vector.");
  endif
  Method = lower (Method);
  if (! any (strcmp (Method, {'varimax', 'quartimax', 'equamax', ...
                              'parsimax', 'orthomax', 'promax', 'procrustes'})))
    error ("rotatefactors: unknown 'Method' '%s'.", Method);
  endif

  if (! (ischar (Normalize) && isrow (Normalize) ...
         && any (strcmpi (Normalize, {'on', 'off'}))))
    error ("rotatefactors: 'Normalize' must be 'on' or 'off'.");
  endif
  normalize = strcmpi (Normalize, 'on');

  if (! (isnumeric (Reltol) && isscalar (Reltol) && isreal (Reltol) ...
         && Reltol > 0))
    error ("rotatefactors: 'Reltol' must be a positive scalar.");
  endif
  if (! (isnumeric (Maxit) && isscalar (Maxit) && isreal (Maxit) ...
         && Maxit >= 1 && Maxit == fix (Maxit)))
    error ("rotatefactors: 'Maxit' must be a positive integer.");
  endif

  ## Dispatch on the requested method
  switch (Method)
    case 'varimax'
      [B, T] = orthomaxrotate (A, 1, normalize, Reltol, Maxit);
    case 'quartimax'
      [B, T] = orthomaxrotate (A, 0, normalize, Reltol, Maxit);
    case 'equamax'
      [B, T] = orthomaxrotate (A, m / 2, normalize, Reltol, Maxit);
    case 'parsimax'
      [B, T] = orthomaxrotate (A, d * (m - 1) / (d + m - 2), ...
                               normalize, Reltol, Maxit);
    case 'orthomax'
      if (! (isnumeric (Coeff) && isscalar (Coeff) && isreal (Coeff)))
        error ("rotatefactors: 'Coeff' must be a real scalar.");
      endif
      [B, T] = orthomaxrotate (A, Coeff, normalize, Reltol, Maxit);
    case 'promax'
      if (! (isnumeric (Power) && isscalar (Power) && isreal (Power) ...
             && Power >= 1))
        error ("rotatefactors: 'Power' must be a scalar >= 1.");
      endif
      [B, T] = promaxrotate (A, Power, normalize, Reltol, Maxit);
    case 'procrustes'
      if (isempty (Target))
        error ("rotatefactors: 'Target' is required for 'procrustes'.");
      endif
      if (! (isnumeric (Target) && isreal (Target) ...
             && isequal (size (Target), [d, m])))
        error (strcat ("rotatefactors: 'Target' must be a real matrix", ...
                       " the same size as A."));
      endif
      if (! (ischar (Type) && isrow (Type) ...
             && any (strcmpi (Type, {'orthogonal', 'oblique'}))))
        error ("rotatefactors: 'Type' must be 'orthogonal' or 'oblique'.");
      endif
      [B, T] = procrustesrotate (A, Target, lower (Type));
  endswitch

endfunction

## Orthomax criterion (to be maximized) of a loadings matrix L
function v = orthomaxcrit (L, gamma, d)
  cs = sum (L .^ 2, 1);
  v = sum (sum (L .^ 4)) - (gamma / d) * sum (cs .^ 2);
endfunction

## Orthomax rotation: successive-SVD iteration with a pairwise fallback for the
## coefficients where the SVD iteration fails to ascend the criterion.
function [B, T] = orthomaxrotate (A, gamma, normalize, reltol, maxit)

  [d, m] = size (A);
  if (m < 2)
    B = A;
    T = eye (m);
    return;
  endif

  if (normalize)
    h = sqrt (sum (A .^ 2, 2));
    An = A ./ h;
  else
    An = A;
  endif

  ## Successive-SVD (gradient-projection) iteration.  Each step replaces T by
  ## the orthonormal polar factor of the criterion gradient.  Whenever a step
  ## would fail to increase the criterion, the iteration is abandoned in favour
  ## of the always-monotone pairwise algorithm.
  T = eye (m);
  B = An;
  sPrev = 0;
  oscillated = false;
  for iter = 1:maxit
    D = diag (sum (B .^ 2, 1));
    G = An' * (B .^ 3 - (gamma / d) * B * D);
    fCur = orthomaxcrit (B, gamma, d);
    [U, S, V] = svd (G);
    Tnew = U * V';
    if (orthomaxcrit (An * Tnew, gamma, d) < fCur - 1e-12 * max (1, abs (fCur)))
      oscillated = true;
      break;
    endif
    T = Tnew;
    B = An * T;
    sCur = sum (diag (S));
    if (iter > 1 && abs (sCur - sPrev) < reltol * sCur)
      break;
    endif
    sPrev = sCur;
  endfor

  if (oscillated)
    [B, T] = pairwiserotate (An, gamma, reltol, maxit);
  endif

  B = A * T;

endfunction

## Pairwise (Jacobi) orthomax rotation: sweep over factor pairs, applying the
## planar rotation that maximizes the criterion, until every angle is negligible.
function [B, T] = pairwiserotate (An, gamma, reltol, maxit)

  [n, m] = size (An);
  T = eye (m);
  B = An;
  for iter = 1:maxit
    converged = true;
    for p = 1:m-1
      for q = p+1:m
        x = B(:,p);
        y = B(:,q);
        u = x .^ 2 - y .^ 2;
        v = 2 * x .* y;
        As = sum (u);
        Bs = sum (v);
        Cs = sum (u .^ 2 - v .^ 2);
        Ds = sum (2 * u .* v);
        num = Ds - (2 * gamma / n) * As * Bs;
        den = Cs - (gamma / n) * (As ^ 2 - Bs ^ 2);
        theta = atan2 (num, den) / 4;
        if (abs (theta) > reltol)
          converged = false;
        endif
        c = cos (theta);
        s = sin (theta);
        B(:,p) =  c * x + s * y;
        B(:,q) = -s * x + c * y;
        R = eye (m);
        R(p,p) = c;
        R(q,q) = c;
        R(p,q) = -s;
        R(q,p) = s;
        T = T * R;
      endfor
    endfor
    if (converged)
      break;
    endif
  endfor

endfunction

## Promax rotation: oblique fit to a power target built from a varimax solution.
function [B, T] = promaxrotate (A, power, normalize, reltol, maxit)

  [d, m] = size (A);
  if (m < 2)
    B = A;
    T = eye (m);
    return;
  endif

  ## Start from a varimax (orthomax coefficient 1) rotation.
  [Bv, Tv] = orthomaxrotate (A, 1, normalize, reltol, maxit);

  ## Build the power target and fit an oblique transformation to it, normalized
  ## so the implied factor correlation matrix has a unit diagonal.
  BStar = sign (Bv) .* abs (Bv) .^ power;
  Q = Bv \ BStar;
  Q = Q .* sqrt (diag (inv (Q' * Q)))';

  T = Tv * Q;
  B = A * T;

endfunction

## Procrustes rotation towards a target, orthogonal or oblique.
function [B, T] = procrustesrotate (A, Target, type)

  if (strcmp (type, 'orthogonal'))
    [U, S, V] = svd (A' * Target);
    T = U * V';
  else
    ## Oblique: least-squares fit, columns normalized so the implied factor
    ## correlation matrix has a unit diagonal.
    T = A \ Target;
    T = T .* sqrt (diag (inv (T' * T)))';
  endif
  B = A * T;

endfunction

%!demo
%! ## Rotate a three-factor loading matrix to the varimax criterion and
%! ## recover the rotation matrix.
%! A = [ 0.8, 0.2, 0.1;  0.7, 0.3, 0.0; ...
%!       0.1, 0.9, 0.2;  0.2, 0.8, 0.1; ...
%!       0.1, 0.2, 0.9;  0.0, 0.1, 0.8];
%! [B, T] = rotatefactors (A, 'Method', 'varimax');
%! B
%! ## T is orthonormal and reconstructs B from A.
%! max (abs (vec (A * T - B)))

## Reference values below are from MATLAB R2023b for
## A = reshape (mod ((1:18)*5, 11), 6, 3) - 5.

%!test
%! A = reshape (mod ((1:18)*5, 11), 6, 3) - 5;
%! B = rotatefactors (A, 'Method', 'varimax');
%! Bref = [ -0.590343326670958, -5.399542278328798,  2.120480592034490; ...
%!           5.301191898681203,  1.349821193720584, -0.274494441363606; ...
%!          -1.789611374449118, -5.388780549360982,  0.870824505437868; ...
%!           4.101923850903043,  1.360582922688400, -1.524150527960228; ...
%!          -2.988879422227277, -5.378018820393167, -0.378831581158754; ...
%!           2.902655803124883,  1.371344651656216, -2.773806614556851];
%! assert_equal (B, Bref, 1e-10);

%!test
%! A = reshape (mod ((1:18)*5, 11), 6, 3) - 5;
%! [B, T] = rotatefactors (A, 'Method', 'varimax');
%! assert_equal (A * T, B, 1e-12);
%! assert_equal (T' * T, eye (3), 1e-12);

%!test
%! A = reshape (mod ((1:18)*5, 11), 6, 3) - 5;
%! B = rotatefactors (A, 'Method', 'quartimax');
%! Bref = [ -1.079871632250411,  1.331921409298459,  5.573137591816050; ...
%!           5.135643864717125,  1.382567861335108, -1.309071504386449; ...
%!          -1.872731015323443, -0.201056474219368,  5.427011593724494; ...
%!           4.342784481644092, -0.150410022182719, -1.455197502478006; ...
%!          -2.665590398396474, -1.734034357737195,  5.280885595632937; ...
%!           3.549925098571060, -1.683387905700547, -1.601323500569563];
%! assert_equal (B, Bref, 1e-10);

%!test
%! A = reshape (mod ((1:18)*5, 11), 6, 3) - 5;
%! B = rotatefactors (A, 'Method', 'equamax');
%! Bref = [ -0.163528146620418, -5.392883888956456,  2.211348436021980; ...
%!           4.856598298899028,  0.925912265845965, -2.357146461102300; ...
%!          -1.752970474391914, -5.249881056555078,  1.538129841051886; ...
%!           3.267155971127525,  1.068915098247334, -3.030365056072396; ...
%!          -3.342412802163413, -5.106878224153711,  0.864911246081789; ...
%!           1.677713643356028,  1.211917930648702, -3.703583651042488];
%! assert_equal (B, Bref, 1e-6);

%!test
%! A = reshape (mod ((1:18)*5, 11), 6, 3) - 5;
%! B = rotatefactors (A, 'Method', 'parsimax');
%! Bref = [ -0.209779706669808, -5.381270194314001,  2.235603625524279; ...
%!           4.824100779801576,  0.835974617312815, -2.455442547795920; ...
%!          -1.809538067356745, -5.214183987228857,  1.593065387896242; ...
%!           3.224342419114638,  1.003060824397959, -3.097980785423958; ...
%!          -3.409296428043682, -5.047097780143712,  0.950527150268203; ...
%!           1.624584058427703,  1.170147031483100, -3.740519023051994];
%! assert_equal (B, Bref, 1e-6);

%!test
%! ## orthomax with Coeff 1 equals varimax.
%! A = reshape (mod ((1:18)*5, 11), 6, 3) - 5;
%! B1 = rotatefactors (A, 'Method', 'orthomax', 'Coeff', 1);
%! B2 = rotatefactors (A, 'Method', 'varimax');
%! assert_equal (B1, B2, 1e-12);

%!test
%! ## orthomax with Coeff 0 equals quartimax.
%! A = reshape (mod ((1:18)*5, 11), 6, 3) - 5;
%! B1 = rotatefactors (A, 'Method', 'orthomax', 'Coeff', 0);
%! B2 = rotatefactors (A, 'Method', 'quartimax');
%! assert_equal (B1, B2, 1e-12);

%!test
%! ## Kaiser normalization 'off'.
%! A = reshape (mod ((1:18)*5, 11), 6, 3) - 5;
%! B = rotatefactors (A, 'Method', 'orthomax', 'Coeff', 0.5, 'Normalize', 'off');
%! Bref = [ -0.018655295494772, -5.634517554136197,  1.500621175407379; ...
%!           5.076579678301866,  2.027215443018021,  0.344581365488290; ...
%!          -1.109307774177549, -5.634304393877173,  0.155081460160809; ...
%!           3.985927199619089,  2.027428603277046, -1.000958349758279; ...
%!          -2.199960252860326, -5.634091233618149, -1.190458255085760; ...
%!           2.895274720936312,  2.027641763536070, -2.346498065004848];
%! assert_equal (B, Bref, 1e-10);

%!test
%! A = reshape (mod ((1:18)*5, 11), 6, 3) - 5;
%! [B, T] = rotatefactors (A, 'Method', 'promax', 'Power', 3);
%! Bref = [  1.224548928356429, -5.243438188632538,  2.071368824506429; ...
%!           5.402710851790721,  0.067482009256478, -0.106948752887073; ...
%!          -0.495720593244403, -5.226782837470801,  0.630342425004671; ...
%!           3.682441330189890,  0.084137360418215, -1.547975152388831; ...
%!          -2.215990114845235, -5.210127486309064, -0.810683974497087; ...
%!           1.962171808589059,  0.100792711579952, -2.989001551890589];
%! assert_equal (B, Bref, 1e-10);
%! assert_equal (A * T, B, 1e-10);

%!test
%! A = reshape (mod ((1:18)*5, 11), 6, 3) - 5;
%! Target = reshape (mod ((1:18)*3, 7), 6, 3) - 3;
%! [B, T] = rotatefactors (A, 'Method', 'procrustes', 'Target', Target);
%! Bref = [ -1.222902114166934,  5.164166701208014, -2.415759239100699; ...
%!           5.403127504320427, -0.893247624515622, -0.091224192807273; ...
%!          -2.313297220899559,  5.148157805250635, -1.070106153619982; ...
%!           4.312732397587801, -0.909256520473001,  1.254428892673443; ...
%!          -3.403692327632185,  5.132148909293257,  0.275546931860735; ...
%!           3.222337290855177, -0.925265416430380,  2.600081978154161];
%! assert_equal (B, Bref, 1e-10);
%! assert_equal (T' * T, eye (3), 1e-12);

%!test
%! A = reshape (mod ((1:18)*5, 11), 6, 3) - 5;
%! Target = reshape (mod ((1:18)*3, 7), 6, 3) - 3;
%! B = rotatefactors (A, 'Method', 'procrustes', 'Target', Target, ...
%!                    'Type', 'oblique');
%! Bref = [   0.000000000000002,  -4.619308411881804, -10.370742303151150; ...
%!           94.640566010004221, -97.005476649517561,  -1.152304700350133; ...
%!          -31.546855336668067,  36.954467295054300,   0.000000000000000; ...
%!           63.093710673336162, -55.431700942581458,   9.218437602801021; ...
%!          -63.093710673336133,  78.528243001990404,  10.370742303151154; ...
%!           31.546855336668081, -13.857925235645340,  19.589179905952172];
%! assert_equal (B, Bref, 1e-8);

%!test
%! ## Single-factor input is returned unchanged.
%! A = [1; 2; 3; 4];
%! [B, T] = rotatefactors (A);
%! assert_equal (B, A);
%! assert_equal (T, 1);

%!test
%! ## Default method is varimax.
%! A = reshape (mod ((1:18)*5, 11), 6, 3) - 5;
%! assert_equal (rotatefactors (A), rotatefactors (A, 'Method', 'varimax'));

## Test input validation
%!error<Invalid call to rotatefactors> rotatefactors ()
%!error<rotatefactors: A must be a numeric matrix.> rotatefactors ({1, 2})
%!error<rotatefactors: A must be a numeric matrix.> rotatefactors (ones (2, 2, 2))
%!error<rotatefactors: A must be real.> rotatefactors ([1+2i; 3])
%!error<rotatefactors: unknown 'Method' 'foo'.> rotatefactors (ones (3), 'Method', 'foo')
%!error<rotatefactors: 'Method' must be a character vector.> ...
%! rotatefactors (ones (3), 'Method', 5)
%!error<rotatefactors: 'Normalize' must be 'on' or 'off'.> ...
%! rotatefactors (ones (3), 'Normalize', 'yes')
%!error<rotatefactors: 'Reltol' must be a positive scalar.> ...
%! rotatefactors (ones (3), 'Reltol', -1)
%!error<rotatefactors: 'Maxit' must be a positive integer.> ...
%! rotatefactors (ones (3), 'Maxit', 2.5)
%!error<rotatefactors: 'Coeff' must be a real scalar.> ...
%! rotatefactors (ones (3), 'Method', 'orthomax', 'Coeff', [1, 2])
%!error<rotatefactors: 'Power' must be a scalar \x3e= 1.> ...
%! rotatefactors (ones (3), 'Method', 'promax', 'Power', 0.5)
%!error<rotatefactors: 'Target' is required for 'procrustes'.> ...
%! rotatefactors (ones (3), 'Method', 'procrustes')
%!error<rotatefactors: 'Target' must be a real matrix the same size as A.> ...
%! rotatefactors (ones (6, 3), 'Method', 'procrustes', 'Target', ones (4, 2))
%!error<rotatefactors: 'Type' must be 'orthogonal' or 'oblique'.> ...
%! rotatefactors (ones (6, 3), 'Method', 'procrustes', 'Target', ones (6, 3), 'Type', 'x')
%!error<rotatefactors: unknown or unpaired optional argument.> ...
%! rotatefactors (ones (3), 'Bogus', 1)
