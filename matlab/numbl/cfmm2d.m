function [U,varargout] = cfmm2d(eps,srcinfo,pg,varargin)
%
% CFMM2D — numbl override of matlab/cfmm2d.m.
%
% Same signature, same return value as the upstream matlab/cfmm2d.m, but
% the actual FMM call goes through the cfmm2d_call JS shim (which calls
% cfmm2d_w in fmm2d.{wasm,so,dylib}) instead of the MEX file. The input
% validation logic mirrors the upstream version one-for-one so the two
% are interchangeable.
%
% cfmm2d is the complex Cauchy-kernel Laplace FMM:
%   u(x) = sum_j c_j / (x - x_j) - d_j / (x - x_j)^2
% with complex charge densities c_j and complex dipole strengths d_j.
% The return values pot/grad/hess are each one complex per source/target
% (grad = d/dz, hess = d^2/dz^2), unlike rfmm2d where grad has 2
% components and hess has 3.

  sources = srcinfo.sources;
  [m,ns] = size(sources);
  assert(m==2,'The first dimension of sources must be 2');
  if(~isfield(srcinfo,'nd'))
    nd = 1;
  end
  if(isfield(srcinfo,'nd'))
    nd = srcinfo.nd;
  end

  if( nargin < 3)
    disp('Not enough input arguments, exiting\n');
    return;
  end
  if( nargin == 3 )
    nt = 0;
    pgt = 0;
    targ = zeros(2,1);
  elseif (nargin == 4)
    nt = 0;
    pgt = 0;
    targ = zeros(2,1);
  elseif (nargin == 5)
    targ = varargin{1};
    pgt = varargin{2};
    [m,nt] = size(targ);
    assert(m==2,'First dimension of targets must be 2');
  elseif (nargin == 6)
    targ = varargin{1};
    pgt = varargin{2};
    [m,nt] = size(targ);
    assert(m==2,'First dimension of targets must be 2');
  end
  ntuse = max(nt,1);

  if((pg ==0 && pgt ==0) || (ns == 0))
    disp('Nothing to compute, set eigher pg or pgt to 1 or 2');
    return;
  end

  if(isfield(srcinfo,'charges'))
    ifcharge = 1;
    charges = srcinfo.charges;
    if(nd==1), assert(length(charges)==ns,'Charges must be same length as second dimension of sources'); end
    if(nd>1), [a,b] = size(charges); assert(a==nd && b==ns,'Charges must be of shape [nd,ns] where nd is the number of densities, and ns is the number of sources'); end
  else
    ifcharge = 0;
    charges = complex(zeros(nd,ns));
  end

  if(isfield(srcinfo,'dipstr'))
    ifdipole = 1;
    dipstr = srcinfo.dipstr;
    if(nd==1), assert(length(dipstr)==ns,'Dipole strength must be same length as second dimension of sources'); end
    if(nd>1), [a,b] = size(dipstr); assert(a==nd && b==ns,'Dipstr must be of shape [nd,ns] where nd is the number of densities, and ns is the number of sources'); end
  else
    ifdipole = 0;
    dipstr = complex(zeros(nd,ns));
  end

  % Make sure targets are full-sized even when nt==0 (Fortran callee
  % still wants a non-null pointer dimensioned (2, ntuse)).
  if nt == 0
    targ = zeros(2, 1);
  end

  % Single packed call. The shim returns one complex column-major flat
  % tensor in this order:
  %   [pot(nd,ns) | grad(nd,ns) | hess(nd,ns) |
  %    pottarg(nd,ntuse) | gradtarg(nd,ntuse) | hesstarg(nd,ntuse)]
  packed = cfmm2d_call(nd, eps, ns, sources, ifcharge, charges, ...
                       ifdipole, dipstr, ...
                       pg, nt, targ, pgt);

  n_pot = nd * ns;
  n_grad = nd * ns;
  n_hess = nd * ns;
  n_pottarg = nd * ntuse;
  n_gradtarg = nd * ntuse;
  n_hesstarg = nd * ntuse;

  off = 0;
  pot = reshape(packed(off+1:off+n_pot), [nd, ns]);                         off = off + n_pot;
  grad = reshape(packed(off+1:off+n_grad), [nd, ns]);                       off = off + n_grad;
  hess = reshape(packed(off+1:off+n_hess), [nd, ns]);                       off = off + n_hess;
  pottarg = reshape(packed(off+1:off+n_pottarg), [nd, ntuse]);              off = off + n_pottarg;
  gradtarg = reshape(packed(off+1:off+n_gradtarg), [nd, ntuse]);            off = off + n_gradtarg;
  hesstarg = reshape(packed(off+1:off+n_hesstarg), [nd, ntuse]);

  U.pot = [];
  U.grad = [];
  U.hess = [];
  U.pottarg = [];
  U.gradtarg = [];
  U.hesstarg = [];
  if(pg >= 1), U.pot = squeeze(reshape(pot,[nd,ns])); end
  if(pg >= 2), U.grad = squeeze(reshape(grad,[nd,ns])); end
  if(pg >= 3), U.hess = squeeze(reshape(hess,[nd,ns])); end
  if(pgt >= 1), U.pottarg = squeeze(reshape(pottarg,[nd,nt])); end
  if(pgt >= 2), U.gradtarg = squeeze(reshape(gradtarg,[nd,nt])); end
  if(pgt >= 3), U.hesstarg = squeeze(reshape(hesstarg,[nd,nt])); end

  varargout{1} = 0;          % ier (the shim throws on nonzero)
  varargout{2} = zeros(8,1); % timeinfo (not exposed by the shim)
end
