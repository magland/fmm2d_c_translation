function [U] = stfmm2d(eps,srcinfo,ifppreg,targ,ifppregtarg)
%
% STFMM2D - numbl override of matlab/stfmm2d.m.
%
% Same signature, same return value as the upstream matlab/stfmm2d.m,
% but the actual FMM call goes through the stfmm2d_call JS shim
% (which calls stfmm2d_w in fmm2d.{wasm,so,dylib}) instead of the MEX
% file.
%
% 2D Stokes FMM with Stokeslet (sigma) and (type-I) stresslet (mu, nu)
% sources. All inputs and outputs are real-valued doubles.

  sources = srcinfo.sources;
  [m,ns] = size(sources);
  assert(m==2,'The first dimension of sources must be 2');
  if(~isfield(srcinfo,'nd'))
    nd = 1;
  end
  if(isfield(srcinfo,'nd'))
    nd = srcinfo.nd;
  end

  if( nargin <= 3 )
    nt = 0;
    ifppregtarg = 0;
    targ = zeros(2,0);
  else
    if( nargin <= 4 ), ifppregtarg = 0; end
    [m,nt] = size(targ);
    assert(m==2,'First dimension of targets must be 2');
  end
  ntuse = max(nt,1);

  if(ifppreg ==0 && ifppregtarg ==0)
    disp('Nothing to compute, set eigher ifppreg or ifppregtarg to 1 or 2 or 3');
    return;
  end

  if(isfield(srcinfo,'stoklet'))
    ifstoklet = 1;
    stoklet = srcinfo.stoklet;
    if(nd == 1), [a,b] = size(squeeze(stoklet)); assert(a==2 && b==ns,'Stoklet must be of shape[2,ns], where ns is the number of sources'); end
    if(nd>1), [a,b,c] = size(stoklet); assert(a==nd && b==2 && c==ns, 'Stoklet must be of shape[nd,2,ns], where nd is number of densities, and ns is the number of sources'); end
    stoklet = reshape(stoklet,[2*nd,ns]);
  else
    ifstoklet = 0;
    stoklet = zeros(2*nd,ns);
  end

  if(isfield(srcinfo,'strslet') && isfield(srcinfo,'strsvec'))
    ifstrslet = 1;
    strslet = srcinfo.strslet;
    strsvec = srcinfo.strsvec;
    if(nd == 1), [a,b] = size(squeeze(strslet)); assert(a==2 && b==ns,'Strslet must be of shape[2,ns], where ns is the number of sources'); end
    if(nd == 1), [a,b] = size(squeeze(strsvec)); assert(a==2 && b==ns,'Strsvec must be of shape[2,ns], where ns is the number of sources'); end
    if(nd>1), [a,b,c] = size(strslet); assert(a==nd && b==2 && c==ns, 'Strslet must be of shape[nd,2,ns], where nd is number of densities, and ns is the number of sources'); end
    if(nd>1), [a,b,c] = size(strsvec); assert(a==nd && b==2 && c==ns, 'Strsvec must be of shape[nd,2,ns], where nd is number of densities, and ns is the number of sources'); end
    strslet = reshape(strslet,[2*nd,ns]);
    strsvec = reshape(strsvec,[2*nd,ns]);
  else
    ifstrslet = 0;
    strslet = zeros(2*nd,ns);
    strsvec = zeros(2*nd,ns);
  end

  if nt == 0
    targ = zeros(2,1);
  end

  % Single packed call. The shim returns one column-major flat
  % Float64 tensor in this order:
  %   [pot(2*nd,ns) | pre(nd,ns) | grad(4*nd,ns) |
  %    pottarg(2*nd,ntuse) | pretarg(nd,ntuse) | gradtarg(4*nd,ntuse)]
  packed = stfmm2d_call(nd, eps, ns, sources, ...
                        ifstoklet, stoklet, ...
                        ifstrslet, strslet, strsvec, ...
                        ifppreg, nt, targ, ifppregtarg);

  n_pot = 2 * nd * ns;
  n_pre = nd * ns;
  n_grad = 4 * nd * ns;
  n_pottarg = 2 * nd * ntuse;
  n_pretarg = nd * ntuse;
  n_gradtarg = 4 * nd * ntuse;

  off = 0;
  pot = reshape(packed(off+1:off+n_pot), [2*nd, ns]);                   off = off + n_pot;
  pre = reshape(packed(off+1:off+n_pre), [nd, ns]);                     off = off + n_pre;
  grad = reshape(packed(off+1:off+n_grad), [4*nd, ns]);                 off = off + n_grad;
  pottarg = reshape(packed(off+1:off+n_pottarg), [2*nd, ntuse]);        off = off + n_pottarg;
  pretarg = reshape(packed(off+1:off+n_pretarg), [nd, ntuse]);          off = off + n_pretarg;
  gradtarg = reshape(packed(off+1:off+n_gradtarg), [4*nd, ntuse]);

  U.pot = [];
  U.pre = [];
  U.grad = [];
  U.pottarg = [];
  U.pretarg = [];
  U.gradtarg = [];
  if(ifppreg >= 1), U.pot = squeeze(reshape(pot,[nd,2,ns])); end
  if(ifppreg >= 2), U.pre = pre; end
  if(ifppreg >= 3), U.grad = squeeze(reshape(grad,[nd,2,2,ns])); end
  if(ifppregtarg >= 1), U.pottarg = squeeze(reshape(pottarg,[nd,2,nt])); end
  if(ifppregtarg >= 2), U.pretarg = pretarg; end
  if(ifppregtarg >= 3), U.gradtarg = squeeze(reshape(gradtarg,[nd,2,2,nt])); end
end
