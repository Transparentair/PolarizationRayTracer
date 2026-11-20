function [Eout,kout,info] = prt_interface(E,kin,N,n1,n2,mode,layers,wavelength)
% Batched isotropic interface, arrays 3-by-N, global Cartesian fields.
% Convention: Re{E exp(i*omega*t-i*k.r)}, passive index n-i*kappa.
% p = cross(k,s); rs/rp signs must be interpreted in these moving bases.
validateattributes(n1,{'numeric'},{'scalar','real','positive','finite'});
validateattributes(n2,{'numeric'},{'scalar','finite'});
assert(real(n2)>0 && imag(n2)<=0,'prt:index','Use passive n-i*kappa.');
% MATLAB can represent a zero imaginary part with complex storage.
if imag(n2)==0, n2=real(n2); end
assert(isequal(size(E),size(kin)) && size(E,1)==3,'Expected 3-by-N arrays.');
assert(isreal(kin) && isreal(N) && all(isfinite(kin(:))) && all(isfinite(N(:))));
assert(all(vecnorm(kin)>0) && all(vecnorm(N)>0),'Zero direction/normal.');
k=kin./vecnorm(kin); N=N./vecnorm(N);
if size(N,2)==1, N=repmat(N,1,size(k,2)); end
assert(isequal(size(N),size(k)));
flip=sum(k.*N,1)>0; N(:,flip)=-N(:,flip);
ci=max(0,min(1,-sum(k.*N,1)));
s=cross(k,N,1); sn=vecnorm(s);
axial=sn<1e-12;
if any(axial)
    ka=k(:,axial); [~,idx]=min(abs(ka),[],1);
    ref=zeros(size(ka)); ref(sub2ind(size(ref),idx,1:size(ka,2)))=1;
    s(:,axial)=cross(ka,ref,1);
end
s=s./vecnorm(s); pin=cross(k,s,1);
% q = n2*cos(theta_t), choose passive/forward square-root branch.
q=sqrt(complex(n2^2-n1^2*(1-ci.^2)));
bad=imag(q)>0 | (imag(q)==0 & real(q)<0); q(bad)=-q(bad);
ds=n1*ci+q; dp=n2^2*ci+n1*q;
rs=(n1*ci-q)./ds; rp=(n2^2*ci-n1*q)./dp;
ts=2*n1*ci./ds; tp=2*n1*n2*ci./dp;
if strcmpi(mode,'ideal-mirror')
    assert(nargin<7 || isempty(layers),'请用真实基底折射率计算镀膜反射。');
    rs(:)=-1; rp(:)=1; ts(:)=0; tp(:)=0;
elseif nargin>=7 && ~isempty(layers)
    [rs,rp,ts,tp]=prt_coating(n1,n2,ci,layers,wavelength);
elseif n2==n1
    rs(:)=0; rp(:)=0; ts(:)=1; tp(:)=1;
elseif any(abs(ds)<eps | abs(dp)<eps)
    error('prt:grazing','Singular grazing interface.');
end
tir=isreal(n2) & n1^2*(1-ci.^2)>real(n2)^2;
switch lower(mode)
    case {'reflect','ideal-mirror'}
        kout=k+2*ci.*N; fs=rs; fp=rp;
    case 'transmit'
        if ~isreal(n2), error('prt:absorbing','Absorbing-medium ray propagation is unsupported.'); end
        if any(tir), error('prt:TIR','No propagating transmitted ray: use reflect mode.'); end
        kout=(n1/n2)*k+((n1/n2)*ci-real(q)/n2).*N;
        fs=ts; fp=tp;
    otherwise
        error('prt:mode','Unknown interface mode.');
end
kout=kout./vecnorm(kout); pout=cross(kout,s,1);
Eout=s.*(fs.*sum(s.*E,1))+pout.*(fp.*sum(pin.*E,1));
info=struct('rs',rs,'rp',rp,'ts',ts,'tp',tp,'cosI',ci,'q',q,'tir',tir);
end
