function [z,zx,zy,valid] = prt_sag(x,y,r,K,A)
% z=c*rho^2/(1+sqrt(1-(1+K)c^2*rho^2))+sum(Aj*rho^(2j)).
% A(1:4) preserve legacy convention: rho^2, rho^4, rho^6, rho^8.
if nargin<5, A=zeros(1,4); end
assert(isscalar(r) && isreal(r) && r~=0 && ~isnan(r));
assert(isscalar(K) && isreal(K) && isfinite(K));
assert(isreal(A) && all(isfinite(A(:))));
c=1/r; u=x.^2+y.^2; rad=1-(1+K)*c*c*u;
valid=rad>0; w=sqrt(max(rad,0));
z=c*u./(1+w); du=c./(2*w);
if c==0, du=zeros(size(u)); end
for j=numel(A):-1:1
    z=z+A(j)*u.^j; du=du+j*A(j)*u.^(j-1);
end
zx=2*x.*du; zy=2*y.*du;
z(~valid)=NaN; zx(~valid)=NaN; zy(~valid)=NaN;
end
