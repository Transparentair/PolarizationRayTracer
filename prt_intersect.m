function [P,N,t,valid,iterations] = prt_intersect(P0,k,surf)
% Forward intersection with vertex-connected axial sag, no fictitious 1e15 plane.
% surf: radius, conic, A, z, aperture (all same length unit).
k=k./vecnorm(k); count=size(k,2); tol=1e-10;
if isfield(surf,'frame') && ~isempty(surf.frame)
    pose=surf.frame;
    local=rmfield(surf,'frame');
    local.z=0;
    R=pose.rotation;
    [P,N,t,valid,iterations]=prt_intersect(R.'*(P0-pose.origin),R.'*k,local);
    P=pose.origin+R*P;
    N=R*N;
    return
end
if isfield(surf,'normal') && ~isempty(surf.normal)
    N=surf.normal(:); N=N/norm(N); den=sum(N.*k,1);
    t=sum(N.*(surf.point(:)-P0),1)./den;
    valid=abs(den)>1e-14 & t>=-tol & isfinite(t);
    t=max(0,t); P=P0+k.*t; N=repmat(N,1,count); iterations=zeros(1,count);
    if isfinite(surf.aperture)
        valid=valid & vecnorm(P-surf.point(:))<=surf.aperture;
    end
    P(:,~valid)=NaN; N(:,~valid)=NaN; t(~valid)=NaN; return
end
% Exact conic roots provide seeds even when k_z=0. Stable quadratic formula.
c=1/surf.radius; v=P0; v(3,:)=v(3,:)-surf.z;
a=c*(k(1,:).^2+k(2,:).^2+(1+surf.conic)*k(3,:).^2);
b=2*c*(v(1,:).*k(1,:)+v(2,:).*k(2,:)+(1+surf.conic)*v(3,:).*k(3,:))-2*k(3,:);
d=c*(v(1,:).^2+v(2,:).^2+(1+surf.conic)*v(3,:).^2)-2*v(3,:);
disc=b.^2-4*a.*d; sg=ones(size(b)); sg(b<0)=-1;
q=-0.5*(b+sg.*sqrt(max(0,disc)));
t1=q./a; t2=d./q; linear=abs(a)<1e-15;
t1(linear)=-d(linear)./b(linear); t2(linear)=NaN;
t1(~linear & disc<0)=NaN; t2(~linear & disc<0)=NaN;
roots=[t1;t2]; roots(roots<-tol)=NaN;
% Reject the opposite sphere/conic sheet, then choose nearest forward seed.
for j=1:2
    R=P0+k.*roots(j,:);
    [z,~,~,ok]=prt_sag(R(1,:),R(2,:),surf.radius,surf.conic,zeros(size(surf.A)));
    roots(j,~ok | abs(R(3,:)-surf.z-z)>1e-7*(1+abs(z)))=NaN;
end
t=min(roots,[],1,'omitnan');
% Polynomial departures can create a hit absent from the base conic.
plane=(surf.z-P0(3,:))./k(3,:);
fallback=~isfinite(t) & isfinite(plane) & plane>=-tol;
t(fallback)=plane(fallback); t(t<0)=0;
valid=isfinite(t); iterations=zeros(1,count); done=false(1,count);
for it=1:50
    P=P0+k.*t;
    [z,zx,zy,ok]=prt_sag(P(1,:),P(2,:),surf.radius,surf.conic,surf.A);
    f=P(3,:)-surf.z-z; df=k(3,:)-zx.*k(1,:)-zy.*k(2,:);
    done=valid & ok & abs(f)<=tol;
    active=valid & ~done; if ~any(active), break; end
    stable=ok & abs(df)>1e-14; valid(active & ~stable)=false;
    active=active & stable; step=f./df;
    % Backtrack until residual decreases and the point stays on the real sag.
    accepted=false(1,count); candidate=t;
    for bt=1:25
        tt=t-step; R=P0+k.*tt;
        [zz,~,~,vv]=prt_sag(R(1,:),R(2,:),surf.radius,surf.conic,surf.A);
        good=active & ~accepted & vv & tt>=-tol & abs(R(3,:)-surf.z-zz)<abs(f);
        candidate(good)=max(0,tt(good)); accepted=accepted | good;
        if all(accepted(active)), break; end
        step=step/2;
    end
    valid(active & ~accepted)=false; t(accepted)=candidate(accepted);
    iterations(active)=it;
end
P=P0+k.*t;
[z,zx,zy,ok]=prt_sag(P(1,:),P(2,:),surf.radius,surf.conic,surf.A);
center=[0;0];
if isfield(surf,'apertureCenter') && ~isempty(surf.apertureCenter)
    center=surf.apertureCenter;
end
valid=valid & ok & abs(P(3,:)-surf.z-z)<=tol & vecnorm(P(1:2,:)-center)<=surf.aperture;
N=[-zx;-zy;ones(1,count)]; N=N./vecnorm(N);
P(:,~valid)=NaN; N(:,~valid)=NaN; t(~valid)=NaN;
end
