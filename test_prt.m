function test_prt
% Analytic and invariant tests; no Symbolic/Parallel toolboxes.
passed=0; k=[0;0;1]; N=[0;0;-1]; E=[1;1i;0]/sqrt(2);
[T,kt,f]=prt_interface(E,k,N,1,1.5,'transmit');
check(norm(T-0.8*E)<1e-12 && norm(kt-k)<1e-12,'normal transmission');
[R,kr]=prt_interface(E,k,N,1,1.5,'reflect');
check(norm(R+0.2*E)<1e-12 && norm(kr+k)<1e-12,'normal reflection');
check(abs(abs(f.rs)^2+1.5*abs(f.ts)^2-1)<1e-12,'normal power');
theta=atan(1.5); kb=[sin(theta);0;cos(theta)];
[~,~,f]=prt_interface([0;1;0],kb,N,1,1.5,'reflect');
check(abs(f.rp)<1e-12,'Brewster zero');
angles=linspace(0,pi/2-1e-5,200); kk=[sin(angles);zeros(size(angles));cos(angles)];
[EE,ko,f]=prt_interface(repmat([0;1;0],1,200),kk,N,1,1.5,'transmit');
flux=real(f.q)./f.cosI;
check(max(abs(abs(f.rs).^2+flux.*abs(f.ts).^2-1))<1e-12 && ...
    max(abs(abs(f.rp).^2+flux.*abs(f.tp).^2-1))<1e-12,'oblique s/p energy');
check(max(abs(sum(EE.*ko,1)))<1e-12 && max(abs(vecnorm(ko)-1))<1e-12,'transverse/unit');
kt=[sin(pi/3);0;cos(pi/3)];
[~,~,f]=prt_interface([0;1;0],kt,N,1.5,1,'reflect');
check(f.tir && abs(abs(f.rs)-1)<1e-12 && abs(abs(f.rp)-1)<1e-12 && imag(f.rs)~=0,'TIR phase');
caught=false; try, prt_interface([0;1;0],kt,N,1.5,1,'transmit'); catch ex, caught=strcmp(ex.identifier,'prt:TIR'); end
check(caught,'TIR transmission rejected');
[~,~,f]=prt_interface([0;1;0],kb,N,1,1.3474-14.133i,'reflect');
check(abs(f.rs)<=1 && abs(f.rp)<=1 && all(isfinite([f.rs,f.rp])),'passive metal');
s=struct('radius',559,'conic',-1,'A',zeros(1,4),'z',0,'aperture',100,'n1',1,'n2',1.5,'mode','transmit');
x=linspace(-75,75,301); P0=[x;zeros(size(x));-10*ones(size(x))]; kk=repmat(k,1,numel(x));
[P,NN,t,v]=prt_intersect(P0,kk,s);
check(all(v) && max(abs(P(3,:)-x.^2/(2*559)))<1e-10 && max(vecnorm(P-P0-kk.*t))<1e-10,'parabola intersection');
check(max(vecnorm(NN-[-x/559;zeros(size(x));ones(size(x))]./vecnorm([-x/559;zeros(size(x));ones(size(x))])))<1e-12,'analytic normal');
s.conic=-2.199; s.A=[1e-6,2e-10,0,0];
[P,~,t,v]=prt_intersect(P0,kk,s); zz=prt_sag(P(1,:),P(2,:),s.radius,s.conic,s.A);
check(all(v) && max(abs(P(3,:)-zz))<1e-10 && max(vecnorm(P-P0-kk.*t))<1e-10,'polynomial intersection');
s.radius=100; s.conic=0; s.A=zeros(1,4);
[P,~,~,v]=prt_intersect([0;0;250],[0;0;-1],s);
check(v && norm(P)<1e-10,'reject opposite sphere sheet');
s.radius=Inf; s.z=5;
[P,~,t,v]=prt_intersect([1;2;0],k,s);
check(v && t==5 && norm(P-[1;2;5])<1e-12,'exact plane');
[~,~,~,v]=prt_intersect([1;2;10],k,s); check(~v,'backward hit rejected');
[local,globalEta]=prt_overlap([1;1i;0],[1;1i;0]*exp(0.7i));
check(abs(local-1)<1e-12 && abs(globalEta-1)<1e-12,'complex overlap');
[local,globalEta]=prt_overlap([1,1;0,0;0,0],[1,-1;0,0;0,0]);
check(all(local==1) && globalEta==0,'coherent cancellation');
[psi,chi]=prt_ellipse([1,1,0;2,1i,0;0,0,0]);
check(abs(psi(1)-atan(2))<1e-12 && abs(chi(2)-pi/4)<1e-12 && isnan(chi(3)),'Stokes ellipse');
% Independent sequential global-matrix reference, non-coplanar interfaces.
normals={[0;0;-1],[1;0;-1],[0;1;-1]}; ns=[1,1.5,1.2,1];
Es=E; ks=k; GM=eye(3);
for j=1:3
    [Es,kn]=prt_interface(Es,ks,normals{j},ns(j),ns(j+1),'transmit');
    [B,~]=prt_interface(eye(3),repmat(ks,1,3),normals{j},ns(j),ns(j+1),'transmit');
    GM=B*GM; ks=kn;
end
check(norm(Es-GM*E)<1e-12,'ordered noncoplanar composition');
out=demo_prt(false); check(all(out.valid) && max(abs(sum(out.E.*out.direction,1)))<1e-12,'full lens batch');
fprintf('PASS: %d analytic/invariant checks\n',passed);
    function check(condition,label)
        assert(condition,'FAILED: %s',label); passed=passed+1;
    end
end
