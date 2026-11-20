function out = demo_prt(doPlot)
% Reuses legacy three-lens radii, but starts BEFORE the first vertex.
% This is a software demonstration, not a reconstruction of the mirror system.
if nargin<1, doPlot=true; end
[x,y]=ndgrid(linspace(-7.594,7.594,101)); mask=x.^2+y.^2<=7.594^2;
P=[x(mask).';y(mask).';-10*ones(1,nnz(mask))];
k=repmat([0;0;1],1,size(P,2)); E=repmat([1;1;0]/sqrt(2),1,size(P,2));
r=[-91.857,-32.689,35.508]; z=[0,2.514,6.914]; n=[1,1.7174,1.5008,1];
s=repmat(struct('radius',Inf,'conic',0,'A',zeros(1,4),'z',0, ...
    'aperture',8,'n1',1,'n2',1,'mode','transmit'),1,3);
for h=1:3
    s(h).radius=r(h); s(h).z=z(h); s(h).n1=n(h); s(h).n2=n(h+1);
end
tic; out=prt_trace(s,P,k,E); out.seconds=toc;
fprintf('Valid %d / %d rays, %.4f s\n',nnz(out.valid),numel(out.valid),out.seconds);
if doPlot
    [psi,chi]=prt_ellipse(out.E);
    tiledlayout(1,2); nexttile;
    scatter(P(1,:),P(2,:),8,psi*180/pi,'filled'); axis equal; colorbar; title('XY azimuth (deg)');
    nexttile; scatter(P(1,:),P(2,:),8,chi*180/pi,'filled'); axis equal; colorbar; title('XY ellipticity (deg)');
end
end
