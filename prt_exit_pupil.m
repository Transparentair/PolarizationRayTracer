function out = prt_exit_pupil(ray,frame,n,wavelength,includePhase)
% Signed extension in the final output medium, including virtual pupils.
R=frame.rotation;
p=R.'*(ray.position-frame.origin);
k=R.'*ray.direction;
t=-p(3,:)./k(3,:);
valid=ray.valid & abs(k(3,:))>1e-12 & isfinite(t);
out.xy=p(1:2,:)+k(1:2,:).*t;
out.position=ray.position+ray.direction.*t;
out.directionGlobal=ray.direction;
out.E=R.'*ray.E;
out.OPL=ray.OPL+n*t;
out.valid=valid;
out.xy(:,~valid)=NaN;
out.position(:,~valid)=NaN;
out.directionGlobal(:,~valid)=NaN;
out.E(:,~valid)=NaN;
out.OPL(~valid)=NaN;
if includePhase && any(valid)
    reference=out.OPL(find(valid,1));
    out.E(:,valid)=out.E(:,valid).*exp(-2i*pi*(out.OPL(valid)-reference)/wavelength);
end
end
