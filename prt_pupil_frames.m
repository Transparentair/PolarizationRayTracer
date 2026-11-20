function [entrance,exit,status] = prt_pupil_frames(surfaces,source,stop,cfg)
% Auto: real or virtual paraxial images of STOP in centered systems.
% Manual: entrance normal follows the source, exit has its own XYZ/RxRyRz.
entrance=source; exit=[]; status='';
if strcmp(cfg.mode,'manual')
    v=cfg.values;
    validateattributes(v,{'numeric'},{'real','finite','size',[2 6]});
    [~,exit]=prt_scene_model(struct('z',0),0,[1 v(2,:) 0 0;zeros(1,9)]);
    entrance.origin=v(1,1:3).';
    return
end
optics=surfaces(1:end-1);
if isempty(stop) || stop>numel(optics)
    status='请设置 STOP，或手动输入出入瞳。'; return
end
if (isfield(surfaces,'frame') && any(~cellfun(@isempty,{surfaces.frame}))) || ...
        norm(source.rotation-eye(3),'fro')>1e-12 || norm(source.origin(1:2))>1e-12
    status='全局位姿系统请手动输入出入瞳。'; return
end
incident=eye(2); emitted=eye(2); d=1;
for h=1:numel(optics)
    s=optics(h);
    n=s.n2; nextD=d;
    if ~strcmp(s.mode,'transmit'), n=s.n1; nextD=-d; end
    c=1/s.radius+2*s.A(1);
    interface=[1 0;(d*s.n1-nextD*n)*c 1];
    if h<stop, incident=interface*incident; else, emitted=interface*emitted; end
    if h<numel(optics)
        distance=(optics(h+1).z-s.z)*nextD;
        if distance < -1e-10
            status='光面轴向顺序与传播方向不一致。'; return
        end
        space=[1 distance/n;0 1];
        if h<stop, incident=space*incident; else, emitted=space*emitted; end
    end
    d=nextD;
end
if abs(incident(1,1))<1e-12 || abs(emitted(2,2))<1e-12
    status='近轴出瞳或入瞳位于无穷远，无法显示有限瞳面。'; return
end
entrance.origin=[0;0;optics(1).z+optics(1).n1*incident(1,2)/incident(1,1)];
exit=struct('origin',[0;0;optics(end).z-d*n*emitted(1,2)/emitted(2,2)], ...
    'rotation',diag([1 d d]));
end
