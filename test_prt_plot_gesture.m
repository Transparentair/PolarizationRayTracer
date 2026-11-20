function test_prt_plot_gesture
state=struct('limits',[-10 10;-5 5;-2 2],'view',[0 18], ...
    'cameraAngle',30,'dims',[1 2],'pixelSize',[400 200],'pivot',[2;1;0]);
p=prt_plot_gesture(state,'pan',[40 20]);
assert(isequal(p.limits,[-12 8;-6 4;-2 2]));
z=prt_plot_gesture(state,'zoom',-1);
assert(max(abs(diff(z.limits(1:2,:),1,2)-[20;10]/1.18))<1e-12);
assert(max(abs((z.limits(1:2,2)-[2;1])-([10;5]-[2;1])/1.18))<1e-12);
b=prt_plot_gesture(state,'box',[-10;-1;0]);
assert(isequal(b.limits,[-10 2;-1 1;-2 2]));
state.dims=[2 3];
p=prt_plot_gesture(state,'pan',[40 20]);
assert(isequal(p.limits,[-10 10;-6 4;-2.4 1.6]));
state.dims=[];
o=prt_plot_gesture(state,'orbit',[50 1000]);
assert(isequal(o.view,[-16 89.5]) && isequal(o.limits,state.limits));
z=prt_plot_gesture(state,'zoom',-1);
assert(abs(z.cameraAngle-30/1.18)<1e-12 && isequal(z.limits,state.limits));
fprintf('PASS: 2D pan, pointer-anchored zoom, box zoom and 3D orbit\n');
end
