function test_prt_pupils
% A plane-parallel plate shifts the apparent image of its front STOP.
rows={'OBJ','Launch',0,5,'AIR',0,Inf,0,0,0,0,0,false; ...
      'Standard','Stop',0,10,'1.5',0,4,0,0,0,0,0,true; ...
      'Standard','Rear',0,5,'AIR',0,4,0,0,0,0,0,false; ...
      'IMA','Receiver',0,0,'',0,0.2,0,0,0,0,0,false};
opts=struct('diameter',4,'samples',9,'wavelength',1550,'polarization',4, ...
    'polarizationName','Circular','phase',true,'pupilConfig',struct('mode','auto'));
[r,m,p,initial]=prt_gui_trace(rows,opts);
assert(abs(m.entrancePupilFrame.origin(3))<1e-12);
assert(abs(m.exitPupilFrame.origin(3)-10/3)<1e-12);
assert(all(r.exitPupil.valid) && nnz(r.image.valid)==1);
assert(all(initial(3,:)==-5) && max(abs(r.exitPupil.xy-p),[],'all')<1e-12);
assert(max(abs(r.exitPupil.OPL-25/3))<1e-12 && abs(r.OPL(r.valid)-20)<1e-12);
assert(all(r.exitPupil.stokesXY(4,:)>0) && r.image.stokesXY(4,r.valid)>0);

% Manual tilted output reference and a differently oriented image plane.
fold=rows([1 2 4],:);
fold(2,3:5)={0,0,'MIRROR'};
fold{3,7}=5;
opts.globalPoses=[1 0 0 -5 0 0 0 0 0;1 0 0 0 0 -45 0 0 0;1 10 0 0 0 90 0 0 0];
opts.pupilConfig=struct('mode','manual','values',[0 0 -2 0 0 0;5 0 0 0 90 0]);
[r,m]=prt_gui_trace(fold,opts);
assert(all(r.valid) && all(r.exitPupil.valid));
assert(all(abs(r.exitPupil.position(1,:)-5)<1e-12));
assert(all(abs(r.position(1,:)-10)<1e-12));
assert(max(abs(r.exitPupil.xy-r.image.xy),[],'all')<1e-12);
assert(max(abs(r.exitPupil.stokesXY-r.image.stokesXY),[],'all')<1e-12);
assert(m.entrancePupilFrame.origin(3)==-2);
fprintf('PASS: entrance reference, virtual exit pupil, independent IMA clipping and tilted detector basis\n');
end
