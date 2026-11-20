function test_prt_lde
% Check the prescription conventions through the same path as the GUI.
rows={'OBJ','Launch',0,10,'AIR',0,Inf,0,0,0,0,0,false; ...
      'Standard','Front',-91.857,2.514,'1.7174',0,8,0,0,0,0,0,true; ...
      'Standard','Cement',-32.689,4.4,'1.5008',0,8,0,0,0,0,0,false; ...
      'Standard','Rear',35.508,27,'AIR',0,8,0,0,0,0,0,false; ...
      'IMA','Receiver',0,0,'',0,Inf,0,0,0,0,0,false};
opts=struct('diameter',15.188,'samples',81,'wavelength',1550,'polarization',3, ...
    'polarizationName','45 deg','phase',false);
[result,model,pupil]=prt_gui_trace(rows,opts);
assert(all(result.valid) && nnz(result.valid)==5025);
assert(max(abs([model.surfaces.z]-[0 2.514 6.914 33.914]))<1e-12);
assert(isequal([model.surfaces.n1],[1 1.7174 1.5008 1]));
assert(isequal([model.surfaces.n2],[1.7174 1.5008 1 1]));
assert(model.stopSurface==1 && all(abs(result.position(3,:)-33.914)<1e-12));
chief=find(all(pupil==0,1),1);
assert(abs(result.OPL(chief)-(10+2.514*1.7174+4.4*1.5008+27))<1e-12);

mirror=rows([1,2,5],:);
mirror(2,3:5)={0,-5,'MIRROR'};
[result,model]=prt_gui_trace(mirror,opts);
assert(all(result.valid) && max(abs(result.OPL-15))<1e-12);
assert(all(result.direction(3,:)==-1) && all(result.position(3,:)==-5));
assert(model.surfaces(2).n1==1 && max(abs(result.E(1,:)+1/sqrt(2)))<1e-12);

mirror{2,5}='MIRROR:1.3474,14.133';
opts.coatings={[],[1.38 0 100],[]};
[result,model]=prt_gui_trace(mirror,opts);
assert(all(result.valid) && isfield(result.history{1}.coating,'rs'));
assert(model.surfaces(1).n2==1.3474-14.133i && model.surfaces(2).n1==1);

rows(2,[1 7 9])={'Even Asphere',4,1e-6};
[surfaces,~,stop]=prt_lde_model(rows);
assert(surfaces(1).A(1)==1e-6 && stop==1);
opts=rmfield(opts,'coatings');
[result,~,pupil]=prt_gui_trace(rows,opts);
assert(any(result.valid) && any(~result.valid));
assert(all(vecnorm(pupil(:,result.valid))<4.1));
fprintf('PASS: LDE thickness, media, image plane, STOP, asphere and coated reflection\n');
end
