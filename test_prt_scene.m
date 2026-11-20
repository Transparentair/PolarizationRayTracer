function test_prt_scene
% A 45-degree plane fold turns +Z into +X, using a tilted IMA.
rows={'OBJ','Source',0,5,'AIR',0,Inf,0,0,0,0,0,false; ...
      'Standard','Fold',0,0,'MIRROR',0,12,0,0,0,0,0,false; ...
      'IMA','Receiver',0,0,'',0,12,0,0,0,0,0,false};
poses=[1 2 3 -5 0 0 90 0 0; ...
       1 0 0  0 0 -45 0 0 0; ...
       1 10 0 0 0 90 0 0 0];
opts=struct('diameter',2,'samples',9,'wavelength',1550,'polarization',1, ...
    'polarizationName','X','phase',false,'globalPoses',poses);
[result,model,pupil,initial]=prt_gui_trace(rows,opts);
assert(all(result.valid));
chief=find(all(pupil==0,1),1);
assert(norm(initial(:,chief)-[2;3;-5])<1e-12);
assert(norm(result.history{1}.position(:,chief)-[2;3;2])<1e-12);
assert(norm(result.position(:,chief)-[10;3;2])<1e-12);
assert(norm(result.direction(:,chief)-[1;0;0])<1e-12);
assert(norm(result.E(:,chief)-[0;-1;0])<1e-12);
assert(abs(result.OPL(chief)-15)<1e-12);
assert(isequal(model.globalPoses,poses));

% Display decimation must preserve the actual central path and gaps.
[paths,count]=prt_scene_rays(model.sourceFrame,opts.diameter,result,initial);
assert(count==17 && size(paths,2)==count*4);
assert(max(abs(paths(:,1:3)-[2 2 10;3 3 3;-5 2 2]),[],'all')<1e-12);
assert(all(isnan(paths(:,4:4:end)),'all'));
[cropped,count]=prt_scene_rays(model.sourceFrame,opts.diameter,result,initial,[1 2]);
assert(count==17 && size(cropped,2)==count*3);
assert(max(abs(cropped(:,1:2)-[2 10;3 3;2 2]),[],'all')<1e-12);
assert(all(isnan(cropped(:,3:3:end)),'all'));
[empty,count]=prt_scene_rays(model.sourceFrame,opts.diameter,result,initial,[1 1]);
assert(count==0 && isempty(empty));
result.history{2}.valid(:)=false;
[clipped,count]=prt_scene_rays(model.sourceFrame,opts.diameter,result,initial);
assert(count==17 && all(isnan(clipped(:,3:4:end)),'all'));
result.history{1}.valid(:)=false;
[empty,count]=prt_scene_rays(model.sourceFrame,opts.diameter,result,initial);
assert(count==0 && isempty(empty));

% Decentered aperture on a parent paraboloid, not a shifted vertex.
rows{2,3}=20; rows{2,8}=-1; rows{2,7}=1;
[surfaces,startZ]=prt_lde_model(rows);
poses(2,:)=[1 0 0 0 0 0 0 3 0];
poses(3,1)=0;
[surfaces,~,resolved]=prt_scene_model(surfaces,startZ,poses);
[hit,~,~,valid]=prt_intersect([0 3;0 0;-5 -5],[0 0;0 0;1 1],surfaces(1));
assert(isequal(valid,[false true]));
assert(abs(hit(3,2)-9/40)<1e-12);
assert(isempty(surfaces(2).frame) && resolved(3,4)==0);
fprintf('PASS: global source, folded reflection, tilted IMA, OPL and off-axis aperture\n');

% A cemented doublet shares one optical interface and leaves air unfilled.
lensRows={'OBJ','Source',0,10,'AIR',0,Inf,0,0,0,0,0,false; ...
    'Standard','Front',50,4,'1.5',0,8,0,0,0,0,0,false; ...
    'Standard','Cement',0,3,'1.6',0,7,0,0,0,0,0,false; ...
    'Standard','Rear',-70,5,'AIR',0,8,0,0,0,0,0,false; ...
    'Standard','Window',0,2,'AIR',0,8,0,0,0,0,0,false; ...
    'IMA','Receiver',0,0,'',0,Inf,0,0,0,0,0,false};
[lens,sourceZ]=prt_lde_model(lensRows);
pose=zeros(6,9);
pose(3,:)=[1 0 0 4 0 0 90 0 0];
lens=prt_scene_model(lens,sourceZ,pose);
[faces,walls]=prt_scene_mesh(lens,8);
assert(numel(faces)==5 && ~isempty(walls{1}) && ~isempty(walls{2}) && isempty(walls{3}));
assert(max(abs(walls{1}(1,:,:)-faces{1}(end,:,:)),[],'all')<1e-12);
assert(max(abs(walls{1}(2,:,1:2)-(7/8)*walls{1}(1,:,1:2)),[],'all')<1e-12);
assert(all(abs(walls{1}(2,:,3)-4)<1e-12));
assert(all(isfinite(walls{2}),'all'));
assert(max(abs(walls{1}(2,1,:)-walls{1}(2,end,:)),[],'all')<1e-12);
lens(2).aperture=Inf;
[~,walls]=prt_scene_mesh(lens,8);
assert(all(cellfun(@isempty,walls)));
fprintf('PASS: glass volumes, rolled edges, air gap and infinite-aperture fallback\n');
end
