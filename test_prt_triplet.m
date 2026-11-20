function test_prt_triplet
[rows,diameter,wavelength]=prt_example_triplet();
options=struct('diameter',diameter,'samples',81,'wavelength',wavelength, ...
    'polarization',3,'polarizationName','45 degrees','phase',false, ...
    'pupilConfig',struct('mode','auto'));
[result,model,~,initial]=prt_gui_trace(rows,options);
assert(model.sourceFrame.atInfinity && all(isfinite(initial),'all'));
assert(numel(model.surfaces)==7 && model.stopSurface==4);
assert(all(result.valid) && ~isempty(result.exitPupil));
[~,walls]=prt_scene_mesh(model.surfaces,diameter/2);
assert(isequal(~cellfun(@isempty,walls),[true false true false true]));
xy=result.image.xy-mean(result.image.xy,2);
rms=sqrt(mean(sum(xy.^2,1)));
assert(rms<0.006);
z=result.direction(1:2,:)./result.direction(3,:);
for defocus=[-0.5,0.5]
    xyDefocused=xy+defocus*(z-mean(z,2));
    assert(sqrt(mean(sum(xyDefocused.^2,1)))>5*rms);
end
fprintf('PASS: Cooke triplet, OBJ infinity, %d rays, image RMS radius %.3f um\n',nnz(result.valid),rms*1000);
end
