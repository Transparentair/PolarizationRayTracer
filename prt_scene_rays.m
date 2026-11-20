function [paths,count] = prt_scene_rays(source,diameter,result,initial,range)
% Pack representative paths in one NaN-separated line, in physical XYZ.
paths=zeros(3,0);
count=0;
history=result.history;
if nargin<5, range=[0,numel(history)]; end
if isempty(initial) || range(2)<=range(1) || isempty(history), return; end
seed=history{max(1,range(1))};
if isempty(seed), return; end
candidates=reshape(seed.rayIndices(seed.valid),1,[]);
if isempty(candidates), return; end

% The closest traced rays to a central sample and two pupil rings.
uv=source.rotation(:,1:2).'*(initial(:,candidates)-source.origin);
theta=(0:7)*pi/4;
unit=[cos(theta);sin(theta)];
targets=(diameter/2)*[zeros(2,1),0.48*unit,0.93*unit];
ids=zeros(1,size(targets,2));
for k=1:numel(ids)
    [~,nearest]=min(sum((uv-targets(:,k)).^2,1));
    ids(k)=candidates(nearest);
end
ids=unique(ids,'stable');
count=numel(ids);
points=nan(3,diff(range)+2,count);
if range(1)==0, points(:,1,:)=reshape(initial(:,ids),3,1,count); end
for h=max(1,range(1)):range(2)
    hit=history{h};
    if isempty(hit), break; end
    [found,where]=ismember(ids,hit.rayIndices);
    match=find(found);
    columns=where(found);
    keep=hit.valid(columns);
    if any(keep)
        points(:,h-range(1)+1,match(keep))=reshape(hit.position(:,columns(keep)),3,1,nnz(keep));
    end
end
drawn=reshape(sum(all(isfinite(points),1),2)>=2,1,[]);
count=nnz(drawn);
points=points(:,:,drawn);
paths=reshape(points,3,[]);
end
