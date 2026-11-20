function [faces,walls] = prt_scene_mesh(surfaces,previewRadius)
% Build display meshes in physical XYZ. The rim uses the clear aperture.
theta=linspace(0,2*pi,49);
count=numel(surfaces);
faces=cell(1,count);
walls=cell(1,max(0,count-2));
for h=1:count
    s=surfaces(h);
    p=[0;0;s.z]; R=eye(3); center=[0;0];
    if isfield(s,'frame') && ~isempty(s.frame)
        p=s.frame.origin; R=s.frame.rotation; center=s.apertureCenter;
    end
    radius=s.aperture;
    if ~isfinite(radius), radius=previewRadius; end
    [t,r]=meshgrid(theta,linspace(0,radius,12));
    x=center(1)+r.*cos(t); y=center(2)+r.*sin(t);
    z=prt_sag(x,y,s.radius,s.conic,s.A);
    q=p+R*[x(:).';y(:).';z(:).'];
    faces{h}=permute(reshape(q,3,size(x,1),size(x,2)),[2 3 1]);
end
for h=1:numel(walls)
    front=surfaces(h); back=surfaces(h+1);
    if ~strcmp(front.mode,'transmit') || imag(front.n2)~=0 || ...
            abs(real(front.n2)-1)<=1e-12 || abs(front.n2-back.n1)>1e-12 || ...
            ~isfinite(front.aperture) || ~isfinite(back.aperture)
        continue
    end
    a=faces{h}(end,1:end-1,:);
    b=faces{h+1}(end,1:end-1,:);
    if ~all(isfinite(a),'all') || ~all(isfinite(b),'all'), continue; end
    b=alignRims(a,b);
    separation=mean(b-a,2);
    if norm(separation(:))<=1e-10 || any(sum((b-a).*separation,3)<=0), continue; end
    walls{h}=cat(1,a(:,[1:end 1],:),b(:,[1:end 1],:));
end
end

function aligned=alignRims(front,back)
% Match angular samples even when a circular surface's local frame is rolled.
front=front-mean(front,2);
backCentered=back-mean(back,2);
count=size(back,2);
best=Inf;
aligned=back;
for direction=[1 -1]
    for offset=0:count-1
        order=mod(offset+direction*(0:count-1),count)+1;
        cost=sum((front-backCentered(:,order,:)).^2,'all');
        if cost<best
            best=cost;
            aligned=back(:,order,:);
        end
    end
end
end
