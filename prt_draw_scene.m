function prt_draw_scene(ax,surfaces,source,diameter,result,initial,range)
% Geometry preview in global coordinates. Dashed outlines are schematic.
if nargin<7, range=[0,numel(surfaces)]; end
[az,el]=view(ax);
cla(ax);
hold(ax,'on');
grid(ax,'on');
radius=diameter/2;
if range(1)==0
    p=source.origin; d=radius*source.rotation(:,3);
    ink=[0.15 0.57 0.26];
    if isfield(source,'atInfinity') && source.atInfinity
        q=p+source.rotation(:,2)*radius*[-0.75 0 0.75];
        arrows=repmat(d,1,size(q,2));
        quiver3(ax,q(3,:),q(1,:),q(2,:),arrows(3,:),arrows(1,:),arrows(2,:),0, ...
            'Color',ink,'LineWidth',1.5);
        label=' OBJ ∞';
    else
        theta=linspace(0,2*pi,49);
        ring=p+radius*source.rotation*[cos(theta);sin(theta);zeros(size(theta))];
        plot3(ax,ring(3,:),ring(1,:),ring(2,:),'Color',ink,'LineWidth',1.5);
        quiver3(ax,p(3),p(1),p(2),d(3),d(1),d(2),0,'Color',ink,'LineWidth',1.5);
        label=' OBJ';
    end
    text(ax,p(3),p(1),p(2),label,'Color',[0.10 0.40 0.19]);
end
[faces,walls]=prt_scene_mesh(surfaces,radius);
for h=max(1,range(1)):min(numel(walls),range(2)-1)
    wall=walls{h};
    if ~isempty(wall)
        surf(ax,wall(:,:,3),wall(:,:,1),wall(:,:,2), ...
            'FaceColor',[0.36 0.67 0.86],'EdgeColor','none','FaceAlpha',0.12);
    end
end
for h=max(1,range(1)):range(2)
    s=surfaces(h);
    face=faces{h};
    X=face(:,:,1); Y=face(:,:,2); Z=face(:,:,3);
    label=sprintf(' %d',h);
    color=[0.20 0.47 0.76];
    if startsWith(s.mode,'reflect') || strcmp(s.mode,'ideal-mirror'), color=[0.80 0.40 0.10]; end
    if h==numel(surfaces), label=' IMA'; color=[0.52 0.23 0.59]; end
    if s.isStop, label=sprintf(' %d STOP',h); end
    if isfinite(s.aperture)
        surf(ax,Z,X,Y,'FaceColor',color,'EdgeColor','none','FaceAlpha',0.29);
        style='-';
    else
        style='--';
    end
    plot3(ax,Z(end,:),X(end,:),Y(end,:),style,'Color',color,'LineWidth',1.2);
    c=reshape(face(1,1,:),3,1);
    if any(~isfinite(c))
        c=[0;0;s.z];
        if isfield(s,'frame') && ~isempty(s.frame), c=s.frame.origin; end
    end
    text(ax,c(3),c(1),c(2),label,'Color',color);
end
caption=sprintf('全局布局 · 面 %d–%d',range);
if ~isempty(result)
    [paths,count]=prt_scene_rays(source,diameter,result,initial,range);
    if count>0
        line(ax,paths(3,:),paths(1,:),paths(2,:),'Color',[0.13 0.29 0.48],'LineWidth',1.15);
    end
    caption=sprintf('%s · 显示 %d 条光线',caption,count);
end
hold(ax,'off');
axis(ax,'equal');
axis(ax,'tight');
xlabel(ax,'Z (mm)'); ylabel(ax,'X (mm)'); zlabel(ax,'Y (mm)');
title(ax,caption);
view(ax,az,el);
end
