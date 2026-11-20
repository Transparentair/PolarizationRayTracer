function fig = prt_gui(visibility)
%PRT_GUI Interactive sequential polarization ray tracing in millimetres.
if nargin==0, visibility='on'; end
fig=uifigure('Name','偏振光追迹','Position',[100 80 1180 760], ...
    'Color',[0.94 0.96 0.98],'Visible',visibility);
root=uigridlayout(fig,[4 1]);
root.RowHeight={42,248,'1x',26};
titleLabel=uilabel(root,'Text','偏振光追迹','FontSize',23,'FontWeight','bold');
titleLabel.FontColor=[0.09 0.20 0.32];
setup=uigridlayout(root,[2 2]);
setup.ColumnWidth={305,'1x'};
setup.RowHeight={'1x',32};
setup.Padding=[0 0 0 0];
inputPanel=uipanel(setup,'Title','入射光束');
inputPanel.Layout.Row=[1 2];
input=uigridlayout(inputPanel,[7 2]);
input.ColumnWidth={114,'1x'};
input.RowHeight={24,24,24,24,24,24,24};
input.RowSpacing=4;
input.Padding=[8 8 8 8];
[exampleRows,exampleDiameter,exampleWave]=prt_example_triplet();
diameter=numberField('光瞳直径 (mm)',exampleDiameter,[eps Inf],false);
samples=numberField('每轴采样点数',81,[3 401],true);
wave=numberField('波长 (nm)',exampleWave,[eps Inf],false);
uilabel(input,'Text','入瞳偏振');
pol=uidropdown(input,'Items',{'X 线偏振','Y 线偏振','45° 线偏振','135° 线偏振','右旋圆偏振','左旋圆偏振','自定义'}, ...
    'Tooltip','面向光源观察：右旋为顺时针，左旋为逆时针。', ...
    'ItemsData',[1 2 3 7 4 5 6],'Value',3,'ValueChangedFcn',@markDirty);
uilabel(input,'Text','偏振参数');
customButton=uibutton(input,'Text','编辑自定义偏振…','ButtonPushedFcn',@editPolarization);
uilabel(input,'Text','逐点偏振');
mapActions=uigridlayout(input,[1 2]);
mapActions.Padding=[0 0 0 0];
mapActions.ColumnWidth={'1x',56};
importTip='支持二维 Jones、三维 Jones（含全局 Kx/Ky/Kz）及 Stokes CSV。';
importButton=uibutton(mapActions,'Text','导入 CSV…','ButtonPushedFcn',@importPupil,'Tooltip',importTip);
clearMap=uibutton(mapActions,'Text','清除','Enable','off','ButtonPushedFcn',@clearPupil);
uilabel(input,'Text','光程相位');
phase=uicheckbox(input,'Text','计入相对 OPL 相位','ValueChangedFcn',@markDirty);
surfaceTable=uitable(setup,'Data',exampleRows, ...
    'ColumnName',{'Type','Comment','Radius','Thickness','Material','Coating','Clear Semi-Dia.','Conic','A1','A2','A3','A4','STOP'}, ...
    'ColumnFormat',{{'OBJ','Standard','Even Asphere','IMA'},'char','numeric','numeric','char', ...
        'numeric','numeric','numeric','numeric','numeric','numeric','numeric','logical'}, ...
    'ColumnWidth',{95,80,69,75,109,59,95,51,46,46,46,46,50}, ...
    'ColumnEditable',[true(1,5),false,true(1,7)],'CellEditCallback',@editLdeCell, ...
    'CellSelectionCallback',@rememberSelection);
updateRowNames();
surfaceTable.Layout.Row=1;
surfaceTable.Layout.Column=2;
surfaceActions=uigridlayout(setup,[1 6]);
surfaceActions.ColumnWidth={76,78,76,90,78,'1x'};
surfaceActions.Padding=[0 0 0 0];
uibutton(surfaceActions,'Text','插入光面','ButtonPushedFcn',@addSurface);
uibutton(surfaceActions,'Text','删除选中','ButtonPushedFcn',@removeSurface);
uibutton(surfaceActions,'Text','载入示例','ButtonPushedFcn',@resetModel);
uibutton(surfaceActions,'Text','编辑膜系…','ButtonPushedFcn',@editCoating);
runButton=uibutton(surfaceActions,'Text','运行追迹','ButtonPushedFcn',@runTrace,'Tag','runTrace');
runButton.BackgroundColor=[0.10 0.32 0.52];
runButton.FontColor=[1 1 1];
uilabel(surfaceActions,'Text','厚度为到下一面 · mm · 反射回程填负值');
tabs=uitabgroup(root);
resultTab=uitab(tabs,'Title','偏振结果');
sceneTab=uitab(tabs,'Title','三维布局');
results=uigridlayout(resultTab,[3 1]);
results.RowHeight={34,'1x',27};
results.Padding=[0 0 0 0];
toolbar=uigridlayout(results,[1 10]);
toolbar.ColumnWidth={55,135,65,92,86,77,78,98,77,'1x'};
toolbar.Padding=[0 0 0 0];
uilabel(toolbar,'Text','结果视图');
viewType=uidropdown(toolbar,'Tag','resultView','Items',{'XY 方位角 (°)','XY 椭圆率角 (°)','场模平方 |E|²','相对光程 (mm)','有效光线','YZ 光线轨迹','XY 偏振度','XY S3/S0','偏振椭圆图','全局 |Ex|²','全局 |Ey|²','全局 |Ez|²'}, ...
    'ItemsData',1:12,'ValueChangedFcn',@(~,~) redraw());
resultPlane=uidropdown(toolbar,'Items',{'出瞳','像面'},'ItemsData',1:2, ...
    'ValueChangedFcn',@(~,~) redraw());
uibutton(toolbar,'Text','设置出入瞳','ButtonPushedFcn',@editPupils);
compare=uicheckbox(toolbar,'Text','输入/输出','ValueChangedFcn',@(~,~) redraw());
exportMat=uibutton(toolbar,'Text','导出 MAT','Enable','off','ButtonPushedFcn',@saveMat,'Tag','exportMat');
exportCsv=uibutton(toolbar,'Text','诊断 CSV','Enable','off','ButtonPushedFcn',@saveCsv);
exportMap=uibutton(toolbar,'Text','偏振态 CSV','Enable','off','ButtonPushedFcn',@savePolarizationCsv, ...
    'Tooltip','导出二维 Jones、全局三维 Jones 或局部 Stokes，表头与导入相同。');
exportPng=uibutton(toolbar,'Text','导出 PNG','Enable','off','ButtonPushedFcn',@saveEllipse);
summary=uilabel(toolbar,'Text','','HorizontalAlignment','right');
charts=uigridlayout(results,[1 2]);
charts.Padding=[0 0 0 0];
charts.ColumnWidth={0,'1x'};
inputAxes=uiaxes(charts,'Tag','inputAxes','Visible','off');
ax=uiaxes(charts,'Tag','resultAxes');
readout=uilabel(results,'Text','','FontColor',[0.27 0.36 0.46]);
fig.WindowButtonDownFcn=@clearReadoutOnBackground;
scene=uigridlayout(sceneTab,[2 2]);
scene.RowHeight={27,'1x'};
scene.ColumnWidth={490,'1x'};
sceneNote=uilabel(scene,'Text','勾选「全局」后独立定位。XYZ/U/V：mm；Rx/Ry/Rz：°。');
sceneControls=uigridlayout(scene,[1 5]);
sceneControls.Padding=[0 0 0 0];
sceneControls.ColumnWidth={100,48,86,48,86};
sceneView=uidropdown(sceneControls,'Items',{'三维视图','XY 平面','XZ 平面','YZ 平面'}, ...
    'ItemsData',1:4,'ValueChangedFcn',@changeSceneView);
uilabel(sceneControls,'Text','起始面');
fromFace=uidropdown(sceneControls,'Items',surfaceTable.RowName,'ItemsData',0:size(exampleRows,1)-1, ...
    'Value',0,'ValueChangedFcn',@changeSceneRange);
uilabel(sceneControls,'Text','结束面');
toFace=uidropdown(sceneControls,'Items',surfaceTable.RowName,'ItemsData',0:size(exampleRows,1)-1, ...
    'Value',size(exampleRows,1)-1,'ValueChangedFcn',@changeSceneRange);
poseTable=uitable(scene,'ColumnName',{'全局','X','Y','Z','Rx','Ry','Rz','口径 U','口径 V'}, ...
    'ColumnFormat',{'logical','numeric','numeric','numeric','numeric','numeric','numeric','numeric','numeric'}, ...
    'ColumnEditable',true(1,9),'ColumnWidth',{49,49,49,49,48,48,48,53,53}, ...
    'CellEditCallback',@editPose,'Tag','globalPoses');
sceneAxes=uiaxes(scene,'Tag','sceneAxes');
status=uilabel(root,'Text','运行追迹以查看结果','FontColor',[0.27 0.36 0.46],'Tag','traceStatus');
selected=[];
coatings=cell(1,size(surfaceTable.Data,1));
globalPoses=zeros(size(surfaceTable.Data,1),9);
customForm='jones';
customValues=[1 0 0 0];
result=[];
model=[];
pupil=[];
initial=[];
pupilData=[];
pupilConfig=struct('mode','auto','values',[]);
changeSceneView([],[]);
redrawScene();

    function field=numberField(label,value,limits,integer)
        uilabel(input,'Text',label);
        field=uieditfield(input,'numeric','Value',value,'Limits',limits,'ValueChangedFcn',@markDirty);
        if integer, field.RoundFractionalValues='on'; end
    end

    function markDirty(~,~)
        result=[];
        exportMat.Enable='off';
        exportCsv.Enable='off';
        exportMap.Enable='off';
        exportPng.Enable='off';
        summary.Text='';
        legend(ax,'off');
        cla(ax,'reset');
        cla(inputAxes,'reset');
        readout.Text='';
        status.Text='设置已更新，请运行追迹';
        redrawScene();
    end

    function importPupil(~,~)
        [file,folder]=uigetfile('*.csv','逐点入瞳偏振 CSV');
        if isequal(file,0), return; end
        try
            pupilData=prt_load_pupil(fullfile(folder,file));
            importButton.Text=sprintf('%d 点…',size(pupilData.xy,2));
            importButton.Tooltip=pupilData.file;
            clearMap.Enable='on';
            pol.Enable='off';
            customButton.Enable='off';
            samples.Enable='off';
            markDirty([],[]);
            status.Text=sprintf('已导入 %d 点 %s。使用 CSV 坐标追迹。',size(pupilData.xy,2),pupilData.format);
        catch problem
            uialert(fig,problem.message,'CSV 导入失败');
        end
    end

    function clearPupil(~,~)
        pupilData=[];
        importButton.Text='导入 CSV…';
        importButton.Tooltip=importTip;
        clearMap.Enable='off';
        pol.Enable='on'; customButton.Enable='on'; samples.Enable='on';
        markDirty([],[]);
    end

    function clearReadoutOnBackground(~,~)
        hit=hittest(fig);
        if isempty(hit) || ~isprop(hit,'Tag') || ~strcmp(hit.Tag,'polarizationSamples')
            readout.Text='';
        end
    end

    function editPose(~,event)
        row=event.Indices(1);
        col=event.Indices(2);
        if col>=8 && row==1
            redrawScene();
            return
        end
        value=event.NewData;
        if (isnumeric(value) || islogical(value)) && isscalar(value) && isreal(value) && isfinite(value)
            % Start from the currently displayed, resolved coordinates.
            globalPoses(row,:)=cellfun(@double,poseTable.Data(row,:));
            if col>1, globalPoses(row,1)=1; end
            markDirty([],[]);
        else
            redrawScene();
        end
    end

    function changeSceneView(~,~)
        % The plot axes are ordered Z, X, Y so +Z is screen-right.
        angles=[0 18;90 0;0 90;0 0];
        view(sceneAxes,angles(sceneView.Value,1),angles(sceneView.Value,2));
    end

    function changeSceneRange(control,~)
        if fromFace.Value>toFace.Value
            if control==fromFace, toFace.Value=fromFace.Value;
            else, fromFace.Value=toFace.Value; end
        end
        redrawScene();
    end

    function syncSceneRange()
        labels=surfaceTable.RowName;
        last=numel(labels)-1;
        oldLast=toFace.ItemsData(end);
        if oldLast==last, return; end
        first=min(fromFace.Value,last);
        finish=min(toFace.Value,last);
        if toFace.Value==oldLast, finish=last; end
        fromFace.Items=labels; fromFace.ItemsData=0:last;
        toFace.Items=labels; toFace.ItemsData=0:last;
        fromFace.Value=first; toFace.Value=finish;
    end

    function editPupils(~,~)
        values=pupilConfig.values;
        if isempty(values)
            values=zeros(2,6);
            try
                [surfaces,startZ,stop]=prt_lde_model(surfaceTable.Data,coatings);
                [surfaces,source]=prt_scene_model(surfaces,startZ,globalPoses);
                [entrance,exit]=prt_pupil_frames(surfaces,source,stop,struct('mode','auto'));
                values(1,1:3)=entrance.origin.';
                if ~isempty(exit)
                    values(2,1:3)=exit.origin.';
                    if exit.rotation(3,3)<0, values(2,4)=180; end
                end
            catch
                % Retain zero placeholders while the prescription is incomplete.
            end
        end
        dialog=uifigure('Name','出入瞳参考面','WindowStyle','modal','Position',[220 160 590 285]);
        panel=uigridlayout(dialog,[4 1]);
        panel.RowHeight={30,66,'1x',32};
        mode=uidropdown(panel,'Items',{'同轴近轴自动','手动输入 XYZ 和出瞳姿态'}, ...
            'ItemsData',{'auto','manual'},'Value',pupilConfig.mode);
        table=uitable(panel,'Data',values,'RowName',{'入瞳','出瞳'}, ...
            'ColumnName',{'X (mm)','Y (mm)','Z (mm)','Rx (°)','Ry (°)','Rz (°)'}, ...
            'ColumnEditable',true(1,6),'CellEditCallback',@lockEntranceAngles);
        uilabel(panel,'Text',['自动模式按 STOP 求近轴瞳。离轴和折转系统请输入手动参考面。' ...
            '入瞳法线和偏振基跟随 OBJ；出瞳旋转定义局部坐标。'], 'WordWrap','on');
        uibutton(panel,'Text','应用设置','ButtonPushedFcn',@applyPupils);
        function lockEntranceAngles(~,event)
            if event.Indices(1)==1 && event.Indices(2)>3
                v=table.Data; v(1,4:6)=0; table.Data=v;
            end
        end
        function applyPupils(~,~)
            if strcmp(mode.Value,'manual') && ~all(isfinite(table.Data),'all')
                uialert(dialog,'手动参考面坐标和角度须为有限实数。','设置无效'); return
            end
            pupilConfig=struct('mode',mode.Value,'values',table.Data);
            delete(dialog);
            markDirty([],[]);
        end
    end

    function redrawScene()
        syncSceneRange();
        try
            [surfaces,startZ,~,atInfinity]=prt_lde_model(surfaceTable.Data,coatings);
            [surfaces,source,poses]=prt_scene_model(surfaces,startZ,globalPoses);
            source.atInfinity=atInfinity;
        catch problem
            cla(sceneAxes);
            title(sceneAxes,'请检查镜组输入');
            sceneNote.Text=problem.message;
            return
        end
        sceneNote.Text='勾选「全局」后独立定位。XYZ/U/V：mm；Rx/Ry/Rz：°。';
        poseTable.RowName=surfaceTable.RowName;
        poseTable.Data=[num2cell(logical(poses(:,1))),num2cell(poses(:,2:end))];
        displayDiameter=diameter.Value;
        if ~isempty(pupilData), displayDiameter=max(displayDiameter,2*max(vecnorm(pupilData.xy))); end
        prt_draw_scene(sceneAxes,surfaces,source,displayDiameter,result,initial, ...
            [fromFace.Value,toFace.Value]);
    end

    function rememberSelection(~,event)
        selected=[];
        if ~isempty(event.Indices), selected=unique(event.Indices(:,1)); end
    end

    function updateRowNames()
        count=size(surfaceTable.Data,1);
        surfaceTable.RowName=[{'0 OBJ'},arrayfun(@num2str,1:count-2,'UniformOutput',false), ...
            {sprintf('%d IMA',count-1)}];
    end

    function editLdeCell(~,event)
        row=event.Indices(1);
        col=event.Indices(2);
        rows=surfaceTable.Data;
        last=size(rows,1);
        endpointLocked=(row==1 && ~ismember(col,[2 4 5])) || ...
            (row==last && ~ismember(col,[2 7]));
        if endpointLocked || (col==1 && any(strcmp(event.NewData,{'OBJ','IMA'})))
            rows{row,col}=event.PreviousData;
        elseif col==13 && isequal(event.NewData,true)
            rows(:,13)=repmat({false},last,1);
            rows{row,13}=true;
        end
        surfaceTable.Data=rows;
        markDirty([],[]);
    end

    function editPolarization(~,~)
        forms={'jones','amplitude','ellipse','stokes','jones3d'};
        labels={{'Re(Ex)','Im(Ex)','Re(Ey)','Im(Ey)'}, ...
            {'|Ex|','|Ey|','φy − φx (°)'}, ...
            {'方位角 ψ (°)','椭圆率角 χ (°)','S0'}, ...
            {'S0','S1','S2','S3'}, ...
            {'全局 Re(Ex)','全局 Im(Ex)','全局 Re(Ey)','全局 Im(Ey)','全局 Re(Ez)','全局 Im(Ez)'}};
        defaults={[1 0 0 0],[1 1 0],[45 0 1],[1 1 0 0],[1 0 0 0 0 0]};
        index=find(strcmp(forms,customForm),1);
        defaults{index}=customValues;
        dialog=uifigure('Name','自定义入射偏振','WindowStyle','modal','Position',[230 160 455 365]);
        grid=uigridlayout(dialog,[9 2]);
        grid.ColumnWidth={130,'1x'};
        grid.RowHeight={28,28,28,28,28,28,28,'1x',30};
        uilabel(grid,'Text','描述形式');
        choice=uidropdown(grid,'Items',{'二维 Jones（局部 XY）','振幅与相位差','偏振椭圆','Stokes 向量','三维 Jones（全局 XYZ）'}, ...
            'ItemsData',1:5,'Value',index,'ValueChangedFcn',@changeForm);
        names=gobjects(1,6); entries=gobjects(1,6);
        for j=1:6
            names(j)=uilabel(grid);
            entries(j)=uieditfield(grid,'numeric');
        end
        note=uilabel(grid,'Text','S3 = 2 Im(conj(Ex) Ey)，S0 表示 XY 场模平方。','WordWrap','on');
        note.Layout.Column=[1 2];
        uibutton(grid,'Text','取消','ButtonPushedFcn',@(~,~) delete(dialog));
        uibutton(grid,'Text','应用偏振态','ButtonPushedFcn',@applyState);
        displayForm(index);

        function changeForm(~,event)
            previous=event.PreviousValue;
            defaults{previous}=[entries(1:numel(labels{previous})).Value];
            displayForm(choice.Value);
        end

        function displayForm(kind)
            count=numel(labels{kind});
            for a=1:6
                if a<=count
                    names(a).Text=labels{kind}{a};
                    entries(a).Value=defaults{kind}(a);
                    names(a).Visible='on';
                    entries(a).Visible='on';
                else
                    names(a).Visible='off';
                    entries(a).Visible='off';
                end
            end
            if kind==5
                note.Text='全局电场须垂直于 OBJ 的发射方向。改变 OBJ 姿态后请重新检查。';
            else
                note.Text='S3 = 2 Im(conj(Ex) Ey)，S0 表示局部 XY 场模平方。';
            end
        end

        function applyState(~,~)
            form=forms{choice.Value};
            values=[entries(1:numel(labels{choice.Value})).Value];
            try
                prt_polarization(form,values);
            catch problem
                uialert(dialog,problem.message,'偏振参数无效');
                return
            end
            customForm=form;
            customValues=values;
            pol.Value=6;
            markDirty([],[]);
            delete(dialog);
        end
    end

    function editCoating(~,~)
        if ~isscalar(selected) || selected<=1 || selected>=size(surfaceTable.Data,1)
            uialert(fig,'请选择一个实际光学面。','选择光面');
            return
        end
        surfaceIndex=selected;
        data=coatings{surfaceIndex};
        if isempty(data), data=zeros(0,3); end
        dialog=uifigure('Name',sprintf('第 %d 面 · 多层膜系',surfaceIndex-1), ...
            'WindowStyle','modal','Position',[250 160 505 355]);
        grid=uigridlayout(dialog,[3 1]);
        grid.RowHeight={44,'1x',32};
        uilabel(grid,'Text',['按入射侧 → 基底填写。折射基底由本行 Material 指定，' ...
            '镀膜反射面请填 MIRROR:n,k。'], 'WordWrap','on');
        selectedLayers=[];
        layerTable=uitable(grid,'Data',data,'ColumnName',{'Re(n)','κ','物理厚度 (nm)'}, ...
            'ColumnEditable',true(1,3),'ColumnWidth',{100,100,155}, ...
            'RowName','numbered','CellSelectionCallback',@selectLayer);
        tools=uigridlayout(grid,[1 4]);
        tools.Padding=[0 0 0 0];
        uibutton(tools,'Text','添加膜层','ButtonPushedFcn',@(~,~) addLayer());
        uibutton(tools,'Text','删除选中','ButtonPushedFcn',@(~,~) removeLayer());
        uibutton(tools,'Text','取消','ButtonPushedFcn',@(~,~) delete(dialog));
        uibutton(tools,'Text','应用膜系','ButtonPushedFcn',@(~,~) applyCoating());

        function selectLayer(~,event)
            selectedLayers=[];
            if ~isempty(event.Indices), selectedLayers=unique(event.Indices(:,1)); end
        end

        function addLayer()
            layerTable.Data=[layerTable.Data;1.38 0 100];
        end

        function removeLayer()
            layers=layerTable.Data;
            layers(selectedLayers(selectedLayers<=size(layers,1)),:)=[];
            layerTable.Data=layers;
            selectedLayers=[];
        end

        function applyCoating()
            layers=layerTable.Data;
            if any(~isfinite(layers(:))) || any(layers(:,1)<=0) || ...
                    any(layers(:,2)<0) || any(layers(:,3)<=0)
                uialert(dialog,'折射率实部和物理厚度须为正，κ 须非负。','膜层参数无效');
                return
            end
            coatings{surfaceIndex}=layers;
            rows=surfaceTable.Data;
            rows{surfaceIndex,6}=size(layers,1);
            surfaceTable.Data=rows;
            markDirty([],[]);
            delete(dialog);
        end
    end

    function addSurface(~,~)
        rows=surfaceTable.Data;
        before=size(rows,1);
        if isscalar(selected), before=min(before,selected+1); end
        material='AIR';
        for previous=before-1:-1:1
            candidate=strtrim(char(string(rows{previous,5})));
            if ~startsWith(candidate,'MIRROR','IgnoreCase',true)
                material=candidate;
                break
            end
        end
        new={'Standard','',0,0,material,0,Inf,0,0,0,0,0,false};
        surfaceTable.Data=[rows(1:before-1,:);new;rows(before:end,:)];
        coatings=[coatings(1:before-1),{[]},coatings(before:end)];
        globalPoses=[globalPoses(1:before-1,:);zeros(1,9);globalPoses(before:end,:)];
        selected=[];
        updateRowNames();
        markDirty([],[]);
    end

    function removeSurface(~,~)
        rows=surfaceTable.Data;
        if isempty(selected), return; end
        remove=selected(selected>1 & selected<size(rows,1));
        if isempty(remove), return; end
        rows(remove,:)=[];
        coatings(remove)=[];
        globalPoses(remove,:)=[];
        surfaceTable.Data=rows;
        updateRowNames();
        selected=[];
        markDirty([],[]);
    end

    function resetModel(~,~)
        surfaceTable.Data=exampleRows;
        updateRowNames();
        syncSceneRange();
        fromFace.Value=0;
        toFace.Value=size(exampleRows,1)-1;
        coatings=cell(1,size(surfaceTable.Data,1));
        globalPoses=zeros(size(surfaceTable.Data,1),9);
        pupilConfig=struct('mode','auto','values',[]);
        diameter.Value=exampleDiameter;
        samples.Value=81;
        wave.Value=exampleWave;
        pol.Value=3;
        phase.Value=false;
        selected=[];
        clearPupil([],[]);
    end

    function runTrace(~,~)
        markDirty([],[]);
        status.Text='正在追迹…';
        drawnow limitrate nocallbacks
        try
            options=struct('diameter',diameter.Value, ...
                'samples',samples.Value,'wavelength',wave.Value,'polarization',pol.Value, ...
                'polarizationName',pol.Items{pol.ItemsData==pol.Value},'phase',phase.Value, ...
                'customForm',customForm,'customValues',customValues,'coatings',{coatings}, ...
                'globalPoses',globalPoses,'pupilConfig',pupilConfig,'pupilData',pupilData);
            if ~isempty(pupilData), options.polarizationName=[pupilData.format ' 逐点数据']; end
            [result,model,pupil,initial,elapsed]=prt_gui_trace(surfaceTable.Data,options);
            count=numel(result.valid);
            redraw();
            redrawScene();
            exportMat.Enable='on';
            exportCsv.Enable='on';
            exportMap.Enable='on';
            if isempty(result.exitPupil)
                summary.Text=sprintf('像面 %d/%d · %.3f s',nnz(result.valid),count,elapsed);
            else
                summary.Text=sprintf('出瞳 %d/%d · 像面 %d/%d · %.3f s', ...
                    nnz(result.exitPupil.valid),count,nnz(result.valid),count,elapsed);
            end
            if any(result.valid)
                status.Text='追迹完成。偏振图使用所选平面的局部坐标。';
            else
                status.Text='没有光线到达 IMA。请检查全局位姿、光面顺序及半口径。';
            end
            if ~isempty(model.pupilStatus)
                status.Text=[model.pupilStatus ' 像面追迹沿用 OBJ 发射面。'];
            end
        catch problem
            markDirty([],[]);
            status.Text='追迹失败，请检查输入参数。';
            uialert(fig,problem.message,'追迹失败');
        end
    end

    function redraw()
        legend(ax,'off');
        cla(ax,'reset');
        legend(inputAxes,'off');
        cla(inputAxes,'reset');
        colorbar(inputAxes,'off');
        readout.Text='';
        inputAxes.Visible='off';
        charts.ColumnWidth={0,'1x'};
        exportPng.Enable='off';
        if isempty(result), return; end
        if viewType.Value==6
            plotPaths();
            enablePlotNavigation(ax);
            return
        end
        data=result.image;
        planeName='像面';
        if resultPlane.Value==1
            data=result.exitPupil;
            planeName='出瞳';
        end
        if isempty(data)
            colorbar(ax,'off');
            title(ax,'未定义有限出瞳');
            text(ax,0.5,0.5,model.pupilStatus,'Units','normalized','HorizontalAlignment','center');
            return
        end
        if compare.Value
            charts.ColumnWidth={'1x','1x'};
            inputAxes.Visible='on';
            drawPlane(inputAxes,result.input,'入瞳');
        end
        drawPlane(ax,data,planeName);
        if compare.Value && any(viewType.Value==[3 4 10:12])
            common=[min([inputAxes.CLim(1),ax.CLim(1)]),max([inputAxes.CLim(2),ax.CLim(2)])];
            clim(inputAxes,common); clim(ax,common);
        end
        if viewType.Value==9, exportPng.Enable='on'; end
    end

    function drawPlane(target,data,planeName)
        if viewType.Value==9
            plotEllipses(target,data,planeName);
            enablePlotNavigation(target);
            return
        end
        ok=data.valid;
        coords=data.xy;
        switch viewType.Value
            case {1,2}
                S=data.stokesXY;
                linear=hypot(S(2,:),S(3,:));
                polarized=hypot(linear,S(4,:));
                psi=0.5*atan2(S(3,:),S(2,:));
                psi(linear<=1e-12*S(1,:))=NaN;
                chi=0.5*atan2(S(4,:),linear);
                chi(polarized<=1e-12*S(1,:))=NaN;
                values=psi*180/pi;
                if viewType.Value==2, values=chi*180/pi; end
            case 3
                values=data.totalFieldSquared;
            case 4
                values=data.OPL;
                if any(ok), values=values-values(find(ok,1)); end
            case 5
                values=double(ok);
                ok=true(size(ok));
                coords=pupil;
            case 7
                S=data.stokesXY;
                values=min(1,vecnorm(S(2:4,:))./S(1,:));
                values(S(1,:)==0)=NaN;
            case 8
                values=data.stokesXY(4,:)./data.stokesXY(1,:);
            case {10,11,12}
                values=data.globalComponentsSquared(viewType.Value-9,:);
        end
        show=ok & isfinite(values) & all(isfinite(coords),1);
        plotIds=find(show);
        scatter(target,coords(1,show),coords(2,show),12,values(show),'filled', ...
            'Tag','polarizationSamples','ButtonDownFcn', ...
            @(~,event) showSample(event,coords,values,plotIds,planeName));
        axis(target,'equal');
        if viewType.Value==5
            xlabel(target,'入瞳 X (mm)'); ylabel(target,'入瞳 Y (mm)');
        else
            xlabel(target,[planeName ' X (mm)']); ylabel(target,[planeName ' Y (mm)']);
        end
        title(target,[planeName ' · ' viewType.Items{viewType.Value}]);
        colormap(target,parula(256));
        c=colorbar(target);
        c.TicksMode='auto';
        c.TickLabelsMode='auto';
        if viewType.Value==5
            clim(target,[0 1]); c.Ticks=[0 1]; c.TickLabels={'未完成','完成'};
        elseif viewType.Value==1
            clim(target,[-90 90]);
            c.Ticks=[-90 -45 0 45 90];
            colormap(target,hsv2rgb([linspace(0,1,257).',repmat([0.72 0.85],257,1)]));
        elseif viewType.Value==2
            clim(target,[-45 45]);
        elseif viewType.Value==7
            clim(target,[0 1]);
        elseif viewType.Value==8
            clim(target,[-1 1]);
        end
        if ~any(show)
            text(target,0.5,0.5,'没有可显示的数值','Units','normalized','HorizontalAlignment','center');
        end
        enablePlotNavigation(target);
    end

    function enablePlotNavigation(target)
        % Keep pan/zoom, without re-enabling sticky built-in data tips.
        target.Interactions=[panInteraction zoomInteraction];
        axtoolbar(target,{'pan','zoomin','zoomout','restoreview'});
    end

    function showSample(event,coords,values,plotIds,planeName)
        if isempty(plotIds), return; end
        pt=event.IntersectionPoint;
        [~,nearest]=min(sum((coords(:,plotIds)-reshape(pt(1:2),2,1)).^2,1));
        id=plotIds(nearest);
        readout.Text=sprintf('%s · 光线 %d · X %.6g mm · Y %.6g mm · %s %.6g', ...
            planeName,id,coords(1,id),coords(2,id),viewType.Items{viewType.Value},values(id));
        if compare.Value && viewType.Value==1
            output=result.image;
            if resultPlane.Value==1, output=result.exitPupil; end
            if strcmp(planeName,'入瞳')
                other=output;
                label='输出';
            else
                other=result.input;
                label='输入';
            end
            S=other.stokesXY(:,id);
            if other.valid(id) && S(1)>0 && hypot(S(2),S(3))>1e-12*S(1)
                readout.Text=sprintf('%s · 对应%s方位角 %.4f°',readout.Text,label, ...
                    0.5*atan2d(S(3),S(2)));
            else
                readout.Text=sprintf('%s · 对应%s方位角无定义',readout.Text,label);
            end
        end
    end

    function plotEllipses(target,data,planeName)
        colorbar(target,'off');
        [x,y,signS3,arrowX,arrowY]=prt_pupil_ellipses(data.xy,data.stokesXY,data.valid,pupil);
        hold(target,'on');
        colors=[0.85 0.32 0.12;0.24 0.30 0.37;0.05 0.43 0.72];
        labels={'S_3 < 0','线偏振','S_3 > 0'};
        for kind=-1:1
            ids=signS3==kind;
            if ~any(ids), continue; end
            xx=[x(:,ids);nan(1,nnz(ids));arrowX(:,ids)];
            yy=[y(:,ids);nan(1,nnz(ids));arrowY(:,ids)];
            plot(target,xx(:),yy(:),'Color',colors(kind+2,:),'LineWidth',1.2, ...
                'DisplayName',labels{kind+2});
        end
        hold(target,'off');
        axis(target,'equal');
        xlabel(target,[planeName ' X (mm)']); ylabel(target,[planeName ' Y (mm)']);
        title(target,{[planeName ' · 局部 XY 偏振椭圆'],'切向箭头表示旋向，双向箭头表示线偏振轴'});
        if isempty(signS3)
            text(target,0.5,0.5,'无可绘制的偏振分量，或光线汇聚于同一点', ...
                'Units','normalized','HorizontalAlignment','center');
        else
            legend(target,'show','Location','northeast');
        end
    end

    function saveEllipse(~,~)
        [file,folder]=uiputfile('*.png','导出偏振椭圆图','polarization_ellipses.png');
        if isequal(file,0), return; end
        try
            [~,~,extension]=fileparts(file);
            if isempty(extension)
                file=[file '.png'];
            else
                assert(strcmpi(extension,'.png'),'请选择 .png 文件名。');
            end
            exportgraphics(ax,fullfile(folder,file),'Resolution',300,'BackgroundColor','white');
            status.Text=['已导出 ' file];
        catch problem
            uialert(fig,problem.message,'导出失败');
        end
    end

    function plotPaths()
        colorbar(ax,'off');
        % Select rays roughly along the pupil Y meridian.
        [~,sortOrder]=sort(abs(pupil(1,:)));
        ids=sortOrder(1:min(model.samplesPerAxis,numel(sortOrder)));
        [~,order]=sort(pupil(2,ids));
        ids=ids(order(round(linspace(1,numel(ids),min(19,numel(ids))))));
        points=NaN(3,numel(ids),1+numel(result.history));
        points(:,:,1)=initial(:,ids);
        for h=1:numel(result.history)
            hit=result.history{h};
            if isempty(hit), break; end
            [~,cols,hitCols]=intersect(ids,hit.rayIndices,'stable');
            points(:,cols,h+1)=hit.position(:,hitCols);
        end
        z=reshape(points(3,:,:),numel(ids),[]).';
        y=reshape(points(2,:,:),numel(ids),[]).';
        plot(ax,z,y,'-o','MarkerSize',3);
        grid(ax,'on');
        xlabel(ax,'全局 Z (mm)'); ylabel(ax,'全局 Y (mm)');
        title(ax,'YZ 光线轨迹 · 光瞳中心附近采样');
    end

    function saveMat(~,~)
        [file,folder]=uiputfile('*.mat','导出追迹结果','polarization_trace.mat');
        if isequal(file,0), return; end
        try
            save(fullfile(folder,file),'result','model','pupil','initial');
            status.Text=['已导出 ' file];
        catch problem
            uialert(fig,problem.message,'导出失败');
        end
    end

    function savePolarizationCsv(~,~)
        data=result.image;
        plane='像面'; stem='image';
        if resultPlane.Value==1
            data=result.exitPupil;
            plane='出瞳'; stem='exit_pupil';
        end
        if isempty(data) || ~any(data.valid)
            uialert(fig,['没有可导出的' plane '偏振采样。'],'无有效数据'); return
        end
        formats={'Stokes','取消'};
        if ~isempty(data.jonesXY), formats=[{'Jones 3D','Jones 2D'},formats]; end
        choice=uiconfirm(fig,['导出' plane '偏振态。三维 Jones 的 E 和 K 使用全局坐标。'],'选择 CSV 格式', ...
            'Options',formats,'DefaultOption',1,'CancelOption',numel(formats));
        if strcmp(choice,'取消'), return; end
        format=choice;
        if strcmp(choice,'Jones 3D'), format='Jones3D'; end
        if strcmp(choice,'Jones 2D'), format='Jones'; end
        [file,folder]=uiputfile('*.csv','导出偏振态 CSV',[stem '_' lower(format) '.csv']);
        if isequal(file,0), return; end
        try
            count=prt_save_pupil(fullfile(folder,file),data,format);
            status.Text=sprintf('已导出 %d 个%s %s 采样点',count,plane,format);
        catch problem
            uialert(fig,problem.message,'导出失败');
        end
    end

    function saveCsv(~,~)
        [file,folder]=uiputfile('*.csv','导出逐光线数据','polarization_rays.csv');
        if isequal(file,0), return; end
        fields=[pupil;result.position;result.direction;real(result.E);imag(result.E); ...
            result.OPL;result.valid;result.stokesXY;result.totalFieldSquared];
        columns={'PupilX_mm','PupilY_mm','X_mm','Y_mm','Z_mm','Kx','Ky','Kz', ...
            'ExReal','EyReal','EzReal','ExImag','EyImag','EzImag','OPL_mm','Valid', ...
            'S0_XY','S1_XY','S2_XY','S3_XY','TotalFieldSquared'};
        fields=[fields;result.image.xy;result.image.stokesXY];
        columns=[columns,{'IMA_X_mm','IMA_Y_mm','IMA_S0','IMA_S1','IMA_S2','IMA_S3'}];
        fields=[fields;result.input.stokesXY];
        columns=[columns,{'EntrancePupil_S0','EntrancePupil_S1','EntrancePupil_S2','EntrancePupil_S3'}];
        if ~isempty(result.exitPupil)
            ep=result.exitPupil;
            fields=[fields;ep.xy;ep.valid;ep.OPL;ep.stokesXY];
            columns=[columns,{'ExitPupil_X_mm','ExitPupil_Y_mm','ExitPupil_Valid', ...
                'ExitPupil_OPL_mm','ExitPupil_S0','ExitPupil_S1','ExitPupil_S2','ExitPupil_S3'}];
        end
        try
            writetable(array2table(fields.','VariableNames',columns),fullfile(folder,file));
            status.Text=['已导出 ' file];
        catch problem
            uialert(fig,problem.message,'导出失败');
        end
    end
end
