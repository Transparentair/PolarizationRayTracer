function state = prt_plot_gesture(state,action,delta)
% Update an axes snapshot without replotting. Limits are 3-by-2, view in deg.
switch action
    case 'pan'
        dims=state.dims;
        span=diff(state.limits(dims,:),1,2);
        shift=-reshape(delta,2,1).*span./max(reshape(state.pixelSize,2,1),1);
        state.limits(dims,:)=state.limits(dims,:)+shift;
    case 'orbit'
        state.view=state.view+[-0.32*delta(1),0.32*delta(2)];
        state.view(2)=max(-89.5,min(89.5,state.view(2)));
    case 'zoom'
        factor=1.18^max(-10,min(10,delta));
        dims=state.dims;
        if isempty(dims)
            state.cameraAngle=max(0.05,min(120,state.cameraAngle*factor));
        else
            current=state.limits(dims,:);
            pivot=min(current(:,2),max(current(:,1),state.pivot(dims)));
            state.limits(dims,:)=pivot+(current-pivot)*factor;
        end
    case 'box'
        dims=state.dims;
        current=state.limits(dims,:);
        corner=min(current(:,2),max(current(:,1),delta(dims)));
        pivot=min(current(:,2),max(current(:,1),state.pivot(dims)));
        state.limits(dims,:)=sort([pivot corner],2);
    otherwise
        error('prt:plotGesture','Unknown plot gesture.');
end
end
