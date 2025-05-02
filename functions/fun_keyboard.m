function fun_keyboard(src,evt)

%KEY FUNCTIONS
switch lower(evt.Key)
    case 'rightarrow'
        fac =  1; %direction
    case 'leftarrow'
        fac = -1; %direction
    case 'uparrow'
        ylim(ylim-0.01*diff(ylim));
        return
    case 'downarrow'
        ylim(ylim+0.01*diff(ylim));
        return
    otherwise
        return
end
%steps
steps = 1*fac; %default step
if ~isempty(evt.Modifier)
    switch lower(evt.Modifier{1}) %1st only
        case 'alt'
            steps = 20*fac;
        case 'control'
            steps = 100*fac;
        case 'shift'
            steps = 1000*fac;
    end
end

%GET DATA
g   = guidata(src);
hp  = findobj(gca,'type','line');
ind = find(g.hp==double(hp)); %handle index
if numel(ind)~=1 %else for sure not made for its purpose
    warning('Wrong axis')
    return
end

%SET DATA
g.shift(ind) = g.shift(ind)+ g.dt*steps;
set(hp,'xdata',g.t{ind}+g.shift(ind))
xlabel({'Time [s]',sprintf('Shift Vm: %.3f',g.shift(ind))})
drawnow

%APPEND NEW DATA
guidata(src,g);
end