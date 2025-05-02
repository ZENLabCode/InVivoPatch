function fun_close(src,~)
g = guidata(src);
assignin('caller','g',g)
delete(findobj(src,'type','figure'));
end

