function prettyAxes(f)
  children=f.Children;
  for i=1:numel(children)
    child=children(i);
    if(strcmp(child.Type,'axes'))
      axes(child);
%       set(gca, 'Color','k');
%       set(gca, 'GridColor','white');
      set(gca, 'TickLabelInterpreter', 'latex');
      set(gca, 'TickDir','both');
      set(gca, 'TickLabelInterpreter', 'latex');
      grid on;
      box on;
    elseif(strcmp(children(i).Type,'legend'))
      set(child,'fontsize', 16);
%       set(child, 'Color',[1,1,1]*0.25);
%       set(child, 'textColor','w');
    end
  end
end