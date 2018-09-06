function [X] = scan_database_file(FILE)
  fclose('all'); f=fopen(FILE,'r');
  stop=0;
  skip=1;
  X=[];
  while(stop~=1)
    line = fgetl(f);
  %   line
    if(skip==0)
      if(not(isempty(regexp(line,"numat ngnod nspec pointsdisp plot_lowerleft_corner_only", "once"))))
        break;
      else
        a=str2num(line);
        X=[X;a(2:3)];
      end
    else
      if(not(isempty(regexp(line,"coorg","once"))))
        skip=0;
        continue
      end
    end
  end
  fclose('all');
end

