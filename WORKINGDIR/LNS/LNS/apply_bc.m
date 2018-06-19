% function [q21_out,q22_out,q1_out,q3_out] = apply_bc(NX,NY,q21,q22,q1,q3,gamma_atmos,c,is_backward,t)
function [q1_out,q21_out,q22_out,q3_out] = apply_bc(q1,q21,q22,q3,gamma_atmos,c,t)

%     %% Reshape fields if needed
%     reshape_bool = false;
%     if(size(vx,3) > 1)
%         reshape_bool = true;
%         vx = reshape(vx,[NX,NY]);
%         vy = reshape(vy,[NX,NY]);
%         rho1 = reshape(rho1,[NX,NY]);
%         p = reshape(p,[NX,NY]);
%     end
  
  % Remove singleton dimensions.
  reshape_bool = false;
  if(size(q1,3) > 1)
    reshape_bool = true;
    q1=squeeze(q1);
    q21=squeeze(q21);
    q22=squeeze(q22);
    q3=squeeze(q3);
  end
  
  NX=size(q1,1);
  NY=size(q1,2);
  
  is_backward=false;
  
  LAPO=1; % Duration of beginning apodisation.
  apot0 = 0.5 .* (1 + erf((t - 0.5 * LAPO) / (0.25 * LAPO))); % Apodisation over LAPO at beginning.
  if(t==0)
    apot0=0;
  end

  % Boundaries at y=y_min and y=y_max (horizontal).
  for i = 1:NX
    for j = [1 NY]
%       % Periodic.
%       if(j==NY)
%         % Map j=NY to j=1.
%         q1(i,j)=q1(i,NY-j+1);
%         q21(i,j)=q21(i,NY-j+1);
%         q22(i,j)=q22(i,NY-j+1);
%         q3(i,j)=q3(i,NY-j+1);
%       end
      
      % Walls.
      coef2 = 0;
      if(j == 1)
        coef1 = 1;
        if(i==1)
          coef2 = 1;
        end
      else
        coef1 = -1;
        if(i==NX)
          coef2 = -1;
        end
      end
      if(is_backward)
        q1(i,j) = 0;%-rho1(i+coef2,j+coef);
        q21(i,j) = 0;%-vx(i+coef2,j+coef);
        q22(i,j) = 0;%-vy(i+coef2,j+coef);
        q3(i,j) = q3(i+coef2,j+coef1);
      else
        q1(i,j) = q1(i+coef2,j+coef1); % \partial_yq_1=0.
%         q21(i,j) = -q21(i+coef2,j+coef1); % ??
        q21(i,j) = q21(i+coef2,j+coef1); % \partial_yq_21=0.
        q22(i,j) = -q22(i+coef2,j+coef1); % q22=0 (wall).
%         pg = (1/22)*(17*q3(i+coef2,j+coef)-9*q3(i+coef2*2,j+coef*2)-5*q3(i+coef2*3,j+coef*3)-q3(i+coef2*4,j+coef*4));
        q3(i,j) = q3(i+coef2,j+coef1); % \partial_yq_3=0.
      end
      
    end
  end

  % Boundaries at x=x_min and x=x_max (vertical).
  for j = 1:NY
    for i = [1 NX]
%       % Periodic.
%       if(i==NX)
%         % Map i=NX to i=1.
%         q1(i,j)=q1(NX-i+1,j);
%         q21(i,j)=q21(NX-i+1,j);
%         q22(i,j)=q22(NX-i+1,j);
%         q3(i,j)=q3(NX-i+1,j);
%       end
      
      % Walls.
      coef2 = 0;
      if(i == 1)
         coef1 = 1;
         if(j==1)
          coef2 = 1;
         end
      else
         coef1 = -1;
         if(j==NY)
          coef2 = -1;
         end
      end
      if(is_backward)
        q1(i,j) = 0;%-rho1(i+coef,j+coef2);
        q21(i,j) = 0;%-vx(i+coef,j+coef2);
        q22(i,j) = 0;%-vy(i+coef,j+coef2);
        q3(i,j) = q3(i+coef1,j+coef2);
      else
        q1(i,j) = q1(i+coef1,j+coef2); % \partial_xq_1=0. %pg/(gamma_atmos(i+coef,j+coef2)*c^2);%rho1(i+coef,j+coef2);
        
%         q21(i,j) = -q21(i+coef1,j+coef2); % q21=0 (wall).
        if(i==1)
          q21(i,j) = apot0*sin(t*2*pi/0.5);%* 0.5*(1 + erf((t - 0.5 * 0.1) / (0.25 * 0.1))); % Oscillation.
        else
          q21(i,j) = -q21(i+coef1,j+coef2); % q21=0 (wall).
        end
        
        q22(i,j) = q22(i+coef1,j+coef2); % \partial_xq22=0.
%         q21(i,j) = -q21(i+coef1,j+coef2); % Wall for v_x (??).
%         q22(i,j) = sin(2*pi*t/0.3);%-vy(i+coef,j+coef2); % Oscillating v_y (??).
%         pg = (1/22)*(17*q3(i+coef,j+coef2)-9*q3(i+coef*2,j+coef2*2)-5*q3(i+coef*3,j+coef2*3)-q3(i+coef*4,j+coef2*4));
        q3(i,j) = q3(i+coef1,j+coef2); % \partial_xq3=0.
      end
      
    end
  end
  
  % Set output variables, eventually adding the singleton dimension for cohesion.
  q1_out = q1;
  q21_out = q21;
  q22_out = q22;
  q3_out = q3;
  if(reshape_bool)
    q1_out = permute(q1_out,[3,1,2]);
    q21_out = permute(q21_out,[3,1,2]);
    q22_out = permute(q22_out,[3,1,2]);
    q3_out = permute(q3_out,[3,1,2]);
  end
%     %% Reshape fields if needed
%     if(reshape_bool)
%         q21_out = reshape(q21,[1,NX,NY]);
%         q22_out = reshape(q22,[1,NX,NY]);
%         q1_out = reshape(q1,[1,NX,NY]);
%         q3_out = reshape(q3,[1,NX,NY]);
%     else
%         q21_out   = q21;
%         q22_out   = q22;
%         q1_out = q1;
%         q3_out    = q3;
%     end
        
end