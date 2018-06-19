function deriv_u = deriv_u_func2(delx, grid_points, Flux,forward)

%Initialization & Value Declaration
df_r=zeros(grid_points,1);
df_n=zeros(grid_points);                            %declaring numerical solution
%Co-efficient values for Boundary points
% eq. (4.1.3) Lele 1992
alpha_b=0;                                          %boundry alpha will be zero to make it explicit
a_b=-((11+2*alpha_b)/6);
b_b=(6-alpha_b)/2;
c_b=(2*alpha_b-3)/2;
d_b=(2-alpha_b)/6;
%Co-efficient values for interior points 
% eq. (2.1.6) Lele 1992
alpha_i=1/4;                                        %input('Input the value of Alpha');
beta_i=0;
a_i=(2/3)*(alpha_i+2);
b_i=(1/3)*((4*alpha_i)-1);
c_i=0;

%---------------------------------------------------------%
%               Right Side Node Calculation               %
%---------------------------------------------------------%

%Formulation for 1st derivative 
% eq. (2.1) Lele 1992
df_r=(a_i/(2*delx))*((diag(ones(grid_points-1,1),1)-(diag(ones(grid_points-1,1),-1))))+(b_i/(4*delx))*((diag(ones(grid_points-2,1),2)-(diag(ones(grid_points-2,1),-2))))+(c_i/(6*delx))*((diag(ones(grid_points-3,1),3)-(diag(ones(grid_points-3,1),-3))));

%Changing Boundary Values as per calculation 
df_r(1,1)=(a_b)/delx;
df_r(1,2)=(b_b)/delx;
df_r(1,3)=(c_b)/delx;
df_r(1,4)=(d_b)/delx;

df_r(grid_points,grid_points)=-(a_b)/delx;
df_r(grid_points,grid_points-1)=-(b_b)/delx;
df_r(grid_points,grid_points-2)=-(c_b)/delx;
df_r(grid_points,grid_points-3)=-(d_b)/delx;

%---------------------------------------------------------%
%               Left Side Node Calculation                %
%---------------------------------------------------------%
% eq. (2.1) Lele 1992
df_l=diag(ones(grid_points,1))+alpha_i*(diag(ones(grid_points-1,1),-1)+diag(ones(grid_points-1,1),+1))+beta_i*(diag(ones(grid_points-2,1),-2)+diag(ones(grid_points-2,1),2));

%Changing Boundary Values & making it explicit
df_l(1,1)=1;
df_l(1,2:grid_points)=0;
df_l(grid_points,grid_points)=1;
df_l(grid_points,1:grid_points-1)=0;

% size((df_l\df_r))
% size(Flux)

deriv_u=(df_l\df_r)*reshape(Flux,[grid_points,1]);

end