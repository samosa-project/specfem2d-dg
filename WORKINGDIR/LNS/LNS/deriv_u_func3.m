function deriv_u = deriv_u_func3(delx, grid_points, Flux, forward)
%Function written by Hales Swift, December 31, 2012
%This is an implementation of a 6th order non-dissipative compact difference as shown in: 
%S. K. Lele. Compact finite difference schemes with spectral-like resolution.
%Journal of Computational Physics, 103(1):16–42, November 1992.
%While functional, this implementation is designed for demonstration 
%and so clarity was pursued at the expense of a more efficient implementation

%Using a tridiagonal solver should allow for a much more efficient solution to the problem. 
%What we have is a system of the form Ax'=Bx=b
N=grid_points;
inarray = reshape(Flux,[grid_points,1]);
delta_scale=delx;

[m,n]=size(inarray);
if n>m
    inarray=inarray';
end
%% Initialize variable
A=zeros(N,N);
B=zeros(N,N);
%% Populate the A-matrix
%left boundary
A(1,1)=1;A(1,2)=2;
%2nd point from left
A(2,1)=1/4; A(2,2)=1; A(2,3)=1/4;
%All interior points
for kk=3:N-2
    A(kk,kk-1)=1/3; A(kk,kk)=1; A(kk,kk+1)=1/3;
end
%2nd point from right
A(N-1,N-2)=1/4; A(N-1,N-1)=1; A(N-1,N)=1/4;
%right boundary
A(N,N-1)=2; A(N,N)=1;

%% Populate the B-matrix
alf=7/9/delta_scale;
bet=1/36/delta_scale;
%left boundary
B(1,1)=-5/2/delta_scale; B(1,2)=2/delta_scale; B(1,3)=1/2/delta_scale;
%2nd point from left
B(2,1)=-3/4/delta_scale; B(2,3)=3/4/delta_scale;
%interior points
for kk=3:N-2
    if(forward)
    B(kk,kk-2)=-bet; B(kk,kk-1)=-alf; B(kk,kk+1)=alf; B(kk,kk+2)=bet;
    else
    B(kk,kk+2)=bet; B(kk,kk+1)=alf; B(kk,kk-1)=-alf; B(kk,kk-2)=-bet;
    end
end
%2nd point from right
B(N-1,N-2)=-3/4/delta_scale; B(N-1,N)=3/4/delta_scale;
%right boundary
B(N,N)=5/2/delta_scale; B(N,N-1)=-2/delta_scale; B(N,N-2)=-1/2/delta_scale;

b=B*inarray;
% deriv_u=linsolve(A,b);
deriv_u=A\b;

end