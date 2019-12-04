% Copyright (c) 2016, Sathyanarayan Rao
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are
% met:
% 
%     * Redistributions of source code must retain the above copyright
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright
%       notice, this list of conditions and the following disclaimer in
%       the documentation and/or other materials provided with the distribution
% 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
% POSSIBILITY OF SUCH DAMAGE.

%------------------------------------------------------------------------%
%------ Gmsh to Matlab script: Import mesh to matlab---------------------%
%------------------------------------------------------------------------%

clc
close all
clear 

%-----------------------------------------------------------------------%
% dlmread(filename,delimiter,[R1 C1 R2 C2]) reads only the range 
% bounded by row offsets R1 and R2 and column offsets C1 and C2.
%-----------------------------------------------------------------------%
file    =  input('your GMsh file in .msh format','s');

% no of nodes is mentioned in 5th row and first column

N_n      = dlmread(file,'',[5-1 1-1 5-1 1-1]);
N_e      = dlmread(file,'',[7+N_n 0 7+N_n 0]);

node_id     = dlmread(file,'',[5 0 4+N_n 0]);
nodes       = dlmread(file,'',[5 1 4+N_n 3]);
elements    = dlmread(file,'',[8+N_n 0 7+N_n+N_e 7]);

%------- 2D Geometry

two_d_nodes = nodes(:,1:2);
elem_type   = elements(:,2);

%--- find the starting indices of 2D elements
two_ind = 1;
for i = 1:N_e
    if(elem_type(i) ~= 2)
        two_ind = two_ind+1;
    end
end
%----------------------------------------------

two_d_elements(1:N_e-two_ind,1:3) = 0;
k = 1;
for i = two_ind:N_e
   two_d_elements(k,1:3) = elements(i,6:8);
   k = k+1;
end

%---- visualize in matlab ---------------------

figure(1)
triplot(two_d_elements,two_d_nodes(:,1),two_d_nodes(:,2))
xlabel('X','fontsize',14)
ylabel('Y','fontsize',14)
title('GMsh to MATLAB import','fontsize',14)
fh = figure(1);
set(fh, 'color', 'white'); 

%-------------------------------------------------------------------------





