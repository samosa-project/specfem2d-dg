function deriv_u = tridag(BB,DD,AA,CC,NN)
  
%   integer :: NN
%   double precision, dimension(NN) :: AA, BB, CC, DD
% 
%   ! local
%   integer :: I, J
% 
% ! this subroutine performs tridiagonal elimination
% ! BB=coefficient behind diagonal
% ! DD=coefficient on diagonal
% ! AA=coefficient ahead of diagonal 
% ! CC=R.H.S vector; contains the solution on kernel_gammareturn
% ! NN=Number of equations
% !
% ! Note=DD and CC are overwritten by the subroutine

%     CC = zeros(NN,1);
        for I=2:NN
                DD(I)=DD(I)-AA(I-1)*(BB(I)/DD(I-1));
                CC(I)=CC(I)-CC(I-1)*(BB(I)/DD(I-1));
        end

        CC(NN) = CC(NN)/DD(NN);
        for I=2:NN
                J=NN-I+1;
                CC(J)=(CC(J)-AA(J)*CC(J+1))/DD(J);
        end
        
        deriv_u = CC;
      
        return
         
end