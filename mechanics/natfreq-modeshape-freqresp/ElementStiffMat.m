% This routine computes the element stiffness matrix for 
% a two-dimensional frame element
%
% function [KM] = ElementStiffMat(EA, EI, X1, Y1, X2, Y2);
% 
% INPUTS
% ------
% EA      - axial rigidity
% EI      - flexural rigidity
% X1, Y1  - X and Y coordinates of node 1
% X2, Y2  - X and Y coordinates of node 2
%
% OUTPUT
% ------
% KM      - 6 x 6 element stiffness matrix 
%
     function [KM] = ElementStiffMat(EA, EI, X1, Y1, X2, Y2);

      L = sqrt((X2-X1)^2 + (Y2-Y1)^2);

      C=(X2-X1)/L;
      S=(Y2-Y1)/L;

      E1=EA/L;
      E2=12.*EI/(L*L*L);
      E3=EI/L;
      E4=6.*EI/(L*L);

      KM = zeros(6,6);

% Define elements of K matrix

      KM(1,1)=C*C*E1+S*S*E2;
      KM(4,4)=KM(1,1);
      KM(1,2)=S*C*(E1-E2);
      KM(2,1)=KM(1,2);
      KM(4,5)=KM(1,2);
      KM(5,4)=KM(4,5);
      KM(1,3)=-S*E4;
      KM(3,1)=KM(1,3);
      KM(1,6)=KM(1,3);
      KM(6,1)=KM(1,6);
      KM(3,4)=S*E4;
      KM(4,3)=KM(3,4);
      KM(4,6)=KM(3,4);
      KM(6,4)=KM(4,6);
      KM(1,4)=-KM(1,1);
      KM(4,1)=KM(1,4);
      KM(1,5)=S*C*(-E1+E2);
      KM(5,1)=KM(1,5);
      KM(2,4)=KM(1,5);
      KM(4,2)=KM(2,4);
      KM(2,2)=S*S*E1+C*C*E2;
      KM(5,5)=KM(2,2);
      KM(2,5)=-KM(2,2);
      KM(5,2)=KM(2,5);
      KM(2,3)=C*E4;
      KM(3,2)=KM(2,3);
      KM(2,6)=KM(2,3);
      KM(6,2)=KM(2,6);
      KM(3,3)=4.*E3;
      KM(6,6)=KM(3,3);
      KM(3,5)=-C*E4;
      KM(5,3)=KM(3,5);
      KM(5,6)=KM(3,5);
      KM(6,5)=KM(5,6);
      KM(3,6)=2.*E3;
      KM(6,3)=KM(3,6);

      

