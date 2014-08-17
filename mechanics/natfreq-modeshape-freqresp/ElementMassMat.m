% This routine computes the element mass matrix for 
% a two-dimensional frame element
%
% function [MM] = ElementMassMat(RHO, X1, Y1, X2, Y2);
%
% INPUTS
% ------
% RHO    - Mass per unit length 
% X1, Y1 - X and Y coordinates of node 1
% X2, Y2 - X and Y coordinates of node 2
%
% OUTPUT
% ------
% MM     - 6 x 6 element mass matrix
% 
      function [MM] = ElementMassMat(RHO, X1, Y1, X2, Y2);

      ELL = sqrt((X2-X1)^2 + (Y2-Y1)^2);

      C=(X2-X1)/ELL;
      S=(Y2-Y1)/ELL;

      MM = zeros(6,6);

% Consistent Mass matrix
     
      MM(1,1)=140.*C*C+156.*S*S;
      MM(1,2)=-16.*C*S;
      MM(1,3)=-22.*ELL*S;
      MM(1,4)=70.*C*C+54.*S*S;
      MM(1,5)=-MM(1,2);
      MM(1,6)=13.*ELL*S;

      MM(2,1)=MM(1,2);
      MM(2,2)=156.*C*C+140.*S*S;
      MM(2,3)=22.*ELL*C;
      MM(2,4)=-MM(1,2);
      MM(2,5)=54.*C*C+70.*S*S;
      MM(2,6)=-13.*ELL*C;
      
      MM(3,1)=MM(1,3);
      MM(3,2)=MM(2,3);
      MM(3,3)=4.*ELL*ELL;
      MM(3,4)=-MM(1,6);
      MM(3,5)=-MM(2,6);
      MM(3,6)=-3.*ELL*ELL;

      MM(4,1)=MM(1,4);
      MM(4,2)=MM(2,4);
      MM(4,3)=MM(3,4);
      MM(4,4)=MM(1,1);
      MM(4,5)=MM(1,2);
      MM(4,6)=-MM(1,3);

      MM(5,1)=MM(1,5);
      MM(5,2)=MM(2,5);
      MM(5,3)=MM(3,5);
      MM(5,4)=MM(4,5);
      MM(5,5)=MM(2,2);
      MM(5,6)=-MM(2,3);

      MM(6,1)=MM(1,6);
      MM(6,2)=MM(2,6);
      MM(6,3)=MM(3,6);
      MM(6,4)=MM(4,6);
      MM(6,5)=MM(5,6);
      MM(6,6)=MM(3,3);

      FAC=RHO*ELL/420.;

      MM = FAC*MM ;
      

