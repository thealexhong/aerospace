  function w = NCweights(m)
% w = NCweights(m)
%
% w is a column m-vector consisting of the weights for the m-point Newton-Cotes rule.
% m is an integer that satisfies 2 <= m <= 11.
 
if      m==2,   w=[1 1]'/2; 
elseif  m==3,   w=[1 4 1]'/6; 
elseif  m==4,   w=[1 3 3 1]'/8; 
elseif  m==5,   w=[7 32 12 32 7]'/90; 
elseif  m==6,   w=[19 75 50 50 75 19]'/288; 
elseif  m==7,   w=[41 216 27 272 27 216 41]'/840; 
elseif  m==8,   w=[751 3577 1323 2989 2989 1323 3577 751]'/17280; 
elseif  m==9,   w=[989 5888 -928 10496 -4540 10496 -928 5888 989]'/28350; 
elseif  m==10,  w=[2857 15741 1080 19344 5778 5778 19344 1080 15741 2857]'/89600; 
else            w=[16067 106300 -48525 272400 -260550 427368 -260550 272400 -48525 106300 16067]'/598752; 
end