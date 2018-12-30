function y = vnorm(y)

% VNORM determine the (2-)norm of the [nxd] matrix 
% (n rows of d element vectors)
%

y = sqrt(sum((y.^2)'))'; 

return