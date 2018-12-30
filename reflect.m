function r = reflect(r, n)

% REFLECT Implement the Law of Reflection.  
% 
% REFLECT(nSlope, nNorm). Here nSlope is the slope of the incident ray,  
% n is the normal at the surface point.  This function will return the 
% reflected ray.
%
% 21/11/06 adapted from SNELL.m
% 

if (any(size(r) ~= size(n)))
    [r, n] = mesh3vec(r, n);
end 

r  = vunit(r); n = vunit(n); n0 = vdot(r, n);
n0 = repmat(n0, 1, 3);
r = -2*n0.*n + r;
return

