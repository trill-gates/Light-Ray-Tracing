function r = refract(r, n, n1, n2)

% REFRACT Implement Snell's Law. 
% 
% REFRACT(r, n, n1, n2). Here the input slope is r, n is the normal and 
% (n1, n2) are indices of refraction. Complex return indicates total 
% internal reflection. 

% 15/12/08     Allow vector of refractive indexes
% 9/3/05       reflection equation fixed 
% 9/3/05       refraction from back surfaces doesnt appear to work 
% 10/3/05      added a delta sign convention 
        
if	((n1 == n2)),
    % nothing to do.
    return 
end

r   = vunit(r);
n   = vunit(n);
in  = vdot(r,n);

% refraction         
d   = 1 - ((n1/n2).^2).*(1 - in.^2);

% 10/3/05 sign convention 
delta   = -sign(n(:,3).*r(:,3));
n0      = (n1/n2).*in + delta.*sqrt(d);
n0      = repmat(n0, 1, 3);
r       = -n0.*n + r*(n1/n2);
return       
