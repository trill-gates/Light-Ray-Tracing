% Date: February 2013
% Title: Basic ray tracing programme
% Author: Deroshan Padotan

%This programme simulates a single light ray passing through a material
%which contains a series of interfaces between two different media. The
%slope of the interfaces are randomized over 180 degrees about the incident
%ray and the main loop is terminated when the ray exits the material.
%NB: In this version of the ray trace programme the only variable that is
%stored at each step is the intercept coordinates of the ray. All other
%variables are overwritten at each iteration. 

clear;

TIR = 0;

in = [1, 0, 1]; %Direction vector for incident ray
ix = [0; 0; 0]; %Position vector for ray on first plane
n = [0, 0, -1]; %Unit normal of first plane
ray_energy = 1; %Preset energy of the incident ray

%The medium variable "med_var" keeps track of which medium the ray is
%currently travelling through. When med_var = 1, the ray is currently in
%medium 1 whereas when med_var = -1 the ray is currently travelling in
%medium 2.
med_var = 1;

%The direction variable "dir_var" stores the z-direction in which the ray 
%is travelling. When dir_var = 1, the ray being followed is travelling in
%the positive z-direction otherwise it is travelling in the negative
%z-direction. 
dir_var = 1;

%The lambda parameter controls the slope of the exponential
%distribution from which the distances between successive interfaces
%are randomly chosen
lambda = 1;

r_z = [0, 0, 1]; %General ray along the z-axis

%Refractive indices of the two media
n1 = 1;
n2 = 1.1;

%Calculating the critical angle
theta_c = asind(n1/n2);

%Normalising the incident ray direction vector
in_unit = in/norm(in);

i = 0;

while ix(3)>=0
    
    i = i+1;
    
    ix_m(i,:) = ix;
    
    if i == 1
        
        r_z = in_unit;
        
    else
        
        r_z = [0,0,dir_var];
        
    end
    
    th_m1 = acosd(dot(r_z, -sign(r_z(3))*n)); %Calculating the incident angle
    
    %Outer if statement checks which media the refracted ray is travelling
    %through wih each iteration.
    %The inner if statements check for total internal reflection by checking
    %whether the refracted ray returned is complex.
    %NB: Total internal reflection can only occur when the light ray is
    %travelling from a more dense medium to a less dense medium
    
    rec_v = 0;
    
    if med_var == -1
        
        rf = refract(r_z, sign(r_z(3))*n, n2, n1);
        
        if(isreal(rf) == 0 && th_m1>theta_c)
            
            rf = [0,0,0];
            
            TIR = i;
            
            rec_v = 1;
            
        end
        
        
    else
        
        rf = refract(r_z, sign(r_z(3))*n, n1, n2);
        
    end
    
    rl = reflect(r_z, sign(r_z(3))*n);
    
    if rec_v == 1
        
        th_m2 = -acosd(dot(rf,-sign(rf(3))*n));
        
    else
        
        th_m2 = acosd(dot(rf,-sign(rf(3))*n)); %Calculating the refraction angle
        
    end
    
    %Calculating the reflection coefficient at each interface
    rs = ((n1*cosd(th_m1) - n2*cosd(th_m2))/(n1*cosd(th_m1) + n2*cosd(th_m2)))^2;
    rp = ((n1*cosd(th_m2) - n2*cosd(th_m1))/(n1*cosd(th_m2) + n2*cosd(th_m1)))^2;
    R_coef = (rp+rs)/2;
    
    ray_energy = ray_energy*R_coef;
    
    %Randomising the distance between successive interfaces over an exponential
    %distribution
    sy = rand(1);
    
    s = (-log(sy))/lambda;
    
    if i == 1
        
        %Calculating the rotation matrix
        R_mat = vecRotMat(rf,r_z);
        
        ix = ix + s*rf';
        
        rf_m(1,:) = rf;
        
        med_var = -med_var;
        
    else
        
        
        %By comparing a randomly generated number to the reflection coefficent at
        %each interface a decision is made as to whether the reflected or the
        %refracted ray willl be followed.
        if ((rand(1) <= R_coef) || (rec_v == 1))
            
            rl_r = R_mat'*rl';
            
            rl = rl_r';
            
            %Using the reflected vector
            ix = ix + s*rl';
            
            R_mat = vecRotMat(rl,r_z);
            
            dir_var = sign(rl(3));
            
        else
            med_var = -med_var;
            
            rf_r = R_mat'*rf';
            
            rf = rf_r';
            
            %Using the refracted vector
            ix = ix + s*rf';
            
            R_mat = vecRotMat(rf,r_z);
            
            dir_var = sign(rf(3));
            
        end
        
    end
    
    %Using the George Marsaglia method to generate random normal vectors
    %distributed evenly over the surface of a half-sphere
    S_n = 2;
    
    while S_n>1
        
        n_3 = -1 + rand(1,1);
        n_1 = -1 + 2*rand(1,1);
        
        S_n = n_3^2 + n_1^2;
        n(2) = 1 - 2*S_n;
        
    end
    
    n(3) = (2*n_3)*sqrt((1-S_n));
    n(1) = (2*n_1)*sqrt((1-S_n));
    n(2) = 1 - 2*S_n;
    
    
end

str3 = sprintf('The ray exited the material after %d interfaces', i);
disp(str3);

ix_m(i+1,:) = [0,0,0];

ix_m2 = zeros(i+1,3);

%NB: ix_m2 is a matrix which stores the coordinates at which the light rays
%intercept each plane. The columns have been re-ordered to comply with the
%altered optic coordinate system.
ix_m2(:,1) = ix_m(:,3);
ix_m2(:,2) = ix_m(:,1);
ix_m2(:,3) = ix_m(:,2);

ix2(1) = ix(3);
ix2(2) = ix(1);
ix2(3) = ix(2);


%Outer if statement checks whether TIR has occurred at interfcae 2 in which
%case no plot is generated
if any(ix_m(2,:))
    figure();
    for k = 1:i
        if any(ix_m2(k+1,:))
            
            %Displaying each of the rays between interfaces
            vectarrow(ix_m2(k,:),ix_m2(k+1,:));
            
            hold on;
            
        else
            
            hold on;
            
            vectarrow(ix_m2(k,:),ix2);
            
        end
        
    end
    hold on;
    scatter3(0,0,0,'filled','g');
    hold on;
    scatter3(ix(3),ix(1),ix(2),'filled','r');
    
    xlabel('z')
    ylabel('x')
    zlabel('y')
    
    
end

%Printing all the data to a text file
fid = fopen('R13Data.txt','w');
fprintf(fid,'Data for test run of Ray_Trace12\r\n\r\nFor this run the ray exited the material after %d interfaces\r\n\r\n',i);
fprintf(fid,'Preset Parameters:\r\n\r\n');
fprintf(fid,'n1 = %.2f\n n2 = %.2f\r\n',n1,n2);
fprintf(fid, 'lambda = %.2f\r\n\r\n',lambda);

fprintf(fid, 'Intercept Points:\r\n\r\nix_m =\r\n\r\n');

for j = 1:(i)
    fprintf(fid, '%.5f  %.5f  %.5f\r\n',ix_m(j,:));
end

fclose(fid);

