%[Vcoord1,Vcoord2]=PixtoCoord(Cam,Vpix1,Vpix2)
function [Vcoord1,Vcoord2]=PixtoCoordGPP(Vpix1,Vpix2,z)
%Modif le 20 Oct ajout de z en input (coordonnée verticale du plan d'eau

clear Vcoord1 Vcoord2
load 'RectGPP';

%Version avec altitude
[X_I] = rectifyGPP(double([Vpix1; Vpix2]),Rckk,Tckk,fc,cc,kc,alpha_c,z);

Vcoord1=X_I(1,:);
Vcoord2=X_I(2,:);


function [X_r] = rectifyGPP(aaa,Rckk,Tckk,fc,cc,kc,alpha_c,Z)

% aaa contains the pixel coordinates (origin top left, horizontal in the
% first row vertical in the second row).
%
% Z is the level you want to recitfy to (can be 
% different for each point, in which case length(Z) must be the same as X_kk))

[n,m]=size(aaa);
if m==2&n~=2
    x_kk=aaa';
end
[n,m]=size(aaa);
if length(Z)==1
    Z=Z.*ones(1,m);
end

%undistort the image coordinates and scale
xn = normalizeGPP(aaa,fc,cc,kc,alpha_c); 


xn=[xn;ones(1,length(xn))];

%rotate the normalised coordinates
R=Rckk'*xn;

%rotate the camera location
T=Rckk'*Tckk*ones(1,length(xn));

%figure out what the scale factor is
z=(Z+T(3,:))./R(3,:);

%apply the scale factor to rotated coordinates, and correct for camera
%position.
X_r=[z.*R(1,:)-T(1,:);z.*R(2,:)-T(2,:);Z];


function [xn] = normalizeGPP(x_kk,fc,cc,kc,alpha_c)

%normalize
%
%[xn] = normalize(x_kk,fc,cc,kc,alpha_c)
%
%Computes the normalized coordinates xn given the pixel coordinates x_kk
%and the intrinsic camera parameters fc, cc and kc.
%
%INPUT: x_kk: Feature locations on the images
%       fc: Camera focal length
%       cc: Principal point coordinates
%       kc: Distortion coefficients
%       alpha_c: Skew coefficient
%
%OUTPUT: xn: Normalized feature locations on the image plane (a 2XN matrix)
%
%Important functions called within that program:
%
%comp_distortion_oulu: undistort pixel coordinates.

if nargin < 5,
   alpha_c = 0;
   if nargin < 4;
      kc = [0;0;0;0;0];
      if nargin < 3;
         cc = [0;0];
         if nargin < 2,
            fc = [1;1];
         end;
      end;
   end;
end;


% First: Subtract principal point, and divide by the focal length:
x_distort = [(x_kk(1,:) - cc(1))/fc(1);(x_kk(2,:) - cc(2))/fc(2)];

% Second: undo skew
x_distort(1,:) = x_distort(1,:) - alpha_c * x_distort(2,:);

if norm(kc) ~= 0,
	% Third: Compensate for lens distortion:
	xn = comp_distortion_ouluGPP(x_distort,kc);
else
   xn = x_distort;
end


function [x] = comp_distortion_ouluGPP(xd,k)

%comp_distortion_oulu.m
%
%[x] = comp_distortion_oulu(xd,k)
%
%Compensates for radial and tangential distortion. Model From Oulu university.
%For more informatino about the distortion model, check the forward projection mapping function:
%project_points.m
%
%INPUT: xd: distorted (normalized) point coordinates in the image plane (2xN matrix)
%       k: Distortion coefficients (radial and tangential) (4x1 vector)
%
%OUTPUT: x: undistorted (normalized) point coordinates in the image plane (2xN matrix)
%
%Method: Iterative method for compensation.
%
%NOTE: This compensation has to be done after the subtraction
%      of the principal point, and division by the focal length.


if length(k) == 1,
    
    [x] = comp_distortionGPP(xd,k);
    
else
    
    k1 = k(1);
    k2 = k(2);
    k3 = k(5);
    p1 = k(3);
    p2 = k(4);
    
    x = xd; 				% initial guess
    
    for kk=1:20,
        
        r_2 = sum(x.^2);
        k_radial =  1 + k1 * r_2 + k2 * r_2.^2 + k3 * r_2.^3;
        delta_x = [2*p1*x(1,:).*x(2,:) + p2*(r_2 + 2*x(1,:).^2);
        p1 * (r_2 + 2*x(2,:).^2)+2*p2*x(1,:).*x(2,:)];
        x = (xd - delta_x)./(ones(2,1)*k_radial);
            
    end
    
end


function [x_comp]  = comp_distortionGPP(x_dist,k2)

%       [x_comp] = comp_distortion(x_dist,k2);
%       
%       compensates the radial distortion of the camera
%       on the image plane.
%       
%       x_dist : the image points got without considering the
%                radial distortion.
%       x : The image plane points after correction for the distortion
%       
%       x and x_dist are 2xN arrays
%
%       NOTE : This compensation has to be done after the substraction
%              of the center of projection, and division by the focal
%              length.
%       
%       (do it up to a second order approximation)


[two,N] = size(x_dist);


if (two ~= 2 ), 
    error('ERROR : The dimension of the points should be 2xN');
end;


if length(k2) > 1,
    
    [x_comp]  = comp_distortion_ouluGPP(x_dist,k2);
    
else
    
    radius_2= x_dist(1,:).^2 + x_dist(2,:).^2;
    radial_distortion = 1 + ones(2,1)*(k2 * radius_2);
    radius_2_comp = (x_dist(1,:).^2 + x_dist(2,:).^2) ./ radial_distortion(1,:);
    radial_distortion = 1 + ones(2,1)*(k2 * radius_2_comp);
    x_comp = x_dist ./ radial_distortion;
    
end
