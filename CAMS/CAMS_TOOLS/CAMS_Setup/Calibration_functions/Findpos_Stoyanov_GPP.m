% Findpos_Stoyanov_GPP.m

function [err]=Findpos_Stoyanov_GPP(focale,dist1,dist2,ImageGPP,PointsControl)

[xx yy cc]=size(imread(ImageGPP));


%make sure it can find the Stoyanov software
%addpath 'D:\UEAcamera\Stoyanov\TOOLBOX_calib\'
%**************************************************************************
%
% Change next two line so Matlab can find files on your path
%
%**************************************************************************
% addpath 'O:\UBF07201\RawData\Rectification\April2007\'
% addpath 'O:\UBF07201\RawData\Rectification\April2007\TOOLBOX_calib'
%addpath 'D:\Raphael\THESE\NZ\codes\RectifyImages'
%if using stoyanov software for lens calibration stoy=1, otherwise need to
%convert from Heikkela results
%**************************************************************************
%
% Change the line below for the camera you want to work on
%
%**************************************************************************

%Which camera will you try?
%**************************************************************************
%
% Change the paths to reflect where you put the camera files, images and
% data files
%
%*************************************************************************

%labels for the plotting
lab=['01';'02';'03';'04';'05';'06';'07';'08';'09';'10';'11';'12';'13';'14';'15';'16';'17';'18';'19';'20';'21';'22';'23';'24';'25';'26';'27';'28';'29';'30'];

%**************************************************************************
%
% At this stage no lens distortion info used, this routine seems powerful
% enough to work without lens distortion

    
	par=[0.9992	%read 'readme'
        
	 focale%2.5	%focal length 1.9
%	  1091.13	%x-center
% 	  757.09	%y-center
%	 -8.703966e-03
%     7.905523e-05	
%	 0.007173
%     -0.020489];
	  round(yy/2)	%x-center
 	  round(xx/2)	%y-center
	 0.0
     	 0.0	
	 0.0
     	 0.0];
 
    [fc,cc,kc,alpha_c]=convert_heikkila_stoyanovGPP(par);
    k1 = dist1; k2 = dist2; k3 = 0; p1 = 0; p2 = 0;   % k1 = -0.8; k2 = 0; 
    kc = [k1,k2,p1,p2,k3];
    
	
%enter GCP matrices. It is assumed that all the GCPs are in one 
%datafile 'GCPs.txt' with xpic, ypic, Xground, Yground, 
%Zground all in columns
I=imread(ImageGPP);
GCPs=load([PointsControl]);  %example GCPS
    x_kk=GCPs(:,5:-1:4)';  
    X_kk=GCPs(:,1:3)';


thresh_cond=1000000;   %used for the iteration routine that 
MaxIter=30;

%calculate the camera position
[omckk,Tckk,Rckk,H,x,ex,JJ] = compute_extrinsicGPP(x_kk,X_kk,fc,cc,kc,alpha_c,MaxIter,thresh_cond);
%omckk, tckk, Rckk are the camera rotation and translation matrices, needed to rectify the iamges.  

save(['./Rectfile'])

[Vcoord1,Vcoord2]=PixtoCoordGPP(x_kk(1,:),x_kk(2,:),zeros(1,length(x_kk(1,:))));
%use camera position to rectify points
X_r(1,:)=Vcoord1;
X_r(2,:)=Vcoord2;

% [X_r] = rectifyGPP(x_kk,Rckk,Tckk,fc,cc,kc,alpha_c,X_kk(3,:));

% plotting the image with GCPs
% figure(4)
% image(I)
% hold on
% plot(x_kk(1,:),x_kk(2,:),'mo')
% plot(x_kk(1,:),x_kk(2,:),'kx')
% plot(cc(1),cc(2),'w*')
% [n,m]=size(x_kk);
% for j=1:m
%     text(x_kk(1,j)+10,x_kk(2,j),num2str(j))
% end
% % % 
figure(1);clf;set(gcf,'Color','w')
image(I)
hold on
plot(x_kk(1,:),x_kk(2,:),'*')
plot(x(1,:),x(2,:),'r*')
title(['comparison of actual and calculated image coordinates'])
legend('Actual','Calculated')
axis('ij')
[n,m]=size(x_kk);
for j=1:m
    text(x_kk(1,j)+10,x_kk(2,j),lab(j,:))
    line([x_kk(1,j) x(1,j)],[x_kk(2,j) x(2,j)])
end

% 
figure(2);clf;set(gcf,'Color','w')
plot(X_kk(2,:),X_kk(1,:),'*')
hold on
plot(X_r(2,:),X_r(1,:),'r*')
title(['comparison of actual and calculated ground coordinates (assuming 1st column is Easting)'])
legend('Actual','Calculated')
ylabel(['Easting (m)'])
xlabel(['Northing (m)'])

size(X_kk(1,:))
size(X_r(1,:))

%Calcul de l'erreur sur les points de control
err=mean(sqrt( (X_kk(2,:)-X_r(2,:)).^2  +  (X_kk(1,:)-X_r(1,:)).^2  ))
    

end



function [omckk,Tckk,Rckk,H,x,ex,JJ] = compute_extrinsicGPP(x_kk,X_kk,fc,cc,kc,alpha_c,MaxIter,thresh_cond)

%compute_extrinsic
%
%[omckk,Tckk,Rckk,H,x,ex] = compute_extrinsic(x_kk,X_kk,fc,cc,kc,alpha_c)
%
%Computes the extrinsic parameters attached to a 3D structure X_kk given its projection
%on the image plane x_kk and the intrinsic camera parameters fc, cc and kc.
%Works with planar and non-planar structures.
%
%INPUT: x_kk: Feature locations on the images
%       X_kk: Corresponding grid coordinates
%       fc: Camera focal length
%       cc: Principal point coordinates
%       kc: Distortion coefficients
%       alpha_c: Skew coefficient
%
%OUTPUT: omckk: 3D rotation vector attached to the grid positions in space
%        Tckk: 3D translation vector attached to the grid positions in space
%        Rckk: 3D rotation matrices corresponding to the omc vectors
%        H: Homography between points on the grid and points on the image plane (in pixel)
%           This makes sense only if the planar that is used in planar.
%        x: Reprojections of the points on the image plane
%        ex: Reprojection error: ex = x_kk - x;
%
%Method: Computes the normalized point coordinates, then computes the 3D pose
%
%Important functions called within that program:
%
%normalize: Computes the normalize image point coordinates.
%
%pose3D: Computes the 3D pose of the structure given the normalized image projection.
%
%project_points.m: Computes the 2D image projections of a set of 3D points



if nargin < 8,
   thresh_cond = inf;
end;


if nargin < 7,
   MaxIter = 20;
end;


if nargin < 6,
   alpha_c = 0;
	if nargin < 5,
   	kc = zeros(5,1);
   	if nargin < 4,
      	cc = zeros(2,1);
      	if nargin < 3,
         	fc = ones(2,1);
         	if nargin < 2,
            	error('Need 2D projections and 3D points (in compute_extrinsic.m)');
            	return;
         	end;
      	end;
   	end;
	end;
end;

% Initialization:

[omckk,Tckk,Rckk] = compute_extrinsic_initGPP(x_kk,X_kk,fc,cc,kc,alpha_c);

% Refinement:
[omckk,Tckk,Rckk,JJ] = compute_extrinsic_refineGPP(omckk,Tckk,x_kk,X_kk,fc,cc,kc,alpha_c,MaxIter,thresh_cond);


% computation of the homography (not useful in the end)

H = [Rckk(:,1:2) Tckk];

% Computes the reprojection error in pixels:

x = project_points2GPP(X_kk,omckk,Tckk,fc,cc,kc,alpha_c);

ex = x_kk - x;


% Converts the homography in pixel units:

KK = [fc(1) alpha_c*fc(1) cc(1);0 fc(2) cc(2); 0 0 1];

H = KK*H;


end


function [omckk,Tckk,Rckk] = compute_extrinsic_initGPP(x_kk,X_kk,fc,cc,kc,alpha_c),

%compute_extrinsic
%
%[omckk,Tckk,Rckk] = compute_extrinsic_init(x_kk,X_kk,fc,cc,kc,alpha_c)
%
%Computes the extrinsic parameters attached to a 3D structure X_kk given its projection
%on the image plane x_kk and the intrinsic camera parameters fc, cc and kc.
%Works with planar and non-planar structures.
%
%INPUT: x_kk: Feature locations on the images
%       X_kk: Corresponding grid coordinates
%       fc: Camera focal length
%       cc: Principal point coordinates
%       kc: Distortion coefficients
%       alpha_c: Skew coefficient
%
%OUTPUT: omckk: 3D rotation vector attached to the grid positions in space
%        Tckk: 3D translation vector attached to the grid positions in space
%        Rckk: 3D rotation matrices corresponding to the omc vectors
%
%Method: Computes the normalized point coordinates, then computes the 3D pose
%
%Important functions called within that program:
%
%normalize: Computes the normalize image point coordinates.
%
%pose3D: Computes the 3D pose of the structure given the normalized image projection.
%
%project_points.m: Computes the 2D image projections of a set of 3D points



if nargin < 6,
   alpha_c = 0;
	if nargin < 5,
   	kc = zeros(5,1);
   	if nargin < 4,
      	cc = zeros(2,1);
      	if nargin < 3,
         	fc = ones(2,1);
         	if nargin < 2,
            	error('Need 2D projections and 3D points (in compute_extrinsic.m)');
            	return;
         	end;
      	end;
   	end;
	end;
end;


%keyboard;

% Compute the normalized coordinates:

xn = normalizeGPP(x_kk,fc,cc,kc,alpha_c);



Np = size(xn,2);

%keyboard;

X_mean = mean(X_kk')';

Y = X_kk - (X_mean*ones(1,Np));

YY = Y*Y';

[U,S,V] = svd(YY);

r = S(3,3)/S(2,2);

%keyboard;


if (r < 1e-3)|(Np < 5), %1e-3, %1e-4, %norm(X_kk(3,:)) < eps, % Test of planarity
   
   %fprintf(1,'Planar structure detected: r=%f\n',r);

   % Transform the plane to bring it in the Z=0 plane:
   
   R_transform = V';
   
   %norm(R_transform(1:2,3))
   
   if norm(R_transform(1:2,3)) < 1e-6,
      R_transform = eye(3);
   end;
   
   if det(R_transform) < 0, R_transform = -R_transform; end;
   
	T_transform = -(R_transform)*X_mean;

	X_new = R_transform*X_kk + T_transform*ones(1,Np);
   
   
   % Compute the planar homography:
   
   H = compute_homographyGPP(xn,X_new(1:2,:));
   
   % De-embed the motion parameters from the homography:
   
   sc = mean([norm(H(:,1));norm(H(:,2))]);
   
   H = H/sc;
   
   % Extra normalization for some reasons...
   %H(:,1) = H(:,1)/norm(H(:,1));
   %H(:,2) = H(:,2)/norm(H(:,2));
   
   if 0, %%% Some tests for myself... the opposite sign solution leads to negative depth!!!
       
       % Case#1: no opposite sign:
       
       omckk1 = rodriguesGPP([H(:,1:2) cross(H(:,1),H(:,2))]);
       Rckk1 = rodriguesGPP(omckk1);
       Tckk1 = H(:,3);
       
       Hs1 = [Rckk1(:,1:2) Tckk1];
       xn1 = Hs1*[X_new(1:2,:);ones(1,Np)];
       xn1 = [xn1(1,:)./xn1(3,:) ; xn1(2,:)./xn1(3,:)];
       e1 = xn1 - xn;
       
       % Case#2: opposite sign:
       
       omckk2 = rodriguesGPP([-H(:,1:2) cross(H(:,1),H(:,2))]);
       Rckk2 = rodriguesGPP(omckk2);
       Tckk2 = -H(:,3);
       
       Hs2 = [Rckk2(:,1:2) Tckk2];
       xn2 = Hs2*[X_new(1:2,:);ones(1,Np)];
       xn2 = [xn2(1,:)./xn2(3,:) ; xn2(2,:)./xn2(3,:)];
       e2 = xn2 - xn;
       
       if 1, %norm(e1) < norm(e2),
           omckk = omckk1;
           Tckk = Tckk1;
           Rckk = Rckk1;
       else
           omckk = omckk2;
           Tckk = Tckk2;
           Rckk = Rckk2;
       end;
       
   else
       
       u1 = H(:,1);
       u1 = u1 / norm(u1);
       u2 = H(:,2) - dot(u1,H(:,2)) * u1;
       u2 = u2 / norm(u2);
       u3 = cross(u1,u2);
       RRR = [u1 u2 u3];
       omckk = rodriguesGPP(RRR);

       %omckk = rodrigues([H(:,1:2) cross(H(:,1),H(:,2))]);
       Rckk = rodriguesGPP(omckk);
       Tckk = H(:,3);
       
   end;
   
      
   
   %If Xc = Rckk * X_new + Tckk, then Xc = Rckk * R_transform * X_kk + Tckk + T_transform
   
   Tckk = Tckk + Rckk* T_transform;
   Rckk = Rckk * R_transform;

   omckk = rodriguesGPP(Rckk);
   Rckk = rodriguesGPP(omckk);
   
   
else
   
   %fprintf(1,'Non planar structure detected: r=%f\n',r);

   % Computes an initial guess for extrinsic parameters (works for general 3d structure, not planar!!!):
   % The DLT method is applied here!!
   
   J = zeros(2*Np,12);
	
	xX = (ones(3,1)*xn(1,:)).*X_kk;
	yX = (ones(3,1)*xn(2,:)).*X_kk;
	
	J(1:2:end,[1 4 7]) = -X_kk';
	J(2:2:end,[2 5 8]) = X_kk';
	J(1:2:end,[3 6 9]) = xX';
	J(2:2:end,[3 6 9]) = -yX';
	J(1:2:end,12) = xn(1,:)';
	J(2:2:end,12) = -xn(2,:)';
	J(1:2:end,10) = -ones(Np,1);
	J(2:2:end,11) = ones(Np,1);
	
	JJ = J'*J;
	[U,S,V] = svd(JJ);
   
   RR = reshape(V(1:9,12),3,3);
   
   if det(RR) < 0,
      V(:,12) = -V(:,12);
      RR = -RR;
   end;
   
   [Ur,Sr,Vr] = svd(RR);
   
   Rckk = Ur*Vr';
   
   sc = norm(V(1:9,12)) / norm(Rckk(:));
   Tckk = V(10:12,12)/sc;
   
	omckk = rodriguesGPP(Rckk);
   Rckk = rodriguesGPP(omckk);
   
end
end

function [omckk,Tckk,Rckk,JJ] = compute_extrinsic_refineGPP(omc_init,Tc_init,x_kk,X_kk,fc,cc,kc,alpha_c,MaxIter,thresh_cond),

%compute_extrinsic
%
%[omckk,Tckk,Rckk] = compute_extrinsic_refine(omc_init,x_kk,X_kk,fc,cc,kc,alpha_c,MaxIter)
%
%Computes the extrinsic parameters attached to a 3D structure X_kk given its projection
%on the image plane x_kk and the intrinsic camera parameters fc, cc and kc.
%Works with planar and non-planar structures.
%
%INPUT: x_kk: Feature locations on the images
%       X_kk: Corresponding grid coordinates
%       fc: Camera focal length
%       cc: Principal point coordinates
%       kc: Distortion coefficients
%       alpha_c: Skew coefficient
%       MaxIter: Maximum number of iterations
%
%OUTPUT: omckk: 3D rotation vector attached to the grid positions in space
%        Tckk: 3D translation vector attached to the grid positions in space
%        Rckk: 3D rotation matrices corresponding to the omc vectors

%
%Method: Computes the normalized point coordinates, then computes the 3D pose
%
%Important functions called within that program:
%
%normalize: Computes the normalize image point coordinates.
%
%pose3D: Computes the 3D pose of the structure given the normalized image projection.
%
%project_points.m: Computes the 2D image projections of a set of 3D points


if nargin < 10,
   thresh_cond = inf;
end;


if nargin < 9,
   MaxIter = 20;
end;

if nargin < 8,
    alpha_c = 0;
    if nargin < 7,
        kc = zeros(5,1);
        if nargin < 6,
            cc = zeros(2,1);
            if nargin < 5,
                fc = ones(2,1);
                if nargin < 4,
                    error('Need 2D projections and 3D points (in compute_extrinsic_refine.m)');
                    return;
                end;
            end;
        end;
    end;
end;


% Initialization:

omckk = omc_init;
Tckk = Tc_init;


% Final optimization (minimize the reprojection error in pixel):
% through Gradient Descent:

param = [omckk;Tckk];

change = 1;

iter = 0;

%keyboard;

%fprintf(1,'Gradient descent iterations: ');

while (change > 1e-10)&(iter < MaxIter),
    
    %fprintf(1,'%d...',iter+1);
    
    [x,dxdom,dxdT] = project_points2GPP(X_kk,omckk,Tckk,fc,cc,kc,alpha_c);
    
    ex = x_kk - x;
    
    %keyboard;
    
    JJ = [dxdom dxdT];
    
    if cond(JJ) > thresh_cond,
        change = 0;
    else
        
        JJ2 = JJ'*JJ;
        
        param_innov = inv(JJ2)*(JJ')*ex(:);
        param_up = param + param_innov;
        change = norm(param_innov)/norm(param_up);
        param = param_up;
        iter = iter + 1;
        
        omckk = param(1:3);
        Tckk = param(4:6);
        
    end;
    
end;

%fprintf(1,'\n');

Rckk = rodriguesGPP(omckk);

end




function [xp,dxpdom,dxpdT,dxpdf,dxpdc,dxpdk,dxpdalpha] = project_points2GPP(X,om,T,f,c,k,alpha)

%project_points2.m
%
%[xp,dxpdom,dxpdT,dxpdf,dxpdc,dxpdk] = project_points2(X,om,T,f,c,k,alpha)
%
%Projects a 3D structure onto the image plane.
%
%INPUT: X: 3D structure in the world coordinate frame (3xN matrix for N points)
%       (om,T): Rigid motion parameters between world coordinate frame and camera reference frame
%               om: rotation vector (3x1 vector); T: translation vector (3x1 vector)
%       f: camera focal length in units of horizontal and vertical pixel units (2x1 vector)
%       c: principal point location in pixel units (2x1 vector)
%       k: Distortion coefficients (radial and tangential) (4x1 vector)
%       alpha: Skew coefficient between x and y pixel (alpha = 0 <=> square pixels)
%
%OUTPUT: xp: Projected pixel coordinates (2xN matrix for N points)
%        dxpdom: Derivative of xp with respect to om ((2N)x3 matrix)
%        dxpdT: Derivative of xp with respect to T ((2N)x3 matrix)
%        dxpdf: Derivative of xp with respect to f ((2N)x2 matrix if f is 2x1, or (2N)x1 matrix is f is a scalar)
%        dxpdc: Derivative of xp with respect to c ((2N)x2 matrix)
%        dxpdk: Derivative of xp with respect to k ((2N)x4 matrix)
%
%Definitions:
%Let P be a point in 3D of coordinates X in the world reference frame (stored in the matrix X)
%The coordinate vector of P in the camera reference frame is: Xc = R*X + T
%where R is the rotation matrix corresponding to the rotation vector om: R = rodrigues(om);
%call x, y and z the 3 coordinates of Xc: x = Xc(1); y = Xc(2); z = Xc(3);
%The pinehole projection coordinates of P is [a;b] where a=x/z and b=y/z.
%call r^2 = a^2 + b^2.
%The distorted point coordinates are: xd = [xx;yy] where:
%
%xx = a * (1 + kc(1)*r^2 + kc(2)*r^4 + kc(5)*r^6)      +      2*kc(3)*a*b + kc(4)*(r^2 + 2*a^2);
%yy = b * (1 + kc(1)*r^2 + kc(2)*r^4 + kc(5)*r^6)      +      kc(3)*(r^2 + 2*b^2) + 2*kc(4)*a*b;
%
%The left terms correspond to radial distortion (6th degree), the right terms correspond to tangential distortion
%
%Finally, convertion into pixel coordinates: The final pixel coordinates vector xp=[xxp;yyp] where:
%
%xxp = f(1)*(xx + alpha*yy) + c(1)
%yyp = f(2)*yy + c(2)
%
%
%NOTE: About 90 percent of the code takes care fo computing the Jacobian matrices
%
%
%Important function called within that program:
%
%rodrigues.m: Computes the rotation matrix corresponding to a rotation vector
%
%rigid_motion.m: Computes the rigid motion transformation of a given structure


if nargin < 7,
   alpha = 0;
   if nargin < 6,
      k = zeros(5,1);
      if nargin < 5,
         c = zeros(2,1);
         if nargin < 4,
            f = ones(2,1);
            if nargin < 3,
               T = zeros(3,1);
               if nargin < 2,
                  om = zeros(3,1);
                  if nargin < 1,
                     error('Need at least a 3D structure to project (in project_points.m)');
                     return;
                  end;
               end;
            end;
         end;
      end;
   end;
end;


[m,n] = size(X);

[Y,dYdom,dYdT] = rigid_motionGPP(X,om,T);


inv_Z = 1./Y(3,:);

x = (Y(1:2,:) .* (ones(2,1) * inv_Z)) ;


bb = (-x(1,:) .* inv_Z)'*ones(1,3);
cc = (-x(2,:) .* inv_Z)'*ones(1,3);


dxdom = zeros(2*n,3);
dxdom(1:2:end,:) = ((inv_Z')*ones(1,3)) .* dYdom(1:3:end,:) + bb .* dYdom(3:3:end,:);
dxdom(2:2:end,:) = ((inv_Z')*ones(1,3)) .* dYdom(2:3:end,:) + cc .* dYdom(3:3:end,:);

dxdT = zeros(2*n,3);
dxdT(1:2:end,:) = ((inv_Z')*ones(1,3)) .* dYdT(1:3:end,:) + bb .* dYdT(3:3:end,:);
dxdT(2:2:end,:) = ((inv_Z')*ones(1,3)) .* dYdT(2:3:end,:) + cc .* dYdT(3:3:end,:);


% Add distortion:

r2 = x(1,:).^2 + x(2,:).^2;

dr2dom = 2*((x(1,:)')*ones(1,3)) .* dxdom(1:2:end,:) + 2*((x(2,:)')*ones(1,3)) .* dxdom(2:2:end,:);
dr2dT = 2*((x(1,:)')*ones(1,3)) .* dxdT(1:2:end,:) + 2*((x(2,:)')*ones(1,3)) .* dxdT(2:2:end,:);


r4 = r2.^2;

dr4dom = 2*((r2')*ones(1,3)) .* dr2dom;
dr4dT = 2*((r2')*ones(1,3)) .* dr2dT;


r6 = r2.^3;

dr6dom = 3*((r2'.^2)*ones(1,3)) .* dr2dom;
dr6dT = 3*((r2'.^2)*ones(1,3)) .* dr2dT;


% Radial distortion:

cdist = 1 + k(1) * r2 + k(2) * r4 + k(5) * r6;

dcdistdom = k(1) * dr2dom + k(2) * dr4dom + k(5) * dr6dom;
dcdistdT = k(1) * dr2dT + k(2) * dr4dT + k(5) * dr6dT;
dcdistdk = [ r2' r4' zeros(n,2) r6'];


xd1 = x .* (ones(2,1)*cdist);

dxd1dom = zeros(2*n,3);
dxd1dom(1:2:end,:) = (x(1,:)'*ones(1,3)) .* dcdistdom;
dxd1dom(2:2:end,:) = (x(2,:)'*ones(1,3)) .* dcdistdom;
coeff = (reshape([cdist;cdist],2*n,1)*ones(1,3));
dxd1dom = dxd1dom + coeff.* dxdom;

dxd1dT = zeros(2*n,3);
dxd1dT(1:2:end,:) = (x(1,:)'*ones(1,3)) .* dcdistdT;
dxd1dT(2:2:end,:) = (x(2,:)'*ones(1,3)) .* dcdistdT;
dxd1dT = dxd1dT + coeff.* dxdT;

dxd1dk = zeros(2*n,5);
dxd1dk(1:2:end,:) = (x(1,:)'*ones(1,5)) .* dcdistdk;
dxd1dk(2:2:end,:) = (x(2,:)'*ones(1,5)) .* dcdistdk;



% tangential distortion:

a1 = 2.*x(1,:).*x(2,:);
a2 = r2 + 2*x(1,:).^2;
a3 = r2 + 2*x(2,:).^2;

delta_x = [k(3)*a1 + k(4)*a2 ;
   k(3) * a3 + k(4)*a1];


%ddelta_xdx = zeros(2*n,2*n);
aa = (2*k(3)*x(2,:)+6*k(4)*x(1,:))'*ones(1,3);
bb = (2*k(3)*x(1,:)+2*k(4)*x(2,:))'*ones(1,3);
cc = (6*k(3)*x(2,:)+2*k(4)*x(1,:))'*ones(1,3);

ddelta_xdom = zeros(2*n,3);
ddelta_xdom(1:2:end,:) = aa .* dxdom(1:2:end,:) + bb .* dxdom(2:2:end,:);
ddelta_xdom(2:2:end,:) = bb .* dxdom(1:2:end,:) + cc .* dxdom(2:2:end,:);

ddelta_xdT = zeros(2*n,3);
ddelta_xdT(1:2:end,:) = aa .* dxdT(1:2:end,:) + bb .* dxdT(2:2:end,:);
ddelta_xdT(2:2:end,:) = bb .* dxdT(1:2:end,:) + cc .* dxdT(2:2:end,:);

ddelta_xdk = zeros(2*n,5);
ddelta_xdk(1:2:end,3) = a1';
ddelta_xdk(1:2:end,4) = a2';
ddelta_xdk(2:2:end,3) = a3';
ddelta_xdk(2:2:end,4) = a1';



xd2 = xd1 + delta_x;

dxd2dom = dxd1dom + ddelta_xdom ;
dxd2dT = dxd1dT + ddelta_xdT;
dxd2dk = dxd1dk + ddelta_xdk ;


% Add Skew:

xd3 = [xd2(1,:) + alpha*xd2(2,:);xd2(2,:)];

% Compute: dxd3dom, dxd3dT, dxd3dk, dxd3dalpha

dxd3dom = zeros(2*n,3);
dxd3dom(1:2:2*n,:) = dxd2dom(1:2:2*n,:) + alpha*dxd2dom(2:2:2*n,:);
dxd3dom(2:2:2*n,:) = dxd2dom(2:2:2*n,:);
dxd3dT = zeros(2*n,3);
dxd3dT(1:2:2*n,:) = dxd2dT(1:2:2*n,:) + alpha*dxd2dT(2:2:2*n,:);
dxd3dT(2:2:2*n,:) = dxd2dT(2:2:2*n,:);
dxd3dk = zeros(2*n,5);
dxd3dk(1:2:2*n,:) = dxd2dk(1:2:2*n,:) + alpha*dxd2dk(2:2:2*n,:);
dxd3dk(2:2:2*n,:) = dxd2dk(2:2:2*n,:);
dxd3dalpha = zeros(2*n,1);
dxd3dalpha(1:2:2*n,:) = xd2(2,:)';



% Pixel coordinates:
if length(f)>1,
    xp = xd3 .* (f * ones(1,n))  +  c*ones(1,n);
    coeff = reshape(f*ones(1,n),2*n,1);
    dxpdom = (coeff*ones(1,3)) .* dxd3dom;
    dxpdT = (coeff*ones(1,3)) .* dxd3dT;
    dxpdk = (coeff*ones(1,5)) .* dxd3dk;
    dxpdalpha = (coeff) .* dxd3dalpha;
    dxpdf = zeros(2*n,2);
    dxpdf(1:2:end,1) = xd3(1,:)';
    dxpdf(2:2:end,2) = xd3(2,:)';
else
    xp = f * xd3 + c*ones(1,n);
    dxpdom = f  * dxd3dom;
    dxpdT = f * dxd3dT;
    dxpdk = f  * dxd3dk;
    dxpdalpha = f .* dxd3dalpha;
    dxpdf = xd3(:);
end;

dxpdc = zeros(2*n,2);
dxpdc(1:2:end,1) = ones(n,1);
dxpdc(2:2:end,2) = ones(n,1);


return;

% Test of the Jacobians:

n = 10;

X = 10*randn(3,n);
om = randn(3,1);
T = [10*randn(2,1);40];
f = 1000*rand(2,1);
c = 1000*randn(2,1);
k = 0.5*randn(5,1);
alpha = 0.01*randn(1,1);

[x,dxdom,dxdT,dxdf,dxdc,dxdk,dxdalpha] = project_points2GPP(X,om,T,f,c,k,alpha);


% Test on om: OK

dom = 0.000000001 * norm(om)*randn(3,1);
om2 = om + dom;

[x2] = project_points2GPP(X,om2,T,f,c,k,alpha);

x_pred = x + reshape(dxdom * dom,2,n);


norm(x2-x)/norm(x2 - x_pred)


% Test on T: OK!!

dT = 0.0001 * norm(T)*randn(3,1);
T2 = T + dT;

[x2] = project_points2GPP(X,om,T2,f,c,k,alpha);

x_pred = x + reshape(dxdT * dT,2,n);


norm(x2-x)/norm(x2 - x_pred)



% Test on f: OK!!

df = 0.001 * norm(f)*randn(2,1);
f2 = f + df;

[x2] = project_points2GPP(X,om,T,f2,c,k,alpha);

x_pred = x + reshape(dxdf * df,2,n);


norm(x2-x)/norm(x2 - x_pred)


% Test on c: OK!!

dc = 0.01 * norm(c)*randn(2,1);
c2 = c + dc;

[x2] = project_points2GPP(X,om,T,f,c2,k,alpha);

x_pred = x + reshape(dxdc * dc,2,n);

norm(x2-x)/norm(x2 - x_pred)

% Test on k: OK!!

dk = 0.001 * norm(k)*randn(5,1);
k2 = k + dk;

[x2] = project_points2GPP(X,om,T,f,c,k2,alpha);

x_pred = x + reshape(dxdk * dk,2,n);

norm(x2-x)/norm(x2 - x_pred)


% Test on alpha: OK!!

dalpha = 0.001 * norm(k)*randn(1,1);
alpha2 = alpha + dalpha;

[x2] = project_points2GPP(X,om,T,f,c,k,alpha2);

x_pred = x + reshape(dxdalpha * dalpha,2,n);

norm(x2-x)/norm(x2 - x_pred)

end




function	[out,dout]=rodriguesGPP(in)

% RODRIGUES	Transform rotation matrix into rotation vector and viceversa.
%		
%		Sintax:  [OUT]=RODRIGUES(IN)
% 		If IN is a 3x3 rotation matrix then OUT is the
%		corresponding 3x1 rotation vector
% 		if IN is a rotation 3-vector then OUT is the 
%		corresponding 3x3 rotation matrix
%

%%
%%		Copyright (c) March 1993 -- Pietro Perona
%%		California Institute of Technology
%%

%% ALL CHECKED BY JEAN-YVES BOUGUET, October 1995.
%% FOR ALL JACOBIAN MATRICES !!! LOOK AT THE TEST AT THE END !!

%% BUG when norm(om)=pi fixed -- April 6th, 1997;
%% Jean-Yves Bouguet

%% Add projection of the 3x3 matrix onto the set of special ortogonal matrices SO(3) by SVD -- February 7th, 2003;
%% Jean-Yves Bouguet

[m,n] = size(in);
%bigeps = 10e+4*eps;
bigeps = 10e+20*eps;

if ((m==1) & (n==3)) | ((m==3) & (n==1)) %% it is a rotation vector
   theta = norm(in);
   if theta < eps
      R = eye(3);
      
      %if nargout > 1,
      
      dRdin = [0 0 0;
	       0 0 1;
	       0 -1 0;
	       0 0 -1;
	       0 0 0;
	       1 0 0;
	       0 1 0;
	       -1 0 0;
          0 0 0];
       
       %end;
	 
   else
      if n==length(in)  in=in'; end; 	%% make it a column vec. if necess.
	 
	 %m3 = [in,theta]

	 dm3din = [eye(3);in'/theta];

	 omega = in/theta;
	 
	 %m2 = [omega;theta]
	 
	 dm2dm3 = [eye(3)/theta -in/theta^2; zeros(1,3) 1];
	 
	 alpha = cos(theta);
	 beta = sin(theta);
	 gamma = 1-cos(theta);
	 omegav=[[0 -omega(3) omega(2)];[omega(3) 0 -omega(1)];[-omega(2) omega(1) 0 ]];
	 A = omega*omega';
	 
	 %m1 = [alpha;beta;gamma;omegav;A];
	 
	 dm1dm2 = zeros(21,4);
	 dm1dm2(1,4) = -sin(theta);
	 dm1dm2(2,4) = cos(theta);
	 dm1dm2(3,4) = sin(theta);
	 dm1dm2(4:12,1:3) = [0 0 0 0 0 1 0 -1 0;
	                     0 0 -1 0 0 0 1 0 0;
			     0 1 0 -1 0 0 0 0 0]';
		       
         w1 = omega(1);
	 w2 = omega(2);
	 w3 = omega(3);
	 
	 dm1dm2(13:21,1) = [2*w1;w2;w3;w2;0;0;w3;0;0];
	 dm1dm2(13: 21,2) = [0;w1;0;w1;2*w2;w3;0;w3;0];
	 dm1dm2(13:21,3) = [0;0;w1;0;0;w2;w1;w2;2*w3];
	 
	 R = eye(3)*alpha + omegav*beta + A*gamma;
	 
	 dRdm1 = zeros(9,21);
	 
	 dRdm1([1 5 9],1) = ones(3,1);
	 dRdm1(:,2) = omegav(:);
	 dRdm1(:,4:12) = beta*eye(9);
	 dRdm1(:,3) = A(:);
	 dRdm1(:,13:21) = gamma*eye(9);
	 
	 dRdin = dRdm1 * dm1dm2 * dm2dm3 * dm3din;
	 
	 
      end;
      out = R;
      dout = dRdin;
      
      %% it is prob. a rot matr.
   elseif ((m==n) & (m==3) & (norm(in' * in - eye(3)) < bigeps)...
	    & (abs(det(in)-1) < bigeps))
      R = in;
      
      % project the rotation matrix to SO(3);
      [U,S,V] = svd(R);
      R = U*V';
      
      tr = (trace(R)-1)/2;
      dtrdR = [1 0 0 0 1 0 0 0 1]/2;
      theta = real(acos(tr));
      
      
      if sin(theta) >= 1e-5,
	 
	 dthetadtr = -1/sqrt(1-tr^2);
	 
	 dthetadR = dthetadtr * dtrdR;
	 % var1 = [vth;theta];
	 vth = 1/(2*sin(theta));
	 dvthdtheta = -vth*cos(theta)/sin(theta);
	 dvar1dtheta = [dvthdtheta;1];
	 
	 dvar1dR =  dvar1dtheta * dthetadR;
	 
	 
	 om1 = [R(3,2)-R(2,3), R(1,3)-R(3,1), R(2,1)-R(1,2)]';
	 
	 dom1dR = [0 0 0 0 0 1 0 -1 0;
	       0 0 -1 0 0 0 1 0 0;
	       0 1 0 -1 0 0 0 0 0];
	 
	 % var = [om1;vth;theta];
	 dvardR = [dom1dR;dvar1dR];
	 
	 % var2 = [om;theta];
	 om = vth*om1;
	 domdvar = [vth*eye(3) om1 zeros(3,1)];
	 dthetadvar = [0 0 0 0 1];
	 dvar2dvar = [domdvar;dthetadvar];
	 
	 
	 out = om*theta;
	 domegadvar2 = [theta*eye(3) om];
	 
	 dout = domegadvar2 * dvar2dvar * dvardR;
	 
	 
      else
	 if tr > 0; 			% case norm(om)=0;
	    
	    out = [0 0 0]';
	    
	    dout = [0 0 0 0 0 1/2 0 -1/2 0;
		  0 0 -1/2 0 0 0 1/2 0 0;
		  0 1/2 0 -1/2 0 0 0 0 0];
	 else 				% case norm(om)=pi; %% fixed April 6th
	    
	    
	    out = theta * (sqrt((diag(R)+1)/2).*[1;2*(R(1,2:3)>=0)'-1]);
	    %keyboard;
	    
	    if nargout > 1,
	       fprintf(1,'WARNING!!!! Jacobian domdR undefined!!!\n');
		 	dout = NaN*ones(3,9);
	    end;
	 end; 
      end;
      
   else
      error('Neither a rotation matrix nor a rotation vector were provided');
   end;

return;

%% test of the Jacobians:

%%%% TEST OF dRdom:
om = randn(3,1);
dom = randn(3,1)/1000000;

[R1,dR1] = rodriguesGPP(om);
R2 = rodriguesGPP(om+dom);

R2a = R1 + reshape(dR1 * dom,3,3);

gain = norm(R2 - R1)/norm(R2 - R2a)

%%% TEST OF dOmdR:
om = randn(3,1);
R = rodriguesGPP(om);
dom = randn(3,1)/10000;
dR = rodriguesGPP(om+dom) - R;

[omc,domdR] = rodriguesGPP(R);
[om2] = rodriguesGPP(R+dR);

om_app = omc + domdR*dR(:);

gain = norm(om2 - omc)/norm(om2 - om_app)


%%% OTHER BUG: (FIXED NOW!!!)

omu = randn(3,1);   
omu = omu/norm(omu)
om = pi*omu;        
[R,dR]= rodriguesGPP(om);
[om2] = rodriguesGPP(R);
[om om2]

%%% NORMAL OPERATION

om = randn(3,1);         
[R,dR]= rodriguesGPP(om);
[om2] = rodriguesGPP(R);
[om om2]

end




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

end



function [fc,cc,kc,alpha_c]=convert_heikkila_stoyanovGPP(par)

% Convert from a Heikkila-style 'par' file to the parameters needed in
% the stoyanov findpos and rectification files. 

            su=par(1);      %scale factor ~1 (Asp)
            f=par(2);       %effective focal length (f)
            u0=par(3);      %principal point (Cpx)
            v0=par(4);      %principal point (Cpy)
            k1=par(5);      %radial distortion coefficients
            k2=par(6);      %radial distortion coefficients
            p1=par(7);      %tangential distortion coefficients
            p2=par(8);      %tangential distortion coefficients
%Entering internal lens parameters calibrated for lens#1
% 	sys=configc(cam);
%   sys = [
%       2016,     %number of pixels in horizontal direction 753
%       720,     %number of pixels in vertical direction 582
%       6.4512,     %effective CCD chip size in horizontal direction 
%       4.8896,     %effective CCD chip size in vertical direction
%       6.0,       %nominal focal length
%       0,       %radius of the circular control points
%       0,       %for future expansions
%       0,
%       0,
%       ];
  sys = [
      1600,     %number of pixels in horizontal direction 753
      728,     %number of pixels in vertical direction 582
      4.536,     %effective CCD chip size in horizontal direction 
      3.416,     %effective CCD chip size in vertical direction
      2.5,       %nominal focal length
      0,       %radius of the circular control points
      0,       %for future expansions
      0,
      0,
      ];
	NDX=sys(1); NDY=sys(2); Sx=sys(3); Sy=sys(4);
 	    
%Scaling the results to the final image coordinates.
    Du=NDX/Sx;
    Dv=NDY/Sy;
    
fc(1,1)=f.*Du.*su;
fc(2,1)=f.*Dv;
cc(1,1)=u0;
cc(2,1)=v0;
alpha_c=0;
kc(1,1)=f.^3.*k1;
kc(2,1)=f.^5.*k2;
kc(3,1)=f.^2.*p1;
kc(4,1)=f.^2.*p2;
kc(5,1)=0;

end


function [Y,dYdom,dYdT] = rigid_motionGPP(X,om,T);

%rigid_motion.m
%
%[Y,dYdom,dYdT] = rigid_motion(X,om,T)
%
%Computes the rigid motion transformation Y = R*X+T, where R = rodrigues(om).
%
%INPUT: X: 3D structure in the world coordinate frame (3xN matrix for N points)
%       (om,T): Rigid motion parameters between world coordinate frame and camera reference frame
%               om: rotation vector (3x1 vector); T: translation vector (3x1 vector)
%
%OUTPUT: Y: 3D coordinates of the structure points in the camera reference frame (3xN matrix for N points)
%        dYdom: Derivative of Y with respect to om ((3N)x3 matrix)
%        dYdT: Derivative of Y with respect to T ((3N)x3 matrix)
%
%Definitions:
%Let P be a point in 3D of coordinates X in the world reference frame (stored in the matrix X)
%The coordinate vector of P in the camera reference frame is: Y = R*X + T
%where R is the rotation matrix corresponding to the rotation vector om: R = rodrigues(om);
%
%Important function called within that program:
%
%rodrigues.m: Computes the rotation matrix corresponding to a rotation vector



if nargin < 3,
   T = zeros(3,1);
   if nargin < 2,
      om = zeros(3,1);
      if nargin < 1,
         error('Need at least a 3D structure as input (in rigid_motion.m)');
         return;
      end;
   end;
end;


[R,dRdom] = rodriguesGPP(om);

[m,n] = size(X);

Y = R*X + repmat(T,[1 n]);

if nargout > 1,
   

dYdR = zeros(3*n,9);
dYdT = zeros(3*n,3);

dYdR(1:3:end,1:3:end) =  X';
dYdR(2:3:end,2:3:end) =  X';
dYdR(3:3:end,3:3:end) =  X';

dYdT(1:3:end,1) =  ones(n,1);
dYdT(2:3:end,2) =  ones(n,1);
dYdT(3:3:end,3) =  ones(n,1);

dYdom = dYdR * dRdom;

end
end


function [H,Hnorm,inv_Hnorm] = compute_homographyGPP(m,M);

%compute_homography
%
%[H,Hnorm,inv_Hnorm] = compute_homography(m,M)
%
%Computes the planar homography between the point coordinates on the plane (M) and the image
%point coordinates (m).
%
%INPUT: m: homogeneous coordinates in the image plane (3xN matrix)
%       M: homogeneous coordinates in the plane in 3D (3xN matrix)
%
%OUTPUT: H: Homography matrix (3x3 homogeneous matrix)
%        Hnorm: Normalization matrix used on the points before homography computation
%               (useful for numerical stability is points in pixel coordinates)
%        inv_Hnorm: The inverse of Hnorm
%
%Definition: m ~ H*M where "~" means equal up to a non zero scalar factor.
%
%Method: First computes an initial guess for the homography through quasi-linear method.
%        Then, if the total number of points is larger than 4, optimize the solution by minimizing
%        the reprojection error (in the least squares sense).
%
%
%Important functions called within that program:
%
%comp_distortion_oulu: Undistorts pixel coordinates.
%
%compute_homography.m: Computes the planar homography between points on the grid in 3D, and the image plane.
%
%project_points.m: Computes the 2D image projections of a set of 3D points, and also returns te Jacobian
%                  matrix (derivative with respect to the intrinsic and extrinsic parameters).
%                  This function is called within the minimization loop.




Np = size(m,2);

if size(m,1)<3,
   m = [m;ones(1,Np)];
end;

if size(M,1)<3,
   M = [M;ones(1,Np)];
end;


m = m ./ (ones(3,1)*m(3,:));
M = M ./ (ones(3,1)*M(3,:));

% Prenormalization of point coordinates (very important):
% (Affine normalization)

ax = m(1,:);
ay = m(2,:);

mxx = mean(ax);
myy = mean(ay);
ax = ax - mxx;
ay = ay - myy;

scxx = mean(abs(ax));
scyy = mean(abs(ay));


Hnorm = [1/scxx 0 -mxx/scxx;0 1/scyy -myy/scyy;0 0 1];
inv_Hnorm = [scxx 0 mxx ; 0 scyy myy; 0 0 1];

mn = Hnorm*m;

% Compute the homography between m and mn:

% Build the matrix:

L = zeros(2*Np,9);

L(1:2:2*Np,1:3) = M';
L(2:2:2*Np,4:6) = M';
L(1:2:2*Np,7:9) = -((ones(3,1)*mn(1,:)).* M)';
L(2:2:2*Np,7:9) = -((ones(3,1)*mn(2,:)).* M)';

if Np > 4,
	L = L'*L;
end;

[U,S,V] = svd(L);

hh = V(:,9);
hh = hh/hh(9);

Hrem = reshape(hh,3,3)';
%Hrem = Hrem / Hrem(3,3);


% Final homography:

H = inv_Hnorm*Hrem;

if 0,
   m2 = H*M;
   m2 = [m2(1,:)./m2(3,:) ; m2(2,:)./m2(3,:)];
   merr = m(1:2,:) - m2;
end;

%keyboard;
 
%%% Homography refinement if there are more than 4 points:

if Np > 4,
   
   % Final refinement:
   hhv = reshape(H',9,1);
   hhv = hhv(1:8);
   
   for iter=1:10,
      

   
		mrep = H * M;

		J = zeros(2*Np,8);

		MMM = (M ./ (ones(3,1)*mrep(3,:)));

		J(1:2:2*Np,1:3) = -MMM';
		J(2:2:2*Np,4:6) = -MMM';
		
		mrep = mrep ./ (ones(3,1)*mrep(3,:));

		m_err = m(1:2,:) - mrep(1:2,:);
		m_err = m_err(:);

		MMM2 = (ones(3,1)*mrep(1,:)) .* MMM;
		MMM3 = (ones(3,1)*mrep(2,:)) .* MMM;

		J(1:2:2*Np,7:8) = MMM2(1:2,:)';
		J(2:2:2*Np,7:8) = MMM3(1:2,:)';

		MMM = (M ./ (ones(3,1)*mrep(3,:)))';

		hh_innov  = inv(J'*J)*J'*m_err;

		hhv_up = hhv - hh_innov;

		H_up = reshape([hhv_up;1],3,3)';

		%norm(m_err)
		%norm(hh_innov)

		hhv = hhv_up;
      H = H_up;
      
   end;
   

end;

if 0,
   m2 = H*M;
   m2 = [m2(1,:)./m2(3,:) ; m2(2,:)./m2(3,:)];
   merr = m(1:2,:) - m2;
end;

return;

%test of Jacobian

mrep = H*M;
mrep = mrep ./ (ones(3,1)*mrep(3,:));

m_err = mrep(1:2,:) - m(1:2,:);
figure(8);
plot(m_err(1,:),m_err(2,:),'r+');
std(m_err')
end

function [x] = comp_distortion_ouluGPP(xd,k);

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
            
    end;
    
end;
end


function [x_comp]  = comp_distortionGPP(x_dist,k2);

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
end