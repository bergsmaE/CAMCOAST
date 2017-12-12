% GPP_direction

function CAMS_N2_Direction(dirN1,dirN2_d,ls_maj_d)
% On liste toutes les images 'A' à traiter 
files=[];
datef=[];
if size(ls_maj_d,1) > 0;
    
for i = 1 : size(ls_maj_d,1);
    lsA = ls([dirN1,ls_maj_d(i,:),'/A*']);
    for j = 1 : size(lsA,1)
        files = [files [dirN1 ls_maj_d(i,:) '\' lsA(j,:)]'];
    end
end

%Datef: liste des dates en jour julien des images à traiter
datef = datenum(files(end-15:end-4,:)','yyyymmddHHMM')';

%Ymin,Ymax, Xmin et Xmax dépendent du Rectfile!
Ymin = -150.0;
Ymax =  -80.0;
Xmin =  0.0;
Xmax =  120;
method = 'linear';     % Method for griddata
Normalisation = 0;     % Normalisation or not of the images
GridRes = 0.5;         % res is defined in zone

% Position du point 0 en lambert 3 et angle de la plage
% Coordonnées video
ang=90+8;
ang=ang*pi/180;
x0=370341;%Easting
y0=694135;%Northing

clear rr gg bb x1 y1
load('RectGPP.mat');

% Application de la zone et passage du tableau en vecteur (1 colonne) pour les
% 3 composantes rouge, vert, bleu
load('ZoneGPP.mat');

x1=x1(1:1:end);
y1=y1(1:1:end);

% [X_I] = rectify0707(double([y1;x1]);,Rckk,Tckk,fc,cc,kc,alpha_c,0*ones(1,length(x1)));
% function [X_r] = rectify0707(aaa,Rckk,Tckk,fc,cc,kc,alpha_c,Z)
% aaa contains the pixel coordinates (origin top left, horizontal in the
% first row vertical in the second row).
% Z is the level you want to recitfy to (can be
% different for each point, in which case length(Z) must be the same as X_kk))

aaa=double([y1;x1]);
Z=0*ones(1,length(x1));

[n,m]=size(aaa);
if m==2&n~=2
    x_kk=aaa';
elseif m~=2&n==2
    x_kk=aaa;
end
[n,m]=size(aaa);
if length(Z)==1
    Z=Z.*ones(1,m);
end


%%%%%%%%%%%% START NORMALIZE

%undistort the image coordinates and scale
% xn = normalize(aaa,fc,cc,kc,alpha_c);

% function [xn] = normalize(x_kk,fc,cc,kc,alpha_c)

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
% if nargin < 5,
%     alpha_c = 0;
%     if nargin < 4;
%         kc = [0;0;0;0;0];
%         if nargin < 3;
%             cc = [0;0];
%             if nargin < 2,
%                 fc = [1;1];
%             end;
%         end;
%     end;
% end;
% First: Subtract principal point, and divide by the focal length:
x_distort = [(x_kk(1,:) - cc(1))/fc(1);(x_kk(2,:) - cc(2))/fc(2)];

% Second: undo skew
x_distort(1,:) = x_distort(1,:) - alpha_c * x_distort(2,:);

if norm(kc) ~= 0,
    % Third: Compensate for lens distortion:
    
    xd=x_distort;
    k=kc;
    
    %     xn = comp_distortion_oulu(x_distort,kc);
    
    %     function [x] = comp_distortion_oulu(xd,k);
    
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
        
        %         [x] = comp_distortion(xd,k);
        x_dist=xd;
        k2=k;
        
        %         function [x_comp]  = comp_distortion(x_dist,k2);
        
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
            
            [x_comp]  = comp_distortion_oulu(x_dist,k2);
            
        else
            
            radius_2= x_dist(1,:).^2 + x_dist(2,:).^2;
            radial_distortion = 1 + ones(2,1)*(k2 * radius_2);
            radius_2_comp = (x_dist(1,:).^2 + x_dist(2,:).^2) ./ radial_distortion(1,:);
            radial_distortion = 1 + ones(2,1)*(k2 * radius_2_comp);
            x_comp = x_dist ./ radial_distortion;
            
        end
        x=x_comp;
        
    else
        
        k1 = k(1);
        k2 = k(2);
        k3 = k(5);
        p1 = k(3);
        p2 = k(4);
        
        x = xd;                 % initial guess
        
        for kk=1:20,
            
            r_2 = sum(x.^2);
            k_radial =  1 + k1 * r_2 + k2 * r_2.^2 + k3 * r_2.^3;
            delta_x = [2*p1*x(1,:).*x(2,:) + p2*(r_2 + 2*x(1,:).^2);
                p1 * (r_2 + 2*x(2,:).^2)+2*p2*x(1,:).*x(2,:)];
            x = (xd - delta_x)./(ones(2,1)*k_radial);
            
        end;
        
    end;
    xn=x;    
    
else
    xn = x_distort;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% END NORMALIZE
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
X_I=X_r;
% Rotation et translation dans le nouveau repère (x0;y0) (pied de
% l'échafaudage = référence)
yr1=(X_I(1,:)-x0)*cos(ang)+(X_I(2,:)-y0)*sin(ang);
xr1=(X_I(2,:)-y0)*cos(ang)-(X_I(1,:)-x0)*sin(ang);
xx=xr1;
yy=yr1;
% Résolution finale de la figure et nouvelle grille
xr=Xmin:GridRes:Xmax;
yr=Ymin:GridRes:Ymax;
[X,Y] = meshgrid(xr,yr);



angi=[];
dati=[];
count=0;bluetot=0;
files1=files;
%Boucle sur les images afin de calculer les directions des vagues. 
for k=1:length(datef)
try
        fA=strfind(files(:,k)','A_');
        files1(fA,k)='I';
        A=imread(files(:,k)');
        I=imread(files1(:,k)');
        A=I-A;
         
        nx=size(A,1);
        ny=size(A,2);
        rr=double(A(x1(:)+(y1(:)-1)*nx)');
        gg=double(A(x1(:)+(y1(:)-1)*nx+nx*ny)');
        bb=double(A(x1(:)+(y1(:)-1)*nx+2*nx*ny)');
        
        % Normalisation
        if Normalisation == 1
            [rr]=StretchBand(rr,1);
            [gg]=StretchBand(gg,1);
            [bb]=StretchBand(bb,1);
            rr=double(rr);
            gg=double(gg);
            bb=double(bb);
        end
        
        red   = rr;
        green = gg;
        blue  = bb;
        meanr=mean(red);
        meang=mean(green);
        meanb=mean(blue);
        FirstImage=0;
        
        % Merging des images
        clear Arect      
        %   Arect(:,:,1) = griddata(xx,yy,red,X,Y,method);
        %   Arect(:,:,2) = griddata(xx,yy,green,X,Y,method);
        bluetot=bluetot+blue;
        count=count+1;
        if count==4 
            % Ici, on compte le nombre d'images à traiter pour
            % avoir un angle. Dans ce cas (count==4), on veut un angle pour
            % quatre (4) images (ce qui correspond à une heure).
            Arect= griddata(xx(1:10:end),yy(1:10:end),bluetot(1:10:end),X,Y,method);
            figure(5);imagesc(Arect)
            count=0;bluetot=0;
            R=radon(detrend(Arect(:,:)));[rd R]=max(nanstd(R));
            disp(['Angle= ', num2str(90-R), '°'])
            angi=[angi R];
            
            %En cas de changement de jour de calcul, on enregistre un
            %fichier journalier et on réinitialise dati et angi
            if size(dati) > 0;
                if str2num(datestr(datef(k),'yyyymmdd'))~=str2num(datestr(dati(end),'yyyymmdd')) 
                    %Sauvegarde d'un fichier 
                    cd(dirN2_d);
                    save(['GPP_Direction_',datestr(dati(end),'yyyymmdd')],'dati','angi');
                    dati=[];
                    angi=[];
                end  
            end
                dati=[dati datef(k)];
                disp(datestr(datef(k))) 
                

end
end
end
if size(dati) > 0;
                cd(dirN2_d);
                save(['GPP_Direction_',datestr(dati(end),'yyyymmdd')],'dati','angi');
                dati=[];
                angi=[];
end
end
end



