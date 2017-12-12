% get_wave_parameters_20160513.m
%get wave parameters from cross-shore video timestack

function [C,depth,hs,hm,Tm,Tp,Er,rolL,nbwave,X1,Y1,dx,Breakstd,Breakmean1,Breakmean2]=GPP_get_wave_parameters(Img,Rectfile,dt,sk1,sk2,CoordCam,BeachAng)

% e.g:  [C,depth,hs,hm,Tm,Tp,Er,X1,Y1]=get_wave_parameters('S_SNT_2_201312051400.jpg','../RectSNT.mat',0.5,[229 1500],[180 5],[1355261.683 303848.647 8], 1);
%Inputs
% Imgfile : image name with full path (absolute/relative) (e.g. 'S_XX.jpg')
% [nt nx nc]; beach on left, waves coming from left upper corner
% Rectfile: rectify file (.mat)
% dt % temporal resolution (in s)
% sk1% coordinates in pixels of stack extreme point 1 [vertical horizontal] 
% sk2%  coordinates in pixels of stack extreme point 1 [vertical horizontal]
% Er % Dissipated energy (in joule)
% CoordCam [x0,y0,z0]: Video_System_Coordinates (same system as in Rectfile)
% BeachAng : Beach Angle from North direction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Outputs
% C : celerity 
% depth : water depth
% hs: significant wave height
% hm : mean wave height 
% Tm: mean wave period
% Tp: peak period
% Er: dissipated energy
% X1: cross-shore coordinate of the stack
% Y1: longshore coordinate of the stack
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Img=imread(Imgfile);
S=double(Img(:,:,1));
% dt temporelle resolution (in s)
%Method for cross-correlation (1==Dt fixé,2==Dx fixé)
%Stack extrema coordinates
% sk1=[229 1500];
% sk2=[180 5];

%Def des coordonnées des stacks
lmax=max(abs(sk2(2)-sk1(2)),abs(sk2(1)-sk1(1)));
pas1=(sk2(1)-sk1(1))/lmax;
pas2=(sk2(2)-sk1(2))/lmax;
s1=round(sk1(1):pas1:sk2(1));
s2=round(sk1(2):pas2:sk2(2));

%resolution (m/pix along the stack)
[Xv,Yv]=PixtoCoord(Rectfile,s2,s1);

% Coordonnées video
ang=90+BeachAng;%90+angle from north direction
ang=ang*pi/180;
% z0=8;%elevation of the video camera
% y0=303848.647;%Easting of the video camera
% x0=1355261.683;%Northing of the video camera
x0=CoordCam(1);
y0=CoordCam(2);
z0=CoordCam(3);

Co1=Xv;Co2=Yv;
nX=nanmedian(Co1);nY=nanmedian(Co2);
if abs(nanmedian(Co1)-x0)>abs(nanmedian(Co2)-x0)
    Co1=Yv;Co2=Xv;
end


        % Rotation et translation dans le nouveau repère (x0;y0) (pied de
        % l'échafaudage = référence)
        X1=(Co1-x0)*cos(ang)+(Co2-y0)*sin(ang);
        Y1=(Co2-y0)*cos(ang)-(Co1-x0)*sin(ang);


if abs(max(X1)-min(X1))<abs(max(Y1)-min(Y1))
   Xt=X1; 
   X1=Y1;
   Y1=Xt;
end
if abs(X1(end))>abs(X1(1))
X1=X1(end:-1:1);
Y1=Y1(end:-1:1);
end


%resolution along the profile (dx)
dx=abs(sqrt(diff(X1).^2+diff(Y1).^2));dx=smooth(smooth(smooth(dx,20),50),100);
% dx=dx.*cosd(atand(abs(X1(end)-X1(1))/abs(Y1(end)-Y1(1))));



[PosX,PosT,Lw,B,Breakstd,Breakmean1,Breakmean2]=RollerPropertiesTaller(S,dt);

round(nanmedian(PosX));
%Timeseries of intensity
I=S(:,round(nanmedian(PosX)));% second argument is the location (in pixel) of period extraction
[Tm,Tp]=Get_Periode(I,dt);

 save('tot')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Cross-correlation computation: extraction of wave celerity
tic
C = smooth(abs(RadonCIndiv_20140903(dt,dx,B,Tm)),3)';
 
%toc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%Compute wave height and show results
[hs,hm]=BreakerHeight(B,PosT,PosX,Lw,X1,Y1,dx,z0);
% 
% disp(['Hs (1/3) = ',num2str(hs, '%6.2f'), ' m, Hm = ',num2str(hm, '%6.2f'), ' m'])
% disp(['Tp (1/3) = ',num2str(Tp, '%6.2f'), ' s, Tm = ',num2str(Tm, '%6.2f'), ' s'])
% 

%Energie dissipée
rho=1030;
Lwm=nanmean(Lw);

%Nombre de vagues deferlées
nbwave=nansum(Lw./Lw);

%re-interpolation
Lwm=interp1(1:length(Lwm),Lwm,1:length(X1));
nbwave=interp1(1:length(nbwave),nbwave,1:length(X1));
C=interp1(1:length(dx),C,1:length(X1));
dx=interp1(1:length(dx),dx,1:length(X1));

%Water depth along the profile h (in m)
% [depth,ct]=LinearC(Tm,C,0.01);
depth=C.^2/9.81;
disp(['Depth (Breakpoint) = ',num2str(depth((round(nanmean(PosX)))), '%6.2f'), ' m, at X = ',num2str(round(X1((round(nanmean(PosX))))), '%6.2f'), ' m'])


rolL=dx.*Lwm;% en m    
Er=abs(rho*rolL.*C./(2/Tm));% energie dissipée en joule

end



function [trms,Tp]=Get_Periode(S,dt)
 Tcoupure=1.5; %frecuencia de corte
 fr=1/dt;
 Val=(1/Tcoupure)*2*(1/fr);
 ord=1000; fil=fir1(ord,Val,'low');
 kk1=conv(fil,S);
 S= detrend(kk1(  (max(ord)/2)+1: length(kk1)-(max(ord)/2)));

 
 Tcoupure=20; %frecuencia de corte
 fr=1/dt;
 Val=(1/Tcoupure)*2*(1/fr);
 ord=1000; fil=fir1(ord,Val,'high');
 kk1=conv(fil,S);
 S= detrend(kk1(  (max(ord)/2)+1: length(kk1)-(max(ord)/2)));

 
%FFT en précisant la taille du vecteur
Pv = pwelch(detrend(S-mean(S)));
[valmx indmx]=max(Pv);
Tp=1./(indmx./(dt*length(Pv)));

%Wave_Char: computation of waves characteristics
%USAGE
% S;% : vecteur de serie temporelle
%dt : pas de temps
filt=0;% : filtrage (1) /pas de filtre (0)
%meth : Mean zero crossing
meth=1;% mean zero crossing
[hs,htiers,hrms,trms,Hmax,h,Tp,t]=Wave_Char(S,dt,filt,meth);

try
Tp=Tp(1);
catch
Tp=NaN;
end
end



function [df,ct]=LinearC(T,c,precision) 

for i=1:length(c)
 ct=0;%mise à zero du compteur
 w=2*pi/T;
 k=w/c(i);
 g=9.81;%gravite
 do=1000; %valeur quelconque grande
 d=c(i)^2/g;
 while(abs(do-d)>precision)
 ct=ct+1;
do=d;
dispe=w^2-g*k*tanh(k*d);
fdispe=-g*(k^2)./(cosh(k*d)^2);
d=d-dispe./fdispe;
 end
df(i)=d;
 
end
end


function [B]=ImagePreProcessing_20090121Taller(A,dt)
%size(A)= n*m
[nt nc ncol]=size(A);
% %Temporal resolution
% rest=2;
%Spatial resolution
% resc=8;
resc=round(nc./100);

clear B
%B=A(:,:,1);
A3=A;

%tic
for ic=1:resc:nc
    
 Tcoupure=1.5; %frecuencia de corte
 fr=1/dt;
 Val=(1/Tcoupure)*2*(1/fr);
 ord=1000; fil=fir1(ord,Val,'low');
 kk1=conv(fil,A3(:,ic));
 S= detrend(kk1(  (max(ord)/2)+1: length(kk1)-(max(ord)/2)));

 
 Tcoupure=20; %frecuencia de corte
 fr=1/dt;
 Val=(1/Tcoupure)*2*(1/fr);
 ord=1000; fil=fir1(ord,Val,'high');
 kk1=conv(fil,S);
 y= detrend(kk1(  (max(ord)/2)+1: length(kk1)-(max(ord)/2)));



 B(:,ic)=y;%Norma(y,2);
end

ii=find(abs(mean(B))>0&isnan(abs(mean(B)))==0);
%Interpolation spatiale entre icmin et icmax
for irt=1:nt
    try
     B(irt,1:nc)=interp1(ii,B(irt,ii),1:nc,'linear');
    end
end
B(isnan(B)==1)=0;

[B]=RadonSeparationmodif(B);


%disp('Pre-Traitement ok : ')
%toc
%  figure(52);
%  clf
% subplot(4,1,1);pcolor(B(1:40,:));shading flat
% % subplot(4,1,2);pcolor(A3);shading flat
%  subplot(4,1,3);plot(double(A(:,1600)));set(gca,'ydir','reverse')
%  subplot(4,1,4);plot(B(:,1600));set(gca,'ydir','reverse')
end


function [lmval,indd]=lmax(xx,filt)
%LMAX 	[lmval, indd]=lmax(xx,filt). Find local maxima in vector XX,where
%	LMVAL is the output vector with maxima values, INDD  is the 
%	corresponding indexes, FILT is the number of passes of the small
%	running average filter in order to get rid of small peaks.  Default
%	value FILT =0 (no filtering). FILT in the range from 1 to 3 is 
%	usially sufficient to remove most of a small peaks
%	For example:
%	xx=0:0.01:35; y=sin(xx) + cos(xx ./3); 
%	plot(xx,y); grid; hold on;
%	[b,a]=lmax(y,2)
%	 plot(xx(a),y(a),'r+')
%	see also LMIN, MAX, MIN
	
%**************************************************|
% 	Serge Koptenko, Guigne International Ltd., |
%	phone (709)895-3819, fax (709)895-3822     |
%--------------06/03/97----------------------------|

x=xx;
len_x = length(x);
	fltr=[1 1 1]/3;
  if nargin <2, filt=0; 
	else
x1=x(1); x2=x(len_x); 
	for jj=1:filt,
	c=conv(fltr,x);
	x=c(2:len_x+1);
	x(1)=x1;  
        x(len_x)=x2; 
	end
  end
lmval=[]; indd=[];
i=2;		% start at second data point in time series
    while i < len_x-1,
	if x(i) > x(i-1)
	   if x(i) > x(i+1)	% definite max
lmval =[lmval x(i)];
indd = [ indd i];
	   elseif x(i)==x(i+1)&x(i)==x(i+2)	% 'long' flat spot
%lmval =[lmval x(i)];  	%1   comment these two lines for strict case
%indd = [ indd i];	%2 when only  definite max included
i = i + 2;  		% skip 2 points
	   elseif x(i)==x(i+1)	% 'short' flat spot
%lmval =[lmval x(i)];	%1   comment these two lines for strict case
%indd = [ indd i];	%2 when only  definite max included
i = i + 1;		% skip one point
	   end
	end
	i = i + 1;
    end
if filt>0 & ~isempty(indd),
	if (indd(1)<= 3)|(indd(length(indd))+2>length(xx)), 
	   rng=1;	%check if index too close to the edge
	else rng=2;
	end
	  for ii=1:length(indd), 	% Find the real maximum value
	    [val(ii) iind(ii)] = max(xx(indd(ii) -rng:indd(ii) +rng));
	    iind(ii)=indd(ii) + iind(ii)  -rng-1;
	  end
  indd=iind; lmval=val;
else
end
end


function [lmval,indd]=lmin(xx,filt)
%LMIN 	function [lmval,indd]=lmin(x,filt)
%	Find local minima in vector X, where LMVAL is the output
%	vector with minima values, INDD is the corresponding indeces 
%	FILT is the number of passes of the small running average filter
%	in order to get rid of small peaks.  Default value FILT =0 (no
%	filtering). FILT in the range from 1 to 3 is usially sufficient to 
%	remove most of a small peaks
%	Examples:
%	xx=0:0.01:35; y=sin(xx) + cos(xx ./3); 
%	plot(xx,y); grid; hold on;
%	[a b]=lmin(y,2)
%	 plot(xx(a),y(a),'r+')
%	see also LMAX, MAX, MIN
	
%
%**************************************************|
% 	Serge Koptenko, Guigne International Ltd., |
%	phone (709)895-3819, fax (709)895-3822     |
%--------------06/03/97----------------------------|

x=xx;
len_x = length(x);
	fltr=[1 1 1]/3;
  if nargin <2, filt=0; 
	else
x1=x(1); x2=x(len_x); 

	for jj=1:filt,
	c=conv(fltr,x);
	x=c(2:len_x+1);
	x(1)=x1;  
        x(len_x)=x2; 
	end
  end

lmval=[];
indd=[];
i=2;		% start at second data point in time series

    while i < len_x-1,
	if x(i) < x(i-1)
	   if x(i) < x(i+1)	% definite min
lmval =[lmval x(i)];
indd = [ indd i];

	   elseif x(i)==x(i+1)&x(i)==x(i+2)	% 'long' flat spot
%lmval =[lmval x(i)];	%1   comment these two lines for strict case 
%indd = [ indd i];	%2 when only  definite min included
i = i + 2;  		% skip 2 points

	   elseif x(i)==x(i+1)	% 'short' flat spot
%lmval =[lmval x(i)];	%1   comment these two lines for strict case
%indd = [ indd i];	%2 when only  definite min included
i = i + 1;		% skip one point
	   end
	end
	i = i + 1;
    end

if filt>0 & ~isempty(indd),
	if (indd(1)<= 3)|(indd(length(indd))+2>length(xx)), 
	   rng=1;	%check if index too close to the edge
	else rng=2;
	end

	   for ii=1:length(indd), 
		[val(ii) iind(ii)] = min(xx(indd(ii) -rng:indd(ii) +rng));
		iind(ii)=indd(ii) + iind(ii)  -rng-1;
	   end
  indd=iind; lmval=val;
else
end
end


function  fsig=FiltreMean(a,b)

version=1;

[n1 n2]=size(a);
if n2==1
    a=a';
end


%Version 1

if(version==1)
la=length(a);
if(la<(2*b))
disp('Attention : vecteur trop petit pour filtrage')
end
fsig(1:la)=0;

fsig2(1:la+2*b)=0;
fsig2(1:b)=mean(a(1:b));

fsig2(b+1:b+la)=a;

fsig2(b+la+1:length(fsig2))=mean(a(la-b:la));



k=0;
    for kk=b+1:b+la
        k=k+1;
    fsig(k)=mean(fsig2(kk-b:kk+b));
    
%     figure(41);clf;hold on;plot(fsig2,'k')
%     hold on;plot(kk-b:kk+b,fsig2(kk-b:kk+b),'r','LineWidth',2)
%     pause
    
    end
%    fsig(la-(b-1):la)=(a(la-1)+a(la-(b-1)))/2;
%    fsig(1:(b-1))=(a(1)+a((b-1)))/2;

elseif(version==2)
    % Version matricielle
    a=double(a);
    
    a2=a';
    
am(1:b)=a(1);
ap(1:b)=a(length(a));
a2=[am' ;a'; ap'];
a3(1:length(a))=a';
c=0;
for i=1:b
a3=a3+a2((b+1:(length(a2))-b)+i)'+a2((b+1:(length(a2))-b)-i)';
c=c+1;
end
fsig=a3/(2*b+1);


end % Version


%On remplace les valeurs NaN (aux bornes) par les valeurs du vecteur
%original
fsig(find(isnan(fsig)==1))=a(find(isnan(fsig)==1));

end    


function [PosX,PosT,Lw,B,Breakstd,Breakmean1,Breakmean2]=RollerPropertiesTaller(S,dt)

I=double(S);
[nc nt nv]=size(S');

 [B]=ImagePreProcessing_20090121Taller(I(:,:,1),dt);
save('totr')


Breakstd=FiltreMean(FiltreMean(nanstd(B(round(nt/20:19*nt/20),:)),round(nc./50)),round(nc./30));
 %Computation of threshold (thresh)
  thresh=2.*nanmax(Breakstd)./3;
Breakstd0=(Breakstd-nanmin(Breakstd))./nanmax(Breakstd-nanmin(Breakstd));
Breakmean1=FiltreMean(nanmax(I),round(nc./50))-FiltreMean(nanmean(I),round(nc./50));Breakmean1=(Breakmean1-nanmin(Breakmean1))./nanmax(Breakmean1-nanmin(Breakmean1));
Breakmean2=FiltreMean(nanmean(nanmean(S,3)),round(nc./50));Breakmean2=(Breakmean2-nanmin(Breakmean2))./nanmax(Breakmean2-nanmin(Breakmean2));



  
      %Normalisation de std
Breakstd=-diff([Breakstd(1) Breakstd]);
Breakstd=(Breakstd-nanmin(Breakstd))./nanmax(Breakstd-nanmin(Breakstd));
B=B.*repmat(Breakstd,nt,1);
B(round([1:5 nt-5:nt]),:)=0;

   B2=zeros(size(B));
   B2(find(B>thresh))=1;
B2(round([1:5 nt-5:nt]),:)=0;
 [nt nc]=size(B2);
 
 
%On recherche les pixels correspondant aux rollers a chaque instant
clear indw Lw
Lw=zeros(nt,nc).*NaN;
for i=nt:-1:1
 %Indices des pixels de roller
   ind= find(B2(i,:)==1);

    if( length(ind)>2)
        
        %on recherche si plusieurs rollers
        fw=find(diff(ind)>10);
         nbw=length(fw)+1;
       
        %On calcule la longueur de chaque roller: limites des rollers
            indw=sort([ind(1) ind(fw) ind(fw+1) ind(length(ind))]);
         
        %Calcul des longueurs
        for v=1:nbw
            Lw(i,round(mean(indw(2*(v-1)+1): indw(2*v))))=length(indw(2*(v-1)+1): indw(2*v));
        end
    end
end
  Lw(find(Lw==0))=NaN;

%Recherche de la position des points de deferlement
clear Def D PosT PosX PosXi
for i=1:nt
    try
    Def(i)=min(find(Lw(i,:)>0));
    catch
    Def(i)=NaN;
        
    end
end

%PosX: position des points de deferlement
% oo=find(abs(Def-nanmedian(Def))<nanstd(Def));
% D=interp1(oo,Def(oo),1:nt);
% [PosXi PosTi]=lmin((FiltreMean(D,1)));
% 
% %Positions sans lissage
% [PosX1 PosT1]=lmin((FiltreMean(D,1)));
% 
% %On reaffine la position et le temps des points de deferlement
% for i=1:length(PosXi)
%     try
%     [az rgt]=min(abs( PosTi(i)-PosT1));
%     PosX(i)=PosX1(rgt);
%     PosT(i)=PosT1(rgt);
%     end
% end
%Version 2
oi=find(isnan(Def)==0);
oo=find(abs(diff(Def(oi)))>10);
oo=oi([1 oo+1]);

PosX=Def(oo);
PosT=oo;

try
PosX=round(PosX);

vec=0.5*(Breakstd0+Breakstd);
% On cherche les points de deferlement prochent du max de variation d'intensité de derlement 
seuil=min(vec)+2*(max(vec)-min(vec))./3;
io=find(vec(PosX)>seuil);
PosX=PosX(io);
PosT=PosT(io);
catch
    PosX=NaN;
    PosT=NaN;
end


for i=1:length(PosX)
    for k=find(1:length(PosX)~=i)
        try
    ix=round(linspace(PosX(i),PosX(k),20));
    it=round(linspace(PosT(i),PosT(k),20));
   if length(find(diag(B(it,ix))<0))<length(ix)/5&PosT(i)<PosT(k);PosT(k)=NaN;PosX(k)=NaN;end
        end
    end
end

PosX=PosX(isnan(PosT)==0);
PosT=PosT(isnan(PosT)==0);

% % On corrige les rollers detectées au dela du point de deferlement
% clear indw Lw
% Lw=zeros(nt,nc).*NaN;
% for i=nt:-1:1
%  %Indices des pixels de roller
%    ind= find(B2(i,:)==1);
%     if length(ind)>2
%         %on recherche si plusieurs rollers
%         fw=find(diff(ind)>10);
%          nbw=length(fw)+1;
%        
%         %On calcule la longueur de chaque roller: limites des rollers
%             indw=sort([ind(1) ind(fw) ind(fw+1) ind(length(ind))]);
%          
%         %Calcul des longueurs
%         for v=1:nbw
%             if abs(indw(2*v)-nanmedian(PosX))<nanstd(PosX)%Si on est proche du point de deferlement moyen
%             Lw(i,round(mean(indw(2*(v-1)+1): indw(2*v))))=length(indw(2*(v-1)+1): indw(2*v));
%             end
%         end
%     end
% end
%   Lw(find(Lw==0))=NaN;






%Enregistrement des parametres moyens
%Position du point de deferlement
% PosM=mean(PosX);
% 
% %Longueur moyenne du roller :LwM
% clear LwM1 LwM
% for i=1:nc
%     LwM1(i)=mean(Lw(find(isnan(Lw(:,i))==0),i));    
% end
% e=find(isnan(LwM1)==0);
% LwM=interp1(e,LwM1(e),1:nc);
% 
% 
% %Affichage

% figure(89);clf;pcolor(Lw);shading flat
% h=colorbar;
% ylabel(h,'Roller length (pix)')
% hold on;plot(PosX,PosT,'o')
% set(gcf,'Color','w')
% xlabel 'distance (pix)'
% ylabel 'time (pix)'
% box on
% title 'Roller characteristics: length and inception'
end


function [hs,htiers,hrms,trms,Hmax,h,Tp,tindiv]=Wave_Char(d,dt,filt,meth)
%Calcul les caractéristiques des vagues à partir d'une série temporelle de
%hauteur d'eau (ou équivalent)
%Inputs
%d : vecteur de serie temporelle
%dt : pas de temps
%filt : filtrage (1) /pas de filtre (0)
%meth : Mean zero crossing
% meth=1 mean zero crossing
% meth=2 up zero crossing
% meth=3 down zero crossing

%Outputs
%h : hauteurs de vague individuelles (ex: pour faire un histogramme)
%tindiv: periodes de vagues individuelles (ex: pour faire un histogramme)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%INITIALISATION%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Chargement du fichier de données :
 %d=load('D:\Documents\THESE\DONNEES_INSITU\DGO41001.txt');
% Numero de la colonne du capteur dans le fichier (ex:10)
  %nc=5;
 % Pas de temps (dt) des mesures des capteurs
 %dt = 0.5; % secondes
 
 %Nombre de pas de temps de l'échantillon (nt)
 [n1 n2]=size(d);
if n1==1
    d=d';
end
 
 
 
 
 sz=size(d);
 nt=sz(1);
 zc=meth;
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %%%%%%%%%%%%%%Fin initialisation%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
% Chargement de la colonne du capteur à analyser :
% ici colonne 10
% durée de donnée à analyser : environ 10 min (soit 1820 mesures)
 k=d(1:nt);
% On detrend le signal brut
%generation d'un vecteur pour les abscisses
 x=[1:nt];
 a=reshape(x,nt,1);
 
 % generation des coefficients du polynome de degré un
  p=polyfit(a,k,1);
 
 %Suppression de la partie trendée de la courbe
 
 k2=k-p(1).*a;
 
 % Suppression de la partie moyenne
 km=k2;
 k2=km-mean(km);
 
 clear k4
if(filt==1)
freq = (1./dt)*(1:nt)/nt;
 %frequences de coupure : %%%%%%%%%%%%%%%%%%%%%
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %fréquence de coupure vers les basses fréqences
 % fcb= 1/Tlimit
  fcb=0.01;
 % fréquence de coupure vers les hautes fréquences
 % fch=1/Tlimit
  fch=1/0.1;
f1 = fft(k2);
Indb = find(freq < fcb);
Indh = find(freq > fch);
f1(Indb) = 0.0;
Sinv = ifft(f1);
k3 = Sinv.*conj(Sinv);
f1(Indh) = 0.0;
Sinv = ifft(f1);
k4 = real(Sinv);
k2=k4;
  k3=FiltreMean(k2,round(dt));
 else%filt==0
     
     

%      k3=k2'-FiltreMean(k2',round(10/(2*dt)));
     k3=FiltreMean(k2,round(dt));
     k4=k2;
     

     
 end %filt
%  size(k2)
%  figure;plot(k3)
%  round(5/(2*dt))
 %k4 est un signal filtré pour les hautes et basses fréquences
 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Mesures de la physique des vagues%%%%%%%%%%%%%%
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
 
 % Localisation des zero du signal :
 clear z j
 j=0;
z(1:floor(nt./2))=0;

%Mean zero crossing
% zc=1 mean zero crossing
% zc=2 up zero crossing
% zc=3 down zero crossing
% for zc=1:2

if(zc==1)
 for i=2:nt-1
     if(k3(i)*k3(i+1)<=0 && k3(i-1)*k3(i)>=0 ) 
         j=j+1;
        z(j)=i;
    end
 end
 
%  figure(36);plot(k2,'k');hold on;plot(k3,'r');plot(z,0,'r*')
elseif(zc==2)
    clear z j
     j=0;
 for i=2:nt-1
     if(k3(i)*k3(i+1)<=0 && mean([k3(i-1) k3(i)])<0 ) 
         j=j+1;
        z(j)=i;
    end
 end   
elseif(zc==3)
     clear z j
     j=0;
 for i=2:nt-1
     if(k3(i)*k3(i+1)<=0 && mean([k3(i-1) k3(i)])>0 ) 
         j=j+1;
        z(j)=i;
        mean([k3(i-1) k3(i)])
    end
 end    
end

%calcul de la taille du vecteur dans lequel sont stockés les indices des zeros
k=0;
sz=size(z);
for i=1:sz(2)
    if(z(i)~=0) 
        k=k+1;
    end
end

% Recherche du max de la hauteur de la vague entre deux zeros

hmax(1:k-1)=0;
hmin(1:k-1)=0;

for i=1:k-1
    if(mean(k4(z(i):z(i+1)))>0.)
        if(max(k3(z(i):z(i+1)))>0)
    hmax(i)=max(k4(z(i):z(i+1)));
        else
    hmax(i)=0;
        end
    hmin(i)=0.;
    
else
    hmax(i)=0.;
        if(min(k3(z(i):z(i+1)))<0)
        hmin(i)=min(k4(z(i):z(i+1)));
    
        else
        hmin(i)=0.; 
        end    
    end
end

% calcul de la hauteur de chaque vague passant par le capteur
h=[];
kh=0;
h(1:k-2)=0;

for i=1:k-2
    
 if(hmin(i)~=0 & hmax(i+1)~=0)
        if(hmin(i+1)==0 & hmax(i)==0)        
        kh=kh+1;
     h(kh)=hmax(i+1)-hmin(i);
        end
 end
end

% Calcul de la hauteur rms moyenne des vagues sur la durée de mesure
%hrms1=0;
clear hrms ord htiers

[sd ord]=sort(h(find(h~=0)),'descend');

%Nombre de vagues pour le calcul
% disp(['Nombre de vagues : ',num2str(length(find(h~=0)))])

% Affichage d'un Histogramme
% figure(85);clf
% bar(h(find(h~=0)));
% set(gcf,'Color','w')
% set(gca,'FontSize',14)
% xlabel('H (m)')
% ylabel ('N°')
% grid on
% [val,nval]=hist(h(find(h~=0)));
% figure(86)
% plot(nval,FiltreMean(val,2),'k')
% set(gcf,'Color','w')
% set(gca,'FontSize',14)
% xlabel('H (m)')
% ylabel ('N°')
% grid on

%Hauteur max
Hmax=max(h);

%Hauteur 1/3
htiers=mean(h(ord(1:round(length(ord)/3))));% on prend le premier tiers pour la hauteur rms
%Hauteur rms
hrms=mean(h(1:kh));

%Hauteur significative
hs=4.*std(k4);
% disp(['Hs=',num2str(hs),'; Hrms=',num2str(hrms),'; H1/3=',num2str(htiers),'; Hmax=',num2str(Hmax)])
%Calcul de la periode de la houle pour chaque vague

clear t
t(1:k-1)=0;
for i=1:k-1
    if(z(i)~=0 | z(i+1)~=0 ) 
        t(i)=2*(z(i+1)-z(i))*dt;
    end
end

tindiv=t;
%Calcul de la periode rms de la houle et affichage

if(zc==1)
trms=mean(t);
elseif(zc==2)
trms=mean(t)/2;
elseif(zc==3)
trms=mean(t)/2;
end    


%%%%%Periode Peak%%%%%%%%%%%
y=detrend(d);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Nombre de points pour la fft
nx=length(y);
%Pas de temps
% x=(1:nx);
%Resolution temporelle (periode d'échantillonnage) (en s)
% dt=1;

%FFT en précisant la taille du vecteur
z=fft(y);
[n1 n2]=size(z);
if n1==1
    z=z';
end

t=1./((1:nx/2)./(dt*nx));
in=find(t<6|t>14);
x=1:nx/2;
k=z(x).*conj(z(x)/nx);
k(in)=0;

[sa ind]=lmax(FiltreMean(FiltreMean(k,3),5));


%  disp(ind./(dt*nx))

% io=find(sa>max(sa)/2);
Tp=1./(ind./(dt*nx));
end


function [hs,hm]=BreakerHeight(B,PosT,PosX,Lw,X1,Y1,dx,z0)
%Calcul de la hauteur des vagues au point de deferlement
%Principe: On prend la distance (horizontale) entre le maximum du point de deferl 
%et le minimum le precedent. On reporte ensuite cette distance en verticale
%connaissant la hauteur relative de la camera / surface de l'eau
%PosT et PosX sont les positions spacio-temporelles des points de
%deferlement
%dx: resolution spatiale
%AngleCam angle avec la camera
save('toth')

try
AngleCam=abs(z0./X1(PosX));%(Rad) (Hcam-Hmean water level/Xcam-Xbreak)
AngleWaveFront=35*pi/180;%(Rad) (or 30°);Wave face angle at breaking inception
dx=abs(dx);
clear L1;
Lwi=nanmax(Lw);Lwi=smooth(interp1(find(isnan(Lwi)==0&Lwi>0),Lwi(find(isnan(Lwi)==0&Lwi>0)),1:length(dx)),80);

for i=1:length(PosT)
    try
%Minimum avant variation (le min le plus proche du
%deferl (max d'intensite)

% vec=FiltreMean(FiltreMean(B(PosT(i),PosX(i):end),5),2);
% %version 1
% [a b]=lmin(vec);b=b(1);

vec=FiltreMean(FiltreMean(nanmax(B(PosT(i)-25:PosT(i)+25,PosX(i)-50:end))-nanmin(B(PosT(i)-25:PosT(i)+25,PosX(i)-50:end)),20),5);
[ngt gft]=nanmax(vec(1:round(end/3)));
% figure(158);clf;plot(vec)

ii=find(vec>nanmin(vec)+0.5*(nanmax(vec)-nanmin(vec)));id=[ii(1) ii(find(diff(ii)>20)) ii(end)];
idl=find(id<gft);idl=idl(end);
idp=find(id>gft);idp=idp(1);
b=length(id(idl):id(idp));
% hold on;plot(id(idl):id(idp),vec(id(idl):id(idp)),'r*')
%  pause

% version2
% b=Lwi(PosX(i));

L1(i)=abs(b).*dx(PosX(i))';
    catch
        L1(i)=NaN;
    end
end
cor=(L1).*tan(AngleCam)/tan(AngleWaveFront);
Lf=(L1-cor).*tan(AngleCam);



ind=find(Lf>0&isnan(Lf)==0);%Only positive wave heights


%On calcul la hauteur des vagues
Lord=sort(Lf(ind),'descend');Lord=Lord(round(max([2 end/10]):min([end-1 9*end/10])));
%Nombre de Vagues sur la duree
nbv=length(ind);
%Hauteur 1/3
hs=nanmedian(Lord(1:round(length(Lord)/3)));% on prend le premier tiers pour la hauteur rms
%Hauteur rms
hm=nanmedian(Lord);


vecm=0;
for i=1:length(PosT)
    try
    vec=FiltreMean(FiltreMean(nanmax(B(PosT(i)-25:PosT(i)+25,PosX(i)-100:end))-nanmin(B(PosT(i)-25:PosT(i)+25,PosX(i)-100:end)),20),5);
    vecm=vecm+vec(1:min([end 350]));
%     figure(147);hold on;plot(vec,'k')
    end
end

[erd gtf]=max(diff(vecm));
[efrd gtd]=min(diff(vecm));
L1=dx(round(nanmedian(PosX))).*abs(gtd-gtf);

% version 1 
% seuil=find(vecm(50:end)>min(vecm)+0.75*(max(vecm)-min(vecm)));ib=find(diff(seuil)>1);
% if length(ib)>0
% L1=dx(round(nanmedian(PosX))).*length(seuil(1:ib(1)));
% else
% L1=dx(round(nanmedian(PosX))).*length(seuil(1:end));
% end




cor=(L1).*tan(nanmean(AngleCam))/tan(AngleWaveFront);
Lf=(L1-cor).*tan(nanmean(AngleCam));
hm=Lf;

% hm=hm.*cosd(atand(nanmean(diff(X1))./nanmean(diff(Y1))));
% hs=hs.*cosd(atand(nanmean(diff(X1))./nanmean(diff(Y1))));

catch
    hm=NaN;
    hs=NaN;
end

end


function [CRadmoy]=RadonCIndiv_20140903(dt,dx,In,Tm)
%Inputs%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%In: Input spatio-temporal elevation matrix, format: (nt,nx)
% In=ZHF;
%pdx: Spatial resolution degradation (in point)
pdx=2;
%pdt: Temporal resolution degradation (in point)
pdt=1;
%Wx: Window size in space for Radon computation (sensity test : better if WX > 10)
% Wx= 30; %
%freq_t : temporal frequency (1/dt) (taking into account for pdt)
freq_t = 1./dt;%
%freq_x : spatial frequency (1/dx)
freq_x = 1./(dx.*pdx);%(here 2 pixel/m)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Outputs%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% an: sorties brutes radon (angle en °) mieux si entre 30 et 70
% Hx: élévation des crêtes détectées pour le calcul de la célérité
% Hm: élévation des creux détectées pour le calcul de la célérité
% Tim: crêtes détectées pour le calcul de la célérité
% C: célérités des crêtes détectées
%Outputs%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 save('totc')

M=detrend(abs(diff(In(1:pdt:length(In(:,1)),1:pdx:length(In(1,:))),1))')';
% M=In(1:pdt:length(In(:,1)),1:pdx:length(In(1,:)));

Wx=round(size(M,2)./10);


clear Mat VER XHFR an Tim Hx Hm C CC2 Tan CRadmoyt CRadmoy
for ix=Wx+1:Wx-1:size(M,2)-(Wx+1)%Boucle sur les positions X

MR=detrend(M(1:length(M(:,1)),ix-Wx:ix+Wx));%selection de la fenetre

nt=size(MR,1);
R=radon(MR,1:1:180);%transformée de Radon
% figure(128);clf;pcolor(R);shading flat
% pause
nr=size(R,1);
amp=nr/nt;
iang=1:1:180;

k=nt;%-round(nt./2);
% trk=round(nr/2+tr*((k-(nr+1)/2))*amp);
trk=floor((size(MR,1))/2)-floor((  0*cosd(iang)+  ((size(MR,1))/2-k*amp)*sind(iang)));trk=trk-min(trk);
res=(nt*dt)./(trk.*2);

R2=R;
for i=iang
   R2(:,i)=smooth(R(:,i)-smooth(R(:,i),round(1+(Tm+2)./(res(i)))),round(max([1 (Tm-2)./res(i)]))); 
end

vec=smooth(nanstd(R2(round(size(R2,1)/4:3*size(R2,1)./4),1:90)),10);
[frd a2]=max(vec);%calcul celerite moyenne

if vec(a2)>1.15*vec(90)
if length(freq_x)==1
C2=(1/mean(freq_x))/(tand(90-a2)*(1/mean(freq_t)));% celerite moyenne
else
C2=(1/nanmedian(freq_x(round(pdx.*(ix-Wx:ix+Wx)))))/(tand(90-a2)*(1/mean(freq_t)));% celerite moyenne
end
else
C2=NaN;
end


 
CRadmoyt(ix-Wx:ix+Wx)=C2; %celerite moyenne
end
notnan=find(isnan(CRadmoyt)==0&CRadmoyt>0);
CRadmoy=dx*NaN;
try
CRadmoy(min(notnan*pdx):max(notnan*pdx))=interp1(pdx.*notnan,CRadmoyt(notnan),min(notnan*pdx):max(notnan*pdx));
end
end


function [Vcoord1,Vcoord2]=PixtoCoord(Rectfile,Vpix1,Vpix2)

clear Vcoord1 Vcoord2
load(Rectfile);
[X_I] = rectify0707(double([Vpix1; Vpix2]),Rckk,Tckk,fc,cc,kc,alpha_c,mean(X_kk(3,:)));

Vcoord1=X_I(1,:);
Vcoord2=X_I(2,:);
end

function [X_r] = rectify0707(aaa,Rckk,Tckk,fc,cc,kc,alpha_c,Z);

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
xn = normalize(aaa,fc,cc,kc,alpha_c); 


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

end



function [omckk,Tckk,Rckk,H,x,ex,JJ] = compute_extrinsic(x_kk,X_kk,fc,cc,kc,alpha_c,MaxIter,thresh_cond)

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

[omckk,Tckk,Rckk] = compute_extrinsic_init(x_kk,X_kk,fc,cc,kc,alpha_c);

% Refinement:
[omckk,Tckk,Rckk,JJ] = compute_extrinsic_refine(omckk,Tckk,x_kk,X_kk,fc,cc,kc,alpha_c,MaxIter,thresh_cond);


% computation of the homography (not useful in the end)

H = [Rckk(:,1:2) Tckk];

% Computes the reprojection error in pixels:

x = project_points2(X_kk,omckk,Tckk,fc,cc,kc,alpha_c);

ex = x_kk - x;


% Converts the homography in pixel units:

KK = [fc(1) alpha_c*fc(1) cc(1);0 fc(2) cc(2); 0 0 1];

H = KK*H;


end


function [omckk,Tckk,Rckk] = compute_extrinsic_init(x_kk,X_kk,fc,cc,kc,alpha_c),

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

xn = normalize(x_kk,fc,cc,kc,alpha_c);



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
   
   H = compute_homography(xn,X_new(1:2,:));
   
   % De-embed the motion parameters from the homography:
   
   sc = mean([norm(H(:,1));norm(H(:,2))]);
   
   H = H/sc;
   
   % Extra normalization for some reasons...
   %H(:,1) = H(:,1)/norm(H(:,1));
   %H(:,2) = H(:,2)/norm(H(:,2));
   
   if 0, %%% Some tests for myself... the opposite sign solution leads to negative depth!!!
       
       % Case#1: no opposite sign:
       
       omckk1 = rodrigues([H(:,1:2) cross(H(:,1),H(:,2))]);
       Rckk1 = rodrigues(omckk1);
       Tckk1 = H(:,3);
       
       Hs1 = [Rckk1(:,1:2) Tckk1];
       xn1 = Hs1*[X_new(1:2,:);ones(1,Np)];
       xn1 = [xn1(1,:)./xn1(3,:) ; xn1(2,:)./xn1(3,:)];
       e1 = xn1 - xn;
       
       % Case#2: opposite sign:
       
       omckk2 = rodrigues([-H(:,1:2) cross(H(:,1),H(:,2))]);
       Rckk2 = rodrigues(omckk2);
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
       omckk = rodrigues(RRR);

       %omckk = rodrigues([H(:,1:2) cross(H(:,1),H(:,2))]);
       Rckk = rodrigues(omckk);
       Tckk = H(:,3);
       
   end;
   
      
   
   %If Xc = Rckk * X_new + Tckk, then Xc = Rckk * R_transform * X_kk + Tckk + T_transform
   
   Tckk = Tckk + Rckk* T_transform;
   Rckk = Rckk * R_transform;

   omckk = rodrigues(Rckk);
   Rckk = rodrigues(omckk);
   
   
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
   
	omckk = rodrigues(Rckk);
   Rckk = rodrigues(omckk);
   
end
end

function [omckk,Tckk,Rckk,JJ] = compute_extrinsic_refine(omc_init,Tc_init,x_kk,X_kk,fc,cc,kc,alpha_c,MaxIter,thresh_cond),

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
    
    [x,dxdom,dxdT] = project_points2(X_kk,omckk,Tckk,fc,cc,kc,alpha_c);
    
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

Rckk = rodrigues(omckk);

end




function [xp,dxpdom,dxpdT,dxpdf,dxpdc,dxpdk,dxpdalpha] = project_points2(X,om,T,f,c,k,alpha)

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

[Y,dYdom,dYdT] = rigid_motion(X,om,T);


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

[x,dxdom,dxdT,dxdf,dxdc,dxdk,dxdalpha] = project_points2(X,om,T,f,c,k,alpha);


% Test on om: OK

dom = 0.000000001 * norm(om)*randn(3,1);
om2 = om + dom;

[x2] = project_points2(X,om2,T,f,c,k,alpha);

x_pred = x + reshape(dxdom * dom,2,n);


norm(x2-x)/norm(x2 - x_pred)


% Test on T: OK!!

dT = 0.0001 * norm(T)*randn(3,1);
T2 = T + dT;

[x2] = project_points2(X,om,T2,f,c,k,alpha);

x_pred = x + reshape(dxdT * dT,2,n);


norm(x2-x)/norm(x2 - x_pred)



% Test on f: OK!!

df = 0.001 * norm(f)*randn(2,1);
f2 = f + df;

[x2] = project_points2(X,om,T,f2,c,k,alpha);

x_pred = x + reshape(dxdf * df,2,n);


norm(x2-x)/norm(x2 - x_pred)


% Test on c: OK!!

dc = 0.01 * norm(c)*randn(2,1);
c2 = c + dc;

[x2] = project_points2(X,om,T,f,c2,k,alpha);

x_pred = x + reshape(dxdc * dc,2,n);

norm(x2-x)/norm(x2 - x_pred)

% Test on k: OK!!

dk = 0.001 * norm(k)*randn(5,1);
k2 = k + dk;

[x2] = project_points2(X,om,T,f,c,k2,alpha);

x_pred = x + reshape(dxdk * dk,2,n);

norm(x2-x)/norm(x2 - x_pred)


% Test on alpha: OK!!

dalpha = 0.001 * norm(k)*randn(1,1);
alpha2 = alpha + dalpha;

[x2] = project_points2(X,om,T,f,c,k,alpha2);

x_pred = x + reshape(dxdalpha * dalpha,2,n);

norm(x2-x)/norm(x2 - x_pred)

end




function	[out,dout]=rodrigues(in)

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

[R1,dR1] = rodrigues(om);
R2 = rodrigues(om+dom);

R2a = R1 + reshape(dR1 * dom,3,3);

gain = norm(R2 - R1)/norm(R2 - R2a)

%%% TEST OF dOmdR:
om = randn(3,1);
R = rodrigues(om);
dom = randn(3,1)/10000;
dR = rodrigues(om+dom) - R;

[omc,domdR] = rodrigues(R);
[om2] = rodrigues(R+dR);

om_app = omc + domdR*dR(:);

gain = norm(om2 - omc)/norm(om2 - om_app)


%%% OTHER BUG: (FIXED NOW!!!)

omu = randn(3,1);   
omu = omu/norm(omu)
om = pi*omu;        
[R,dR]= rodrigues(om);
[om2] = rodrigues(R);
[om om2]

%%% NORMAL OPERATION

om = randn(3,1);         
[R,dR]= rodrigues(om);
[om2] = rodrigues(R);
[om om2]

end




function [xn] = normalize(x_kk,fc,cc,kc,alpha_c)

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
	xn = comp_distortion_oulu(x_distort,kc);
else
   xn = x_distort;
end

end



function [fc,cc,kc,alpha_c]=convert_heikkila_stoyanov(par)

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
  sys = [
      2016,     %number of pixels in horizontal direction 753
      720,     %number of pixels in vertical direction 582
      6.4512,     %effective CCD chip size in horizontal direction 
      4.8896,     %effective CCD chip size in vertical direction
      6.0,       %nominal focal length
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


function [Y,dYdom,dYdT] = rigid_motion(X,om,T);

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


[R,dRdom] = rodrigues(om);

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


function [H,Hnorm,inv_Hnorm] = compute_homography(m,M);

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

function [x] = comp_distortion_oulu(xd,k);

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
    
    [x] = comp_distortion(xd,k);
    
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


function [x_comp]  = comp_distortion(x_dist,k2);

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
end



% RadonSeparation.m
function [Sin]=RadonSeparationmodif(M)
%rafael.almar@ird.fr

%Input
%M=Timestack [nt,nx]

%Output
%Sin=Timestack [nt,nx]: incident component
%Sout=Timestack [nt,nx]: reflected component
    
%Matrix dimensions
[nt,nx]=size(M);
if nx>nt
    disp ('Warning !! nx>nt, rotate matrix (time,space) ?')
Sin=M;
else
% Pre-Traitement
for i=1:nx
	M(1:nt,i)=M(1:nt,i)-nanmean(M(1:nt,i));
end

[Sin]=FiltreRadon(M(1:nt,1:nx),1,89);
%disp('Radon ok: composante incidente extraite')
    

end
%petites valeurs: horizontal
%valeurs proches de 90: vertical   
end


function [S2]=FiltreRadon(S,Lm,Lx)
%rafael.almar@ird.fr

%Inputs
%S: Image originale
%Lm: angle min accepte (en deg: 0:179)
%Lx: angle max accepte (en deg: 0:179)
%Outputs
%Sout: Image filtree entre les angles Lm et Lx

%Taille de l'image
[c1 c2 c]=size(S);
    
%Transformation de Radon
R=radon(double(S(:,:,1)),0:179);

%On reconstruit l'image 
I = iradon(R(:,Lm:Lx),(Lm:Lx)-1,'linear','Hann',c1); 

%Indices reellement filtres
mx=min(size(I));
mn=min(size(S));
ind=1+round(mx/2-mn/2):round(mx/2-mn/2)+mn;

%Outputs
S2=S.*NaN;
if c2>c1
S2(:,1+(c2-c1)/2:c2-(c2-c1)/2)=I(:,ind);
else
S2(:,:)=I(:,ind);
    
end
S2=S2*0.5;
end