%DetectSeuil_V20130725_GPP.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 function [x,y]=DetectSeuil_V20130725_GPP(Img,cmin,cmax,nbmin,path)
% Description : Fonction qui calcule la position de la ligne d'eau sur une
% image video non-rectifiée
%tic
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Inputs :
% - Img : image sur laquelle va être calculée la ligne d'eau
%ex:
%Img=imread('E:\TV1\Cam-Era\Camera2\20080308\13/CamB_2008_03_08_13_02_36_93.jpg');
% 
% - Cam : numero de la camera
%ex: pour Truc Vert : Cam=2 ou Cam=1
%
% - Site : nom du site : Truc Vert(TV) ou Biscarrosse (BI)
%ex : Site='TV';
%ex: Site='Sete';
%  - cmin critère (entre pic max et min) pour début recherche de seuil
% ex : 1/5
%  - cmax critère (entre pic max et min) pour début recherche de seuil
% ex : 1/10
%
%  - nbmin nombre de min maximal
% ex : 5
%  - path chemin d'accès aux fichier zone
% ex : 5

% Outputs :
%
% (x2,y2): coordonnées lissées des points de la ligne d'eau calculée (Ne
% pas utiliser - variables vides)
% 
% (x,y) : coordonnées des points de la ligne d'eau calculée (a utiliser)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%Affichage%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Usage : Aff=''yes'/'no'
   Aff1='no';
   Aff='no';
   Aff2='no';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   preTraitement de l'image
%Correction des effets optiques de variation d'intensité sur les bords de l'image
%      a3=IntensityCorrection(Img,Site,Cam);
     a3=Img;
% end
% a3=StretchBand(a3,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%AFFICHAGE%%%%%%%%%%%%%%%%%
% figure(12);set(gcf,'Color','w')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% clf
 
%Paramètre de convergence 1
Conv1=28;
 


bande1=double(a3(:,:,1));%rouge
% bande2=double(a3(:,:,3));%vert
bande2=0.5*double((a3(:,:,2)+a3(:,:,3)));%vert et bleu

%Definition du ration rouge/vert
 rv=bande1./bande2;
 
rv=reshape(bande1./bande2,numel(rv),1);
 
%On enlève les cas de saturation (r=v ou le noir (r=v)) ainsi que les
%valeurs extrèmes
dif=reshape(bande1-bande2,numel(rv),1);
  
%Sur une zone sélectionnée
%  rv2=find((rv(indo))<1.4&(rv(indo))>0.7&dif(indo)~=0);
% load('D:\Raphael\THESE\NZ\codes\Shoreline/ZoneHistC3.mat');
 
load('ZoneGPP');
xt=x1;x1=y1;y1=xt;

rv2=find((rv(indo))<3&(rv(indo))>0.4&dif(indo)~=0);
%Genere un histogramme des intensitées de pixel
 
  [val0 nval]=hist([0.4;rv(indo(rv2));1.8],1000);
% [val0 nval]=hist(rv(rv2),1000);
%Filtrage de l'histogramme
val0=FiltreMean(val0,5);
val0=FiltreMean(val0,10);
val=FiltreMean(val0,90);
%%%AFFICHAGE%%%%%%%%%%%%%%%%%
if strcmp(Aff1,'yes')
%figure(9);clf; hold on;plot(nval,val,'g','LineWidth',2);set(gcf,'Color','w')
%figure(10);clf; hold on;plot(nval,val,'g','LineWidth',2);set(gcf,'Color','w')
end

 %Recherche des minima locaux
[vmin indmin]=lmin(val);
 %Recherche des maxima locaux et classement par ordre descendant des max
vmax=ones(40,1)*0;
vmin=ones(40,1)*0;
 
 [vx ix]=lmax(val);%maxima
 vmax(1:length(ix))=vx;
 [vmax indx]=sort(vmax,'descend');
 
 [vn in]=lmin(val);%Minima
 vmin(1:length(in))=vn;
 [vmin indx]=sort(vmax,'descend');
 seuil=nval(in);
 
    ith=0;
   % disp('   Début itération sur Histogramme')
  while   length(vmax(find(vmax>10)))<2|vmax(2)<0.1*vmax(1)|vmin(1)>0.9*vmax(2)
%  while (length(vmax(find(vmax>10)))==2&vmax(2)<0.5*vmax(1)&vmin(1)<0.8*vmax(2))~=1
vmax=ones(40,1)*0;
    ith=ith+1;
%      disp(['   ',num2str(ith)])
    val=FiltreMean(val0,10);
    val=FiltreMean(val,150-5*ith);
    val=FiltreMean(val,150-5*ith);
    [vx ix]=lmax(val);
    vmax(1:length(ix))=vx;
    [vmax indx]=sort(vmax,'descend');
  
    %Test de mauvaise image
     if(ith>=Conv1) ;break  ; end
    
if length(ix)>1 %Boucle pour localiser le minimum
    
      ix=ix(indx(1:2));
 
    ind=min(ix(1:2)):max(ix(1:2));
    [vn in]=lmin(val(ind));%Minima
    if numel(in)<1; in=round(length(ind)/2);vn=val(in);end
    vmin(1:length(in))=vn;
    [sv si]=min(vn);
    seuil=nval(ind(in(si)));
    
%     disp(['Minimum insuffisant : ',num2str(length(ix))])
end
%%%AFFICHAGE%%%%%%%%%%%%%%%%%
% if strcmp(Aff1,'yes'); 
%     figure(10);
%     hold on;
%     plot(nval,val,'k','LineWidth',2);
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
%%%AFFICHAGE%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if strcmp(Aff1,'yes'); %figure(10);
% if(ith<Conv1&ith>0)
%     hold on;plot(nval,val,'k','LineWidth',2);
%     hold on;plot(nval(ix(1:2)),vmax(1:2),'ks','LineWidth',2);
%     hold on;plot(seuil(1),vmin(1),'rs','LineWidth',2);
%     set(gca,'FontSize',14,'LineWidth',1);box on;grid on
% end
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
%Si trop sombre, pas de ligne d'eau
if max(nval(ix))<0.6; ith=0; end
 
% if diff(nval(ix))>0.13; seuil(1)=min(nval(ix)); end
clear dif rv2
 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Détection 2D de la ligne d'eau%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Initialisation de it
count=0;
NbPts=[];
DP=0;

%Critère de qualité de ligne d'eau : crit
crit='no';
 
if(ith<Conv1&ith>0)  
    
seuil0=min(nval(ix(1:2)))+cmin*abs(diff(nval(ix(1:2)))); 
seuil(1)=seuil0;
 
%disp('      Début Itération sur détection de la ligne d''eau')
 
%%%%AFFICHAGE%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%
if strcmp(Aff,'yes'); %figure(12);clf
%     subplot(2,2,1);box on;hold on;plot(it,new,'ks');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
%Resolution utilisée pr la définition de la zone
res=1; 
clear rv3
rv=double(a3(:,:,1))./double(a3(:,:,3));

[rv]=PreTraitShor_20080113_ter(rv);

% rv=smoothc(rv, 5, 5);
[m n]=size(rv);
rv(find(isnan(rv)==1))=0;

load('ZoneGPP');
xt=x1;x1=y1;y1=xt;
x1=round(x1/res);y1=round(y1/res);
Zon=ones(round(m/res),round(n/res)).*0;
%Matrice qui définit la zone avec des 1
for i=1:length(x1)
Zon(y1(i),x1(i))=1;
end

while strcmp(crit,'no')==1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%AFFICHAGE%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if strcmp(Aff,'yes'); figure(12);end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
count=count+1;

%Initialisation
  rv3=rv(1:res:m,1:res:n);
  rv3(find(rv3>=seuil(1)))=10;
  rv3(find(rv3<seuil(1)))=0;
    se1 = strel('octagon',15);
    se2 = strel('octagon',9);

%      rv3=imclose(rv3,se1);
%      rv3=imopen(rv3,se1);

     rv4=bwperim(rv3);
    [a b]=find(rv4==1);
    a=a.*res+1;
    b=b.*res+1;
clear x y
x=NaN;y=NaN;
 
% vo=0;
% for i=1:length(a)
% if(isempty(find(a(i)==y1&b(i)==x1, 1))==0); vo=vo+1; x(i)=a(i);y(i)=b(i); end;
% end
% 

%Methode 2 matricielle
vo=length(find(rv4.*Zon==1));



%%%AFFICHAGE%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(Aff,'yes');
    figure(12);clf;
subplot(2,2,1);hold on;plot(nval,val,'k','LineWidth',2);plot(seuil(1),mean(val)+2,'k.','LineWidth',2);grid on;box on;set(gca,'FontSize',14,'LineWidth',1)
h=legend(['c=',num2str(count),' s=',num2str(seuil(1))]);seuil(1)
subplot(2,2,2);image(a3(:,:,:))
hold on;plot(y.*res,x.*res,'g*','LineWidth',0.01);set(gca,'FontSize',14,'LineWidth',1)
zoom on
end

NbPts(count)=vo;
%Nombre de points sur y

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%AFFICHAGE%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(Aff,'yes'); figure(12);set(gcf,'Color','w')
subplot(2,2,3);hold on;box on;plot(count,st(count),'ks');grid on;box on;set(gca,'FontSize',14,'LineWidth',1)
hold on;plot(1:count,15*ones(count,1),'k')
subplot(2,2,4);hold on;box on;plot(count,DP(count),'ks');grid on;box on;set(gca,'FontSize',14,'LineWidth',1)
hold on;plot(1:count,0.020*n/res*ones(count,1),'k')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

seuil(1)=seuil(1)+0.01;
 
    
 if(count==25|count==50|count==75)
 %disp(num2str(count));
 end
 
if(seuil(1)>max(nval(ix(1:2)))-cmax*abs(diff(nval(ix(1:2))))); % On a fini le calcul des lignes d'eau
%disp('        Calcul de la ligne d''eau ');
crit='yes';
end
 
end %Boucle while

%Calcul de la variation du nombre de points
DP=[0 abs(diff(NbPts))];

clear DPAve
%Methode 2 % Sans forcer les bords de l'analyse de stabilité
in=min(find(DP~=0)):max(find(DP~=0));
sizein=0.25*length(in);
c=1;
ind=[];
while (length(ind)>nbmin | length(ind)<1 ) & c<sizein
    c=c+1;
DPAve(in)=FiltreMean(FiltreMean(FiltreMean(DP(in)./(NbPts(in)+1),2),3),c);
[val ind]=lmin(DPAve(in));
ind;
end

if length(ind)>0
    [frd tfh]=min(val);
indok=in(ind(tfh));
else
%recherche de la premiere valeur lors de la convergence
DPAve(in)=FiltreMean(FiltreMean(DP(in)./(NbPts(in)+1),2),round(sizein/3));
io=find(abs(DPAve(in)-min(DPAve(in)))<1/10*abs(max(DPAve(in))-min(DPAve(in))));
indok=in(io(1));
end


if(DPAve(indok)<0.1*n/res)
    
seuilok=seuil0+(indok+35)*0.01;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Recalcul de la ligne d'eau pour la valeur de seuil sélectionnée%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('ZoneGPP');
xt=x1;x1=y1;y1=xt;

count=count+1;
[m n]=size(rv);
%Resolution utilisée pr la définition de la zone
 
clear rv3
rv=double(a3(:,:,1))./double(a3(:,:,3));

[rv]=PreTraitShor_20080113_ter(rv);

rv(find(isnan(rv)==1))=0;
clear rv3

%Initialisation
  rv3=rv(1:res:m,1:res:n);
  rv3(find(rv3>=seuilok))=10;
  rv3(find(rv3<seuilok))=0;

    se4 = strel('octagon',60);
    se0 = strel('octagon',30);
    se1 = strel('octagon',15);
    se2 = strel('octagon',9);

%      
    %tic
    rv3=imopen(rv3,se4);
    %disp('Imopen')
    %toc;

    rv3=imopen(rv3,se0);
      %disp('Imopen')

    %toc;

    rv3=imopen(rv3,se0);
      %disp('Imopen')

    %toc;

    rv3=imopen(rv3,se1);
      %disp('Imopen')

    %toc;

    rv3=imopen(rv3,se2);
     % disp('Imopen')

    %toc;

    rv3=imopen(rv3,se2);
      %disp('Imopen')

    %toc;

    rv3=imclose(rv3,se0);
      %disp('Imclose')

    %toc;
     %disp('Finish')
    rv4=bwperim(rv3);
    [a b]=find(rv4==1);
    a=a.*res+1;
    b=b.*res+1;
clear x y
x=NaN;y=NaN;
 
vo=0;
for i=1:length(a)
if(isempty(find(a(i)==y1&b(i)==x1))==0); vo=vo+1; x(i)=a(i);y(i)=b(i); end;
end

ind=find(x>0&isnan(x)==0);
x=x(ind);
y=y(ind);

%toc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
%%%AFFICHAGE%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(Aff2,'yes')
vec=seuil0+(1:count-2).*0.002;
figure(14);clf;set(gcf,'Color','w')
subplot(2,1,1)

image(a3(:,:,:));axis image; axis ij
hold on
plot(y,x,'r*','LineWidth',0.01,'MarkerSize',2);set(gca,'FontSize',14,'LineWidth',1)
subplot(2,1,2);hold on
plot(vec(in(1):min([length(vec) max(in)])),DPAve(in(1):min([length(vec) max(in)]))/max(DPAve(in(1):min([length(vec) max(in)]))),'k')
plot(vec(indok),DPAve(indok)/max(DPAve),'ks');set(gca,'FontSize',14,'LineWidth',1)
set(gca,'XTick',vec(1:10:length(vec)))
xlim([seuil0 vec(count-1)])
ylabel ('\DeltaL/L');
xlabel ('RG');
box on
grid on
end
%%%Fin calcul de la ligne d'eau ok%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%disp('   Calcul de la ligne d''eau OK')

x2=(1:10);y2=(1:10).*0;

else
disp('   Pas de ligne d''eau - MAUVAISE IMAGE')
x2=(1:10);y2=(1:10).*0;
x=(1:10);y=(1:10).*0; 
end


else
disp('   Pas de ligne d''eau - MAUVAISE IMAGE')
x2=(1:10);y2=(1:10).*0;
x=(1:10);y=(1:10).*0; 

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
end

function [lmval,indd]=lmax_pw(xx, dx)
%   Find  piece-wise  local maxima in vector XX,where
%	LMVAL is the output vector with maxima values, INDD  is the 
%	corresponding indexes, DX is length of piece where maxima is searched, 
%   IMPORTANT:      FIRST and LAST point in vector are excluded
%   IMPORTANT:      XX must be single column vector
%   IMPORTANT:      Length of DX must be very carefully selected 
%	For example compare dx=30; and dx=1000;
%
%   dx=150; xx=[0:0.01:35]'; y=sin(xx .* cos(xx /4.5)) + cos(xx); 
%    plot(xx,y); grid; hold on;
%   %   Excluding first and last points
%   [b,a]=lmax_pw(y,dx); plot(xx(a),y(a),'r+')
%   % Way to include first and last points can be as:
%   y(1)=1.5; yy=[0; y; -1;];   % padd with smaller values
%   [b,a]=lmax_pw(yy,dx); a=a-1; plot(xx(a),y(a),'go')
%
%	see also LMIN, LMAX,  LMIN_PW, MATCH

% 	Sergei Koptenko, Applied Acoustic Technologies, Toronto, Canada
%   sergei.koptenko@sympatico.ca,  March/11/2003  

if nargin <2, 
	disp('Not enough arguments'); return
end

len_x = length(xx);
xx = [xx; xx(len_x); xx(len_x)]; 
nn=floor(len_x/dx);
ncount=1; lmval=[]; indd=[];
	for ii=1:nn,
        [lm,ind] = max(xx(ncount: ii*dx+2)) ;
        ind=ind+(ii-1)*dx;
                 if (ind ~=ncount) & (ind~=ii*dx+2),    
                    lmval=[lmval, lm]; indd=[indd, ind]; 
                end      
        ncount=ncount +dx;
	end
[lm,ind] = max(xx(ii*dx:len_x));
        if (ind ~=len_x) & (ind~=ii*dx),    
            lmval=[lmval, lm]; indd=[indd, (ind+ii*dx-1)]; 
        end
    
       if indd(end)==len_x,  
           indd=  indd(1:end-1); 
           lmval=lmval(1:end-1);    
       end
return
end

function [lmval,indd]=lmin_pw(xx, dx)
%   Find  piece-wise local minima in vector XX,where
%	LMVAL is the output vector with minima values, INDD  is the 
%	corresponding indexes, DX is scalar length of piece where minima is searched, 
%   IMPORTANT:      FIRST and LAST point in vector are excluded
%   IMPORTANT:      XX must be single column vector
%   IMPORTANT:      Length of DX must be very carefully selected 
%	For example compare dx=10; and dx=1000;
%
%   dx=150; xx=[0:0.01:35]'; y=sin(xx .* cos(xx /4.5)) + cos(xx); 
%   y(length(y))=-2; plot(xx,y); grid; hold on;
%   %   Excluding first and last points
%   [b,a]=lmin_pw(y,dx); plot(xx(a),y(a),'r+')
%   % Way to include first and last points can be as:
%   yy=[1.5; y; 0];         % padd with values higher than end values
%   [b,a]=lmin_pw(yy,dx); a=a-1; plot(xx(a),y(a),'go')
%
%	see also LMIN,LMAX, LMAX_PW, MATCH

% 	Sergei Koptenko, Applied Acoustic Technologies, Toronto, Canada
%   sergei.koptenko@sympatico.ca,  March/11/2003  

if nargin <2, 
	disp('Not enough arguments'); return
end

len_x = length(xx);
xx = [xx; xx(end)]; 
nn=floor(len_x/dx);
ncount=1; lmval=[]; indd=[];
for ii=1:nn,
    [lm,ind] = min(xx(ncount: ii*dx+1)) ;
        ind=ind+(ii-1)*dx;
         if (ind ~=ncount) & (ind~=ii*dx+1),    
         lmval=[lmval, lm]; indd=[indd, ind]; 
end      
ncount=ncount +dx;
end
[lm,ind] = min(xx(ii*dx:len_x));
    if (ind ~=len_x) & (ind~=ii*dx),    
    lmval=[lmval, lm]; indd=[indd, (ind+ii*dx-1)]; 
    end
    
     if indd(end)==len_x,
    indd=indd(1:end-1); lmval=lmval(1:end-1);
    end
return
end

function mO = smoothc(mI, Nr, Nc)

% SMOOTHC.M: Smooths matrix data, cosine taper.
% MO=SMOOTHC(MI,Nr,Nc) smooths the data in MI
% using a cosine taper over 2*N+1 successive points, Nr, Nc points on
% each side of the current point.
%
% Inputs: mI - original matrix
% Nr - number of points used to smooth rows
% Nc - number of points to smooth columns
% Outputs:mO - smoothed version of original matrix
%
%
if nargin<2, error('Not enough input arguments!'), end;

% Determine convolution kernel k
Nr=Nr+1;
Nc=Nc+1;
kr=2*Nr+1;
kc=2*Nc+1;
midr=Nr+1;
midc=Nc+1;
maxD=(Nr.^2+Nc.^2).^0.5;
for irow=1:kr;
for icol=1:kc;
D=((midr-irow).^2+(midc-icol).^2).^(0.5);
k(irow,icol)=cos(D*pi/2./maxD);
end;
end;

k = k./sum(k(:));

% Perform convolution
mO=conv2(mI,k,'valid');
end


function [rvout]=PreTraitShor_20080113_ter(rvin)


%addpath '/ProgSecondaires'
[m n]=size(rvin);
res=2;
%Generation de la correction à appliquer + rouge sur la plage et + foncé ds
%l'eau
%Correction linéaire haut/bas
b=-0.02;
a=(2/m)*(-b);
%Angle de la ligne d'eau (en pix)
%Axe diagonal : ang=m/n /axe perp : ang=0/n
ang=0/n;
[X1,Y1]=meshgrid(b+(1:m)*a,b+(1:n)*a*ang);
%Correction dépendante de la distance au centre de l'image
% x=0.001*(abs(m/2-(1:m)).^2)/max(abs(m/2-(1:m)).^2);
% y=0.1*(abs(n/2-(1:n)).^2)/max(abs(n/2-(1:n)).^2);

x=0.001*(abs(m/2-(1:m)))/max(abs(m/2-(1:m)));
y=0.08*(abs(n/2-(1:n)))/max(abs(n/2-(1:n)));

[X2,Y2]=meshgrid(x,y);

%Smoothing
%Facteur horizontal
facth=2;
%Facteur vertical
factv=2;

rvout=rvin;

rvout(factv+1:m-factv-2,facth+1:n-facth-2)=smoothc(rvin,factv,facth);

%   figure(512);surfc(rvout(50:res:m-50,80:res:n-80));shading flat;set(gcf,'Color','w');box on
% % 
% % 
% %   figure(513);surfc(rvout(50:res:m-50,50:res:n-50)+X1(50:res:n-50,50:res:m-50)'+Y1(50:res:n-50,50:res:m-50)');shading flat;set(gcf,'Color','w');box on
% % 
%  figure(514);surfc(rvout(50:res:m-50,50:res:n-50)+X1(50:res:n-50,50:res:m-50)'+Y1(50:res:n-50,50:res:m-50)'+Y2(50:res:n-50,50:res:m-50)');shading flat;set(gcf,'Color','w');box on

rvout=rvout+X1'+Y1'+Y2';
end