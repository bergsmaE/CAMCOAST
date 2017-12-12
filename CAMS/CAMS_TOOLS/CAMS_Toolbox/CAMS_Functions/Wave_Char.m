function [hs,htiers,hrms,trms,Hmax,h,Tp,sl,tcrete,t]=Wave_Char(d,dt,filt,meth)
%Author: Rafael Almar (rafael.almar@ird.fr)
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
%z=temps des vagues
%Outputs
%h : hauteurs de vague individuelles (ex: pour faire un histogramme)
%htiers;
%hrms;
%trms;
%Hmax;
%h; hauteurs individuelles des vagues (vecteur)
%Tp; periode peak
%sl: slope (m/s) (vecteur)
%tcrete= temps des cretes individuelles  (vecteur)
%t:periode individuelle (vecteur)

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
%  x=[1:nt];
%  a=reshape(x,nt,1);
%  
%  % generation des coefficients du polynome de degré un
%   p=polyfit(a,k,1);
%  
%  %Suppression de la partie trendée de la courbe
%  
%  k2=k-p(1).*a;
%  
%  % Suppression de la partie moyenne
%  km=k2;
%  k2=km-mean(km);
 
k2=detrend(d);

 clear k4
if(filt==1)
freq = (1./dt)*(1:nt)/nt;
 %frequences de coupure : %%%%%%%%%%%%%%%%%%%%%
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %fréquence de coupure vers les basses fréqences
 % fcb= 1/Tlimit
  fcb=0.05;
 % fréquence de coupure vers les hautes fréquences
 % fch=1/Tlimit
  fch=0.17;
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
t=[];
kh=0;
%hauteurs: h
% h(1:k-2)=0;
%pentes:sl
sl=[];
for i=1:k-2
    
 if(hmin(i)~=0 & hmax(i+1)~=0)
        if(hmin(i+1)==0 & hmax(i)==0)        
        kh=kh+1;
     h(kh)=hmax(i+1)-hmin(i);
     sl(kh)=hmax(i+1)-hmin(i)./((z(i+1)-z(i)));
     try
     tcrete(kh)=z(i+1)+0.1*abs((z(i)-z(i+1)));
     t(kh)=2*(z(i+1)-z(i))*dt;     
     end
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
Hmax=nanmax(h);

%Hauteur 1/3
htiers=nanmean(h(ord(1:round(length(ord)/3))));% on prend le premier tiers pour la hauteur rms
%Hauteur rms
hrms=nanmean(h(1:kh));

%Hauteur significative
hs=4.*std(k4);
% disp(['Hs=',num2str(hs),'; Hrms=',num2str(hrms),'; H1/3=',num2str(htiers),'; Hmax=',num2str(Hmax)])
%Calcul de la periode de la houle pour chaque vague

% clear t
% t(1:k-1)=0;
% for i=1:k-1
%     if(z(i)~=0 | z(i+1)~=0 ) 
%         t(i)=2*(z(i+1)-z(i))*dt;
%     end
% end

%Calcul de la periode rms de la houle et affichage

if(zc==1)
trms=nanmean(t);
elseif(zc==2)
trms=nanmean(t)/2;
elseif(zc==3)
trms=nanmean(t)/2;
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
Pv = pwelch(detrend(d));
[valmx indmx]=max(Pv);
Tp=1./(indmx./(dt*length(Pv)));

%Affichage de la periode
% disp(['Periode = ',num2str(trms)])

% try
% %Positin des cretes des vagues
% iii=find(tt>0);
% tc=tt(iii)+0.5*diff([tt(iii) max(tt(iii))]);
% id=round(tc);
% tcrete=id(find(k4(id)>0));
% catch
%     tcrete=[];
% end

h=h(1:kh);
% end%zc

% figure;
% subplot(2,1,1)
% plot((1:nt)./(dt*nt),f1.*conj(f1))
% legend('Spectre de houle filtré')
% subplot(2,1,2)
% plot((1:nt)./(dt*nt),fft(k2).*conj(fft(k2)))
% legend('Spectre de houle initial')