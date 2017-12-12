function [Hs_coef1,Hs_coef2,Hs_erms,Tm_coef1,Tm_coef2,Tm_erms,dir_erms] = GPP_calibr_ADCP_VID_2(path)
cd(path)
liste=dir('GPP_param*.mat');
 Hsn=[];Tmn=[];date=[];
 for i=1:length(liste);
 load(liste(i).name)
 Hsn=[Hsn hsT];
 Tmn=[Tmn TmT];
 date=[date dateT]; 
 [date1 ordX]=sort(date);
  Hs=Hsn(ordX);
    Tm=Tmn(ordX);
 end 

clear date
load('D:\M2_2015\matlab\data\stage\ADCP_GPP2014_wave_tide.mat'); 

%%% Mes modifs
mondatef=str2num(datestr(date,'yyyymmdd'));
days=unique(mondatef);

mondatefv=str2num(datestr(date1,'yyyymmdd'));
daysv=unique(mondatefv);

Hs1=[];
HsADCP=[];
Tm1=[];
TpADCP=[];
date2=[];

for i=1:length(days)
    ind_v=find(mondatefv==days(i)); % indices de "date2" (vidéo) correspondantes au jour i
%     ind_v1=find(max(date1(ind_v))); % Indice de la date max de "date2" (vidéo) du jour i
    date_maxvid=max(date1(ind_v)); % Date max de "date2" (vidéo) du jour i
%     date_maxvid=date1(ind_v(ind_v1)); 
%     ind_v2=find(min(date1(ind_v))); % Indice de la date min de "date2" (vidéo) du jour i
    date_minvid=min(date1(ind_v));
%     date_minvid=date1(ind_v(ind_v2)); % Date min de "date2" (vidéo) du jour i
    ind_d=find(mondatef==days(i)); % indices de "date" (ADCP) correspondantes au jour i
    date_test=date(ind_d); % Dates ADCP du jour i
    ind_d1=find(date_test>=date_minvid&date_test<=date_maxvid); % Indices
    % ADCP du jour i correspondants à la période de fonctionnement de la
    % caméra (la camera ne fonctionne pas la nuit)
    
    date_int=date(ind_d(ind_d1)); % Dates ADCP du jour i pendant la période
    % de fonctionnement de la caméra
    date2=[date2 date_int]; % Vecteur contenant toutes les Dates ADCP du
    % jour i pendant la période de fonctionnement de la caméra
    
    Hs2v=interp1(date1(ind_v),Hs(ind_v),date_int); % Hs vidéo interpolés
    % sur les temps ADCP durant le fonctionnement de la camera
    HsADCP=[HsADCP;Hs_swell(ind_d(ind_d1))]; % Hs ADCP correspondants aux
    % temps ADCP durant le fonctionnement de la camera
    Hs1=[Hs1 Hs2v];
    
    Tm2v=interp1(date1(ind_v),Tm(ind_v),date_int); % Tm vidéo interpolés
    % sur les temps ADCP durant le fonctionnement de la camera
    TpADCP=[TpADCP;Tp(ind_d(ind_d1))]; % Tm ADCP correspondants aux
    % temps ADCP durant le fonctionnement de la camera
    Tm1=[Tm1 Tm2v];
end

 Hs_swell1=Hs_swell;
 Tp1=Tp;
clear Tp;
 s1=nanstd(Hs_swell1);
 s2=nanstd(Hs1);
  s3=nanstd(Tp1);
 s4=nanstd(Tm1);
 Hs2=(s1/s2).*Hs1;
 Hs_coef1=round((s1/s2)*10)/10;
 Tm_coef1=round((s4/s3)*10)/10;
  
  %%%% correction de la hauteur significative (Hs-vidéo)%%%%%%%%%
H1=nanmean(Hs2);
H2=nanmean(Hs_swell1);
Hs_coef2=round((H2-H1)*10)/10;
%Hs3=Hs2+H3;
   Hs3=Hs_coef1.*Hs1+Hs_coef2;

T1=nanmean(Tm1);
T2=nanmean(Tp1);
Tm_coef2=round((T2-T1)*10)/10;
%Tm1=Tm1+T3;
  Tm1=Tm_coef1.*Tm1+Tm_coef2;

%%%%%%%%%%%%%%% pour la direction des vagues%%%%%%%%%%%%%%%%%%%%%
Dir_cor=abs(169-Dir);
%load('D:\M2_2015\matlab\data\stage\Wave_Angle_GPP.mat');
load('Wave_Angle_GPP.mat')
angi=angi-90;
ind_st=find(dati>=min(date)&dati<floor(min(date))+1);
for i=1:length(ind_st)
    if dati(ind_st(i))==min(dati(ind_st))
        ind_star=ind_st(i);
    end
end
ind_en=find(dati>=floor(max(date))&dati<max(date));
if isempty(ind_en)
    ind_en=find(dati>=floor(max(date))-1&dati<floor(max(date)));
    for i=1:length(ind_en)
        if dati(ind_en(i))==max(dati(ind_en))
            ind_end=ind_en(i);
        end
    end
else
    for i=1:length(ind_en)
        if dati(ind_en(i))==max(dati(ind_en))
            ind_end=ind_en(i);
        end
    end
end
ind_interp=ind_star;
ind_dir=1;
dir_interp=zeros(1,length(dati(ind_star:ind_end)));

while (ind_interp<=ind_end)
    dir_interp(ind_dir)=interp1(date,Dir_cor,dati(ind_interp));
    ind_interp=ind_interp+1;
    ind_dir=ind_dir+1;
end

angi_vid=angi(ind_star:ind_end);
dif=angi_vid-dir_interp;
angi_vid=angi(ind_star:ind_end);
dif=angi_vid-dir_interp;

% Erreur moyenne
err_moy=nanmean(dif);

formatSpec = '    Erreur moyenne : %d degrés\n';
sprintf(formatSpec,err_moy)

%%%%%%%%%%%%%%%%%%% les figures %%%%%%%%%%%%%%%%%%%%%%
figure;
subplot(3,1,1);
plot(date2,Hs3,'sk','LineWidth',2,'MarkerEdgeColor','k','MarkerFaceColor','k'); hold on; 
plot(date,Hs_swell1,'sr','LineWidth',2,'MarkerEdgeColor','r','MarkerFaceColor','r');
% ylim([0.2 1.3]);
set(gca,'XTick',date(1):1:date(end));
 set(gca,'XTickLabel',datestr(date(1):1:date(end),'dd mmm'),'fontsize',12);
 legend('video','ADCP');
  x=datenum(2014,03,16);
  text(datenum(2014,03,11,07,45,00),0.7,'(a)','fontsize',13);
  y=datenum(2014,03,12,12,00,00);
  Hs_erms=round(rmse(Hs3,HsADCP')*100)/100;
  %title({['Validation Hauteur Hs video par ADCP (Campagne Grand Popo Mars 2014'];[' RMSE ',num2str(r),'m']},'FontSize',12);%['Hsvidcor=(Hsvid.*',num2str(ss1),') + ',num2str(H3)]})
  title({['validation Hs-vidéo'];[' RMSE ',num2str(Hs_erms),'m']})%'FontSize',12)
%   text(y,1.2,'HsADCP=(Hsvid.*0.17) + 0.48','fontsize',12);
 %xlabel('mois/jour/heure-2014','FontSize',12);
 ylabel('Hauteur (m)','FontSize',12);
 hold off;
subplot(3,1,2);
plot(date2,Tm1,'sk','LineWidth',2,'MarkerEdgeColor','k','MarkerFaceColor','k'); hold on ; plot(date,Tp1,'sr','LineWidth',2,'LineWidth',2,'MarkerEdgeColor','r','MarkerFaceColor','r'); hold off
%ylim([10 15]);
set(gca,'XTick',date(1):1:date(end));
 set(gca,'XTickLabel',datestr(date(1):1:date(end),'dd mmm'),'fontsize',12);
 legend('video','ADCP');
 ylabel('Période(s)','FontSize',12);
  x=datenum(2014,03,15,21,00,00);
  y=datenum(2014,03,12);
  text(datenum(2014,03,11,07,00,00),10,'(b)','fontsize',13);
  
  Tm_erms=round(rmse(Tm1,TpADCP')*100)/100;
%title({['Validation Période-video par ADCP (Campagne Grand Popo Mars 2014'];[' RMSE ',num2str(r1),'s']},'FontSize',12)%;['Tmvidcor=(Tmvid.*',num2str(ss2),') + ',num2str(T3)]})
title({['validation Tm-vidéo'];[' RMSE ',num2str(Tm_erms),'s']})
  
  subplot(3,1,3);
  %plot(date,Dir_cor,'*r'); hold on
  plot(dati(ind_star:ind_end),angi_vid,'sk','LineWidth',2,'MarkerEdgeColor','k','MarkerFaceColor','k'); hold on; plot(dati(ind_star:ind_end),dir_interp,'sr','LineWidth',2,'LineWidth',2,'MarkerEdgeColor','r','MarkerFaceColor','r'); hold off
set(gca,'XTick',date(1):1:date(end));
set(gca,'XTickLabel',datestr(date(1):1:date(end),'dd mmm'),'fontsize',12);
xlabel('mois/jour/Mars-2014','FontSize',12);
ylabel('Incidence des vagues (°)','FontSize',12) 
title('validation Direction de la houle','FontSize',12)
legend('video','ADCP')
patch([dati(ind_star) dati(ind_end)], [0 0], [1 0 0],'LineWidth',3);%%tracer droite
text(datenum(2014,03,17,16,0,0),0,'Incidence')
text(datenum(2014,03,11,10,45,00),55,'(c)','fontsize',13);
ylim([-10 60])
hgexport(gcf, 'figure1.jpg', hgexport('factorystyle'), 'Format', 'jpeg');%% pour enregistrer une figure
%EVALUATION du RMSE Dir_interp (ADCP) et angi (video)
dir_erms=rmse(angi_vid,dir_interp);

cor=corrcoef(angi_vid',dir_interp'); % Calcul de la correlation
% entre les directions de houle de la video et celles de l'ADCP

formatSpec = '    Root Mean Square Error : %d \n';
% sprintf(formatSpec,RMSE)

formatSpec = '    Correlation : %d \n';
sprintf(formatSpec,cor(1,2))

text(datenum(2014,03,12,0,0,0),50,'RMSE : 11.02,  Erreur moyenne : 3°','FontSize',12)
end