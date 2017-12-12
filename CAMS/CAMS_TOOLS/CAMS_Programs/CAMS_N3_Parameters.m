function CAMS_N3_Parameters(dirN2_p,dirN2_d,dirN2_s,dirN3)

%Extraction de Hm, Tm, et dateT:
cd(dirN2_p);
liste=dir('GPP_Param*.mat');
Hsn=[];Tmn=[];date=[];
for i=1:length(liste)   
    load(liste(i).name)
    Hsn=[Hsn hsT];
    Tmn=[Tmn TmT];
    date=[date dateT]; 
end 
[date_p ordX]=sort(date);
Hs=Hsn(ordX);
Tm=Tmn(ordX); 

%Extraction de la direction et des dates associées:
cd(dirN2_d);
liste=dir('*Dir*.mat');
Dir=[]; date_d=[];

for i=1:length(liste)   
load(liste(i).name)
Dir =[Dir angi];
date_d=[date_d  dati]; 
end 
    [date_d ordX]=sort(date_d);
    Dir = Dir(ordX); 

%Extraction du trait de cote:
cd(dirN2_s);
liste=dir('GPP_ligne-eau_*.mat');
Yn=[];date=[];
for i=1:length(liste)
    load(liste(i).name)
    Yn=[Yn Y];
    date=[date dateavg]; 
    [MT ordX]=sort(date);
    MY=Yn(:,ordX);
end 

%Elimination des valeurs abérantes 
indHs = find(Hs > (2.*nanstd(Hs)+nanmean(Hs)));
Hs(indHs) = NaN;

%Pour le trait de cote 
ind=find(nanmean(MY)>=50&nanmean(MY)<=80); % Indices de lignes d'eau exploitables
MYe=nanmean(MY(:,ind)); % vecteur instantané (toutes les 15 min) de ligne d'eau exploitable moyenné long-shore
date_s=MT(ind);

%LES SERIES TEMPOELLES Hs, Tp, Dir et MY
%Création des séries temporelles avec des NaN aux instants où il n'y a pas de données

% Série temporelle avec des intervalles de 15 min = 1/96 jour Julien
temps = min([date_p(1),date_d(1),date_s(1)]):1/96:max([date_p(end),date_d(end),date_s(end)]);

% Extraction des valeurs de "Hs_gpp" et "Tp_gpp" (vidéo) correspondantes à "temps". 
Hs_temps = nan(1,length(temps));
Tm_temps = nan(1,length(temps));
Dir_temps = nan(1,length(temps));
MYe_temps = nan(1,length(temps));
%On veut obtenir une série temporelle avec des "NaN" aux instants du
%vecteur temps où on ne dispose pas de données.
for i=1:size(temps,2)  
%     On cherche les dates pour lesquelles on a une valeur de Hs. Le
%     +0.00007 (= 06 min) permet d'avoir une marge d'erreur sur les dates.         
    ind_t = find(date_p >= temps(i)-0.00007 & date_p <= temps(i)+0.00007);
    ind_d = find(date_d >= temps(i)-0.00007 & date_d <= temps(i)+0.00007);
    ind_s = find(date_s>= temps(i)-0.00007 & date_s <= temps(i)+0.00007);

    if length(ind_t)==1
        Hs_temps(i) = Hs(ind_t);
        Tm_temps(i) = Tm(ind_t);       
%     else
%         Hs_temps(i) = NaN;
%         Tm_temps(i) = NaN;
    end
    
    if length(ind_d)==1
        Dir_temps(i) = Dir(ind_d);
%     else
%         Dir_temps(i) = NaN;
    end
    
    if length(ind_s)==1
        MYe_temps(i) = MYe(ind_s);      
%     else
%          MYe_temps(:,i) = NaN;
    end
end

Hs = Hs_temps;
Tm = Tm_temps;
Dir = Dir_temps;
MY = MYe_temps;
    
%CALIBRATION
%Calcul des coefficients de calibration
[C1_Hs,C2_Hs,C1_Tm,C2_Tm,C1_Dir,C2_Dir]  = GPP_calibr_ADCP_VID(Hs,Tm,temps,Dir);

cd(dirN3);
save('GPP_N3_Coeff_Calibration','C1_Hs','C2_Hs','C1_Tm','C2_Tm','C1_Dir','C2_Dir');

%Hs Calibrée
Hs =(Hs.*C1_Hs) + C2_Hs; 

%Tm Calibrée 
Tm =(Tm.*C1_Tm) + C2_Tm; 

%Dir Calibrée
Dir = (Dir.*C1_Dir) + C2_Dir;
 
%MOYENNES JOURNALIERES ET ECARTS TYPES DES PARAMETRES DE VAGUES
%Moyenne journalière et écarts types de Hs et Tm:

% Vecteurs temps, liste des jours.
liste_j = datenum((datestr(temps,'yyyymmdd')),'yyyymmdd');%liste des jours du vecteur temps, un jour  apparait plusieurs fois. 
temps_j = unique(liste_j); %liste des jours du vecteur temps, chaque jour apparait une fois

Hs_j = []; 
Std_Hs_j = [];
Tm_j = [];
Std_Tm_j = [];
Dir_j=[];
Std_dir_j=[];
MY_j = []; 
Std_MY_j = [];    
for j = 1 : length(temps_j)
        ind_mj = find(liste_j == temps_j(j));
        
        Hs_j = [Hs_j nanmean(Hs(ind_mj))];
        Std_Hs_j = [Std_Hs_j nanstd(Hs(ind_mj))];
        
        Tm_j = [Tm_j nanmean(Tm(ind_mj))];
        Std_Tm_j = [Std_Tm_j nanstd(Tm(ind_mj))];
        
        Dir_j=[Dir_j nanmean(Dir(ind_mj))];
        Std_dir_j=[Std_dir_j nanstd(Dir(ind_mj))];
        
        MY_j = [MY_j nanmean(MY(ind_mj))];
        Std_MY_j = [Std_MY_j nanstd(MY(ind_mj))];
                       
end

%MOYENNES MENSUELLES ET ECARTS TYPES DES PARAMETRES DE VAGUES    

% Vecteurs temps, liste des jours.
liste_m = datenum((datestr(temps,'yyyymm')),'yyyymm');%liste des mois du vecteur temps, un mois  apparait plusieurs fois. 
temps_m = unique(liste_m); %liste des mois du vecteur temps, chaque mois apparait une fois

Hs_m = []; 
Std_Hs_m = [];
Tm_m = [];
Std_Tm_m = [];
Dir_m=[];
Std_dir_m=[];
MY_m = []; 
Std_MY_m = [];

for j = 1 : length(temps_m)
        ind_mj = find(liste_m == temps_m(j));
        
        Hs_m = [Hs_m nanmean(Hs(ind_mj))];
        Std_Hs_m = [Std_Hs_m nanstd(Hs(ind_mj))];
        
        Tm_m = [Tm_m nanmean(Tm(ind_mj))];
        Std_Tm_m = [Std_Tm_m nanstd(Tm(ind_mj))];
        
        Dir_m=[Dir_m nanmean(Dir(ind_mj))];
        Std_dir_m=[Std_dir_m nanstd(Dir(ind_mj))];
        
        MY_m = [MY_m nanmean(MY(ind_mj))];
        Std_MY_m = [Std_MY_m nanstd(MY(ind_mj))];        
end


%FIGURES MOYENNES JOURNALIERES 
fig_j = figure;;
    subplot(3,1,1);
     plot(temps_j,Hs_j,'LineWidth',2); 
     title ('Evolution journalière de la hauteur significative');
     %xlim([min(date2) max(date2)]);
     xlabel('time (YY/mm)');
     ylabel ('<Hs_j>_Y (m)'); set(gca,'XTick',temps_j(1):60:temps_j(end));
     set(gca,'XTick',temps_j(1):60:temps_j(end));
     set(gca,'XTickLabel',datestr(temps_j(1):60:temps_j(end),'YY/mm'));
     legend('hauteur significative'); 

    subplot(3,1,2);
     plot(temps_j,Tm_j,'LineWidth',2); 

     %xlim([min(date2) max(date2)]);
     title ('evolution journalière de la periode moyenne');
     xlabel('time (YY/mm)');
     ylabel ('<Tmjr>_Y (m)');
     set(gca,'XTick',temps_j(1):60:temps_j(end));
     set(gca,'XTickLabel',datestr(temps_j(1):60:temps_j(end),'YY/mm'));
     legend('periode moyenne'); 

     subplot(3,1,3)

    plot(temps_j,Dir_j,'LineWidth',2);
    legend('Direction journalière','Location','NorthEast');
    ylabel('Direction(°)')
    set(gca,'XTick',temps_j(1:31:end))
    set(gca,'XTickLabel',datestr(temps_j(1:31:end),'yy/mm'));
    xlabel('time (YY/mm)')
    set(gcf,'Color','w')
    title('Evolution de la direction des vagues de GPP par vidéo')
    %patch([temps_j(1) temps_j(end)], [0 0], [1 0 0],'LineWidth',2);
    a=Dir_j(1);
    b=Dir_j(end);
    text(temps_j(end-100),-3,'Incidence cross-shore');
    ylim([min(Dir_j)-5 max(Dir_j)+5])
%xlim([(temps_j(1)-1) (temps_j(end)+1)])


fig_m = figure;
    subplot(3,1,1);
        plot(temps_m,Hs_m,'LineWidth',2); hold on; plot(temps_m,Hs_m-Std_Hs_m,'r','LineWidth',2); plot(temps_m,Hs_m+Std_Hs_m,'r','LineWidth',2)
        title ('Evolution mensuelle de la hauteur significative');
        %xlim([min(date2) max(date2)]);
        xlabel('time (YY/mm)');
        ylabel ('<Hs_j>_Y (m)'); set(gca,'XTick',temps_m(1):60:temps_m(end));
        set(gca,'XTick',temps_m(1):60:temps_m(end));
        set(gca,'XTickLabel',datestr(temps_m(1):60:temps_m(end),'YY/mm'));
        legend('hauteur significative','ecart-type'); hold off

        subplot(3,1,2);
        plot(temps_m,Tm_m,'LineWidth',2); 
        hold on; 
        plot(temps_m,Tm_m-Std_Tm_m,'r','LineWidth',2);
        plot(temps_m,Tm_m+Std_Tm_m,'r','LineWidth',2);
        %xlim([min(date2) max(date2)]);
        title ('Evolution journalière de la periode moyenne');
        xlabel('time (YY/mm)');
        ylabel ('<Tmjr>_Y (m)');
        set(gca,'XTick',temps_m(1):60:temps_m(end));
        set(gca,'XTickLabel',datestr(temps_m(1):60:temps_m(end),'YY/mm'));
        legend('periode moyenne','ecart-type'); hold off

    subplot(3,1,3)
        plot(temps_m,Dir_m,'xb','LineWidth',2);
        legend('Direction journalière','Location','NorthEast');
        ylabel('Direction(°)')
        set(gca,'XTick',temps_j(1:31:end))
        set(gca,'XTickLabel',datestr(temps_j(1:31:end),'yy/mm'));
        xlabel('time (YY/mm)')
        set(gcf,'Color','w')
        title('Evolution de la direction des vagues de GPP par vidéo')
        %patch([temps_j(1) temps_j(end)], [0 0], [1 0 0],'LineWidth',2);
        a=Dir_j(1);
        b=Dir_j(end);
        text(temps_j(end-100),-3,'Incidence cross-shore');
        ylim([min(Dir_j)-5 max(Dir_j)+5])
        %xlim([(temps_j(1)-1) (temps_j(end)+1)])
        

fig_s = figure;
    subplot(2,1,1);
        plot(temps_j,MY_j-MY_j(1),'LineWidth',2);
        title ('Evolution journalière de la ligne d''eau');
        %xlim([min(date2) max(date2)]);
        xlabel('time (YY/mm)');
        ylabel ('<Position trait de côte (m)'); set(gca,'XTick',temps_j(1):60:temps_j(end));
        set(gca,'XTick',temps_j(1):60:temps_j(end));
        set(gca,'XTickLabel',datestr(temps_j(1):60:temps_j(end),'YY/mm'));
        legend('Ligne d''eau');

    subplot(2,1,2);
        plot(temps_m,MY_m,'LineWidth',2); hold on; plot(temps_m,MY_m-Std_MY_m,'--','color','red'); plot(temps_m,MY_m+Std_MY_m,'--','color','red')
        title ('Evolution mensuelle de la ligne d''eau');
        %xlim([min(date2) max(date2)]);
        xlabel('time (YY/mm)');
        ylabel ('Position trait de côte (m)'); set(gca,'XTick',temps_j(1):60:temps_j(end));
        set(gca,'XTick',temps_j(1):60:temps_j(end));
        set(gca,'XTickLabel',datestr(temps_j(1):60:temps_j(end),'YY/mm'));
        legend('Ligne d''eau','ecart-type'); hold off

        
cd(dirN3);

saveas(fig_j,'GPP_N3_Evolution_parametres_jours.jpg');
saveas(fig_m,'GPP_N3_Evolution_parametres_mois.jpg');
saveas(fig_s,'GPP_N3_Evolution_Shoreline_jours.jpg');

filename_jour = ['GPP_N3_JOURS_' datestr(temps(1),'yyyymmdd') '-' datestr(temps(end),'yyyymmdd') '.txt'];
fid1 = fopen(filename_jour,'wt');
fprintf(fid1,'%s\t',' Date  ');
fprintf(fid1,'%s\t',' Hs(m)  ');
fprintf(fid1,'%s\t',' Tp(s)  ');
fprintf(fid1,'%s\t','Dir(°)');
fprintf(fid1,'%s\n','Trait de côte(m)');

for i=1:length(Hs_j)
    fprintf(fid1,'%s\t',datestr(temps_j(i)));
    fprintf(fid1,'%8.4f\t',Hs_j(i));   
    fprintf(fid1,'%8.4f\t',Tm_j(i));
    fprintf(fid1,'%1.0f\t',Dir_j(i));
    fprintf(fid1,'%8.4f\n',MY_j(i));
end

filename_mois = ['GPP_N3_MOIS_' datestr(temps(1),'yyyymmdd') '-' datestr(temps(end),'yyyymmdd') '.txt'];
fid2 = fopen(filename_mois,'wt');
fprintf(fid2,'%s\t',' Date  ');
fprintf(fid2,'%s\t',' Hs(m)  ');
fprintf(fid2,'%s\t',' Std(Hs)(m)  ');
fprintf(fid2,'%s\t',' Tp(s)  ');
fprintf(fid2,'%s\t',' Std(Tp)(s) ');
fprintf(fid2,'%s\t',' Dir(°)  ');
fprintf(fid2,'%s\t',' Std(Dir)(°)  ');
fprintf(fid2,'%s\t','Trait de côte(m)');
fprintf(fid2,'%s\n',' Std(trait côte)(m)');

for i=1:length(Hs_m)
    fprintf(fid2,'%s\t',datestr(temps_m(i)));
    fprintf(fid2,'%8.4f\t',Hs_m(i));
    fprintf(fid2,'%8.4f\t',Std_Hs_m(i));
    fprintf(fid2,'%8.4f\t',Tm_m(i));
    fprintf(fid2,'%8.4f\t',Std_Tm_m(i));
    fprintf(fid2,'%8.4f\t',Dir_m(i));
    fprintf(fid2,'%8.4f\t',Std_dir_m(i));
    fprintf(fid2,'%8.4f\t',MY_m(i));
    fprintf(fid2,'%8.4f\n',Std_MY_m(i));
    
end

cd(dirN3);
filename = ['GPP_N3_Parameters_' datestr(temps(1),'yyyymmdd') '-' datestr(temps(end),'yyyymmdd')];
save(filename,'temps_j','Hs_j','Tm_j','Dir_j','MY_j','Std_Hs_j','Std_MY_j','Std_Tm_j','Std_dir_j',...
    'temps_m','Hs_m','Tm_m','Dir_m','MY_m','Std_Hs_m','Std_Tm_m','Std_dir_m','Std_MY_m');
end



