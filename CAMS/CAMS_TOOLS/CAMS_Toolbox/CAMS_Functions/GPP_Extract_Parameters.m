%Fonction permettant l'extraction et création des vecteurs Hauteurs significatives, Hauteurs max, Periodes Pic et
%Periodes significatives.
%Cette fonction est utilisée afin de comparer les valeurs obtenues à GPP avec celles émises par la bouée. 
%[Date_gpp,Hs_gpp,Hm_gpp,Tp_gpp] = ExtractGPP(Directory,Date_deb,Date_f);
%eg. ExtractGPP('C:\Users\HP\Desktop\GrandPopo_Data','20160203','20160303')

function [Date_gpp,Hs_gpp,Hm_gpp,Tp_gpp] = GPP_Extract_Parameters(Directory,date_deb,date_f);

cd(Directory);
l = dir([Directory '\GPP_Parameter*']);

Date_debut =  datenum(date_deb,'yyyymmdd');
Date_fin = datenum(date_f,'yyyymmdd');

liste = [];

for i=1 : length(l);
    file_name = l(i).name;
    
    d = datenum(file_name(end-11:end),'yyyymmdd');
    
    if d >= Date_debut & d <= Date_fin ;
        
        liste = [liste; file_name];
        
    end
    
end

Date_gpp = []; Hs_gpp = []; Hm_gpp = []; Tp_gpp = []; 

ll=size(liste,1);

for j = 1 : ll
    
    load(liste(j,:));
    
    Date_gpp = [Date_gpp dateT];
    Hs_gpp = [Hs_gpp hsT];
    Hm_gpp = [Hm_gpp hmT]; 
    Tp_gpp = [Tp_gpp TpT]; 
   
   
end

%Data_gpp = [Date_gpp; Hs_gpp; Hm_gpp; Tp_gpp]';
%Data_gpp_text = {'Date_gpp' 'Hs_gpp' 'Hm_gpp' 'Tp_gpp'};

%cd('C:\Users\HP\Desktop\GrandPopo_Data\03-GrandPopo_Postprocess_Data\03.1-Wave_Parameters_extracted');
%filename = strcat('GPP_Hs-Hm-Tp_',datestr(Date_debut,'yyyymmdd'),'-',datestr(Date_fin,'yyyymmdd'));

%save(filename,'Date_gpp','Hs_gpp','Hm_gpp','Tp_gpp','Data_gpp','Data_gpp_text');

end



    