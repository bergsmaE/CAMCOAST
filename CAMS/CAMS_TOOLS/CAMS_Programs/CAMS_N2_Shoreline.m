%PostProcess_GPP.m
%%% Calcule le trait de côte.
%%% Chemin conduisant aux images à traiter %%%%%%%%%%%%

function CAMS_N2_Shoreline(dirN1,dirN2,ls_maj)


GPP_path = dirN2(1:36);

%j_maj = numéro de la liste_maj
if size(ls_maj,1) > 0
    
for j_maj=1:size(ls_maj(:,1)) %Boucle sur les dates de la liste_maj
X=[];Y=[];dateavg=[];Error = [];    
    try
    
    disp(ls_maj(j_maj,:));
    
    %Liste des images 'A' de la journée j_maj
    lsimg=ls([dirN1,ls_maj(j_maj,:),'\A*']);   
    
    for i=1:1:size(lsimg,1) %Boucle sur les images 'A' de la journée j_maj
        im_location = [dirN1,ls_maj(j_maj,:),'/',lsimg(i,:)];
        Img=imread(im_location);
        

        cmin=1/5;
        cmax=1/10;
        nbmin=5;


        [x,y]=DetectSeuil_V20130725_GPP(Img,cmin,cmax,'GPP',[GPP_path,'GPP_Toolbox\GPP_Functions\']);

        ii=find(x>0&y>0&isnan(x)==0&isnan(y)==0);kl=length(ii);

        if kl>40&length(find(abs(diff(x(ii)))>30))<1
            x2=x(ii);
            [fd gtf]=unique(y(ii));
            xi=interp1(fd,x2(gtf),1:size(Img,2));
            %resolution (m/pix along the stack)

            z=zeros(1,1600);

            gg=1:size(Img,2);

            [Xv,Yv]=PixtoCoordGPP(gg,xi,z);
            %[Xv,Yv]=PixtoCoordGPP('RectGPP',1:size(Img,2),xi');

            % Coordonnées video
            ang=90-10;
            ang=ang*pi/180;
            y0=370341;%Easting
            x0=694135;%Northing


            Co2=Xv;Co1=Yv;
            % Rotation et translation dans le nouveau repère (x0;y0) (pied de
            % l'échafaudage = référence)
            X1=(Co1-x0)*cos(ang)+(Co2-y0)*sin(ang);
            Y1=(Co2-y0)*cos(ang)-(Co1-x0)*sin(ang);


            %coordonnées shoreline
            X=[X X1'];
            Y=[Y Y1'];

            dateavg=[dateavg datenum(lsimg(i,3:end-4),'yyyymmddHHMM')];

        end    
    end

 

    catch
        warning('Problem computing PostProcess_GPP_ligne_eau_201511111 at 40')
        Error = 'Problem computing PostProcess_GPP_ligne_eau_201511111 at 40'
    end        


%Sauvegarde du fichier journalier
File_name=strcat('GPP_ligne-eau_',ls_maj(j_maj,:)); 
cd(dirN2);
save(File_name,'X','Y','dateavg','Error');

end
end
end




