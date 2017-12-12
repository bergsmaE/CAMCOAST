% MAJBDD1
% Comparaison des dates du produit de niveau 1 (images) avc celles du
% produit de niveau 2 (ligne eau, parametres).
% Elaboration d'une liste des dates des images qui n'ont pas encore été
% tratées: date_maj

function [liste_maj] = MAJBDD1(dirN1,dirN2);
%liste des dates des images, NIVEAU 1, format:'YYYYMMDD' 
lsN1=str2num(ls([dirN1,'*20*']));

%liste des dates du NIVEAU 2, format:'YYYYMMDD' 
lsN2= ls([dirN2,'*_20*']);
if length(lsN2)>0
    lsN2 = str2num(lsN2(:,[end-11:end-4]));
else
    lsN2 = [];
end

    %Comparaison des listes
    liste_maj = [];

    for i=1:length(lsN1)
        if length(find(lsN2 == lsN1(i)))==0;
            liste_maj = [liste_maj lsN1(i)];
        end
    end

    liste_maj = num2str(sort(liste_maj)');

end


