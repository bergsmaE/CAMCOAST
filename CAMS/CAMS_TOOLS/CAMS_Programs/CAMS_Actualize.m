%Actualiser
clear;
%Ce spcript permet d'actualiser la banque de données du système video
%A chaques nouvelles images sauvegardées dans le dossier
%du NIVEAU 1 de la banque de données, le script détecte 
%ces nouvelles images et lance les programmmes permettant d'obtenir:
%
% NIVEAU 2 de la banque de données
%   - les coordonnées de la ligne d'eau via CAMS_N2_Shoreline.mat
%   - les paramètres de vague via CAMS_N2_parameters.mat
%   - les direcions des vagues via CAMS_N2_direction.mat
%
% NIVEAU 3 de la banque de données
%   - les moyennes journalières et mensuelles de HS,TP,Dir,Ligne_eau
%   - les graphes associés

%% CONFIGURATION:

SITENAME = 'GRANDPOPO';
SITECODE= 'GPP';

%Emplacement du dossier "CAMS_DATA":
dirdata = 'D:\Tout le stage\CAMS\CAMS_DATA\';

%Emplacement du dossier "CAMS_TOOLS":
dirtools = 'D:\Tout le stage\CAMS\CAMS_TOOLS\';

%% IDENTIFICATION DES CHEMINS DES DOSSIERS DE LA BASE DE DONNEES

%Ajout des programmes de CAMS_TOOLS au Path de Matlab:
addpath(genpath(dirtools));
%Chemin du dossier NIVEAU 1 de la base de données: dossier images
dirN1 = [dirdata '02 - DATA\' SITENAME '\' SITECODE '_NIVEAU 1\' ];
%Chemin du dossier NIVEAU 2
dirN2 = [dirdata '02 - DATA\' SITENAME '\' SITECODE '_NIVEAU 2\' ];
%Chemin du dossier NIVEAU 3
dirN3 = [dirdata '02 - DATA\' SITENAME '\' SITECODE '_NIVEAU 3\' ];


%% ACTUALISATION DU NIVEAU 2
disp('ACTUALISATION DU NIVEAU 2')
%SHORELINE
% disp('CALCUL DES POSITIONS DU TRAIT DE COTE')
  dirN2_s = [dirN2 SITECODE '_Shoreline_Data\'];

%On inventorie les dates pour lesquelles le trait de côte n'a pas été
%calculé:
[ls_maj_s] = MAJBDD1(dirN1,dirN2_s);

%On lance le calcul du trait de côte pour les jours listés dans ls_maj_s 
CAMS_N2_Shoreline(dirN1,dirN2_s,ls_maj_s)
disp('trait de côte à jour')

%LES PARAMETRES DE VAGUE
disp('CALCUL DES PARAMETRES DE VAGUE')
 dirN2_p = [dirN2 SITECODE '_Parameters_Data\'];
[ls_maj_p] = MAJBDD1(dirN1,dirN2_p);
CAMS_N2_Parameters(dirN1,dirN2_p,ls_maj_p,0)
disp('Paramètres de vague à jour');


%LA DIRECTION DES VAGUES
disp('CALCUL DE LA DIRECTION DES VAGUES')
 dirN2_d = [dirN2 SITECODE '_Direction_Data\'];
[ls_maj_d] = MAJBDD1(dirN1,dirN2_d);
CAMS_N2_Direction(dirN1,dirN2_d,ls_maj_d);
disp('Directions à jour');

%% ACTUALISATION DU NIVEAU 3
disp('ACTUALISATION DU NIVEAU 3')
CAMS_N3_Parameters(dirN2_p,dirN2_d,dirN2_s,dirN3);


