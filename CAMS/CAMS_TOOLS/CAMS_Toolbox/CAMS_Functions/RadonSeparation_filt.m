% RadonSeparation_filt.m
function [Sin,Sout]=RadonSeparation_filt(M)
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
end

% Pre-Traitement
for i=1:nx
	M(1:nt,i)=M(1:nt,i)-nanmean(M(1:nt,i));
end

%Taille de l'image
[c1 c2 c]=size(M);
    
%Transformation de Radon
R=radon(double(M(:,:,1)),0:179);

%Filtre temporel
nr=size(R,1);
amp=nr/c1;
iang=1:1:180;
dt=0.5;

k=nt;%-round(nt./2);
% trk=round(nr/2+tr*((k-(nr+1)/2))*amp);
trk=floor((size(M,1))/2)-floor((  0*cosd(iang)+  ((size(M,1))/2-k*amp)*sind(iang)));trk=trk-min(trk);
res=(nt*dt)./(trk.*2);
for i=iang
   R(:,i)=smooth(R(:,i)-smooth(R(:,i),round(1+(18)./(res(i)))),round(max([1 1./res(i)]))); 
end
%Fin filtre temporel

%Incident
Lm=1;Lx=89;
%On reconstruit l'image 
I = iradon(R(:,Lm:Lx),(Lm:Lx)-1,'linear','Hann',c1); 
%Indices reellement filtres
mx=min(size(I));
mn=min(size(M));
ind=1+round(mx/2-mn/2):round(mx/2-mn/2)+mn;

%Outputs
Sin=M.*NaN;
if c2>c1
Sin(:,1+(c2-c1)/2:c2-(c2-c1)/2)=I(:,ind);
else
Sin(:,:)=I(:,ind);
    
end
Sin=Sin*0.5;
disp('Radon ok: composante incidente extraite')


%Incident
Lm=91;Lx=179;
%On reconstruit l'image 
I = iradon(R(:,Lm:Lx),(Lm:Lx)-1,'linear','Hann',c1); 
%Indices reellement filtres
mx=min(size(I));
mn=min(size(M));
ind=1+round(mx/2-mn/2):round(mx/2-mn/2)+mn;

%Outputs
Sout=M.*NaN;
if c2>c1
Sout(:,1+(c2-c1)/2:c2-(c2-c1)/2)=I(:,ind);
else
Sout(:,:)=I(:,ind);
    
end
Sout=Sout*0.5;
disp('Radon ok: composante réfléchie extraite')
    
% Sin=Sin.*((abs(Sin)-abs(Sout))./abs(Sin));
% Sout=Sout.*((abs(Sin)-abs(Sout))./abs(Sout));
% 
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

%%%%%%%%%%%Filtre%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nt=size(S,1);
nr=size(R,1);
amp=nr/nt;
iang=1:1:180;
dt=1;
k=nt;%-round(nt./2);
% trk=round(nr/2+tr*((k-(nr+1)/2))*amp);
trk=floor(nt/2)-floor((  0*cosd(iang)+  (nt/2-k*amp)*sind(iang)));trk=trk-min(trk);
res=(nt*dt)./(trk.*2);

for i=iang
   R(:,i)=smooth(R(:,i)-smooth(R(:,i),round(1+18./(res(i)))),3); 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


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