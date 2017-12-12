% ComputeGPPCharLongTerm_20150127.m
%Fonction permettant de calculer les paramètres de vague
function CAMS_N2_Parameters(dirN1,dirN2_p,ls_maj_p,swash)
if length(ls_maj_p) > 0
for j_maj = 1 : length(ls_maj_p)
    lsrep = ls_maj_p(j_maj,:);
    disp(lsrep);
 try
files=[];
datef=[];

for i=1:size(lsrep,1)
    lsjour=ls([dirN1,lsrep(i,:),'/S_3_*']);
    for j=1:size(lsjour,1)
        files=[files [dirN1,lsrep(i,:),'/',lsjour(j,:)]'];
        datef=[datef datenum(files(end-15:end-4,end)','yyyymmddHHMM')];
    end
end

mondatef=str2num(datestr(datef,'yyyymm'));
months=unique(mondatef);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SystemVid=[370341 694135 10 -10];
dt = 0.5;
sk1=[130 520];
sk2=[360 1600];

% CoordCam = [1355261.683 303848.647 8] ;
%Def des coordonnées des stacks
sk1=double(sk1);
sk2=double(sk2);
lmax=max(abs(sk2(2)-sk1(2)),abs(sk2(1)-sk1(1)));
pas1=(sk2(1)-sk1(1))/lmax;
pas2=(sk2(2)-sk1(2))/lmax;
hor=round(sk1(1):pas1:sk2(1));
vert=round(sk1(2):pas2:sk2(2));


freq=2;%htz
periode=1.5;%s

for m=1:length(months)
    
    %Célérité des vagues (m/s)
    CT=[];
    %profondeur (m)
    depthT=[];
    %Hauteur significative (m)
    hsT=[];
    %Hauteur moyenne (m)
    hmT=[];
    %Période rms (m)
    TmT=[];
    %Periode peak (m)
    TpT=[];
    % Dissipation d'énergie due au déferlement (j)
    ErT=[];
    %Longueur roller (m)
    rolLT=[];
    % Nombre de vagues détectées sur la durée du stack
    nbwaveT=[];
    % Coordonées cross-shore (m)
    X1T=[ ];
    %Coordonées longshore (m)
    Y1T=[ ];
    %Resolution sur le stack (m/pix)
    dxT=[ ];
    
    BreakstdT=[ ];
    Breakmean1T=[ ];
    Breakmean2T=[ ];
    
    %Date
    dateT=[];
    %dadcp=[];
    %hadcp=[];
    %tadcp=[];
    %timeadcp=[];
    
    %Swash
    T=[];Tup=[];Tout=[];Tot=[];T_out=[];
    SwashTIn=[];
    SwashTOut=[];
    Asym=[];
    In=[];
    Ou=[];
    To=[];
    Stacktot=[];
    PosT=[];
    frac_up=[];
    frac_tot=[];
    datefrac=[];
    SwashT=[];
    time_inc=[];time_out=[];time_tot=[];
    Skew=[ ];
    ExtensionSWT1=[ ];
    ExtensionSWT2=[ ];
    TimeSwash=[];
    htiersOut=[ ];hrmsOut=[ ];trmsOut=[ ];
    htiersIn=[ ];hrmsIn=[ ];trmsIn=[ ];
    
    
    
    fmonth=find(mondatef==months(m));
    
    for fi=1:length(fmonth)
        
        datevid=datef(fmonth(fi));
        %          try
        S=imread(files(:,fi)');
        
        S=S(:,:,:);
%%%        
        [C,depth,hs,hm,tm,tp,Er,rolL,nbwave,X1,Y1,dx,Breakstd,Breakmean1,Breakmean2]=GPP_get_wave_parameters(S,['RectGPP.mat'],dt,sk1,sk2,SystemVid(1:3),SystemVid(4));
%%%        
        veci=-350:1:-25; %vecteur de projection des données
        clear idx
        for pk=1:length(veci)
            [frd gtf]=min(abs(veci(pk)-X1));
            idx(pk)=gtf;
        end
        CT=[CT C(idx)'];
        depthT=[depthT depth(idx)'];
        ErT=[ErT Er(idx)'];
        rolLT=[rolLT rolL(idx)'];
        hsT=[hsT hs];
        hmT=[hmT hm];
        TmT=[TmT tm];
        TpT=[TpT tp];
        dateT=[dateT datevid];
        nbwaveT=[nbwaveT nbwave'];
        X1T=[X1T X1(idx)'];
        Y1T=[Y1T Y1(idx)'];
        dxT=[dxT dx(idx)'];
        BreakstdT=[BreakstdT Breakstd(idx)'];
        Breakmean1T=[Breakmean1T Breakmean1(idx)'];
        Breakmean2T=[Breakmean2T Breakmean2(idx)'];
        [frd ft]=min(abs(date-datevid));
        
           
        %Swash%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if swash == 1
            
        I=S;        
        tims=datenum(files([end-15:end-4],fmonth(fi))','yyyymmddHHMM');
        
        vec=slidefun(@mean,20,slidefun(@mean,50,(nanmean(double(I(:,:,1)))./(0.5*nanmean(double(I(:,:,2))+double(I(:,:,3)))))));
        vecs=sort(vec);mx=nanmedian(vecs( round(3*length(vec)./4:length(vec))));mn=nanmedian(vecs( round(1:1*length(vec)./4)));
        op=find(vec>mn+0.66*(mx-mn));ShorelinePos=op(1);%swash location (in pix)
        ii=find(vec>mn+0.25*(mx-mn)&vec<mn+0.75*(mx-mn));id=find(diff(ii)>1);
        if numel(id)>1
            widthSZ=length(find(vec(1:ii(id(2)))>mn+0.25*(mx-mn)&vec(1:ii(id(2)))<mn+0.75*(mx-mn)));%Swash length (in pix)
        else
            widthSZ=length(find(vec>mn+0.25*(mx-mn)&vec<mn+0.75*(mx-mn)));%Swash length (in pix)
            % [frd gtf]=min(diff(smooth(smooth(nanstd(detrend(detrend(double(rgb2gray(I(:,:,:))))')'),size(I,2)./5),20)))
        end
        for is=1:1
            ind=(is-1)*(floor(mean(diff(linspace(1,size(I,1),2)))))+(1:floor(mean(diff(linspace(1,size(I,1),2)))));
            A=detrend(detrend(double(rgb2gray(I(ind,max([1 ShorelinePos-widthSZ]):min([size(I,2) ShorelinePos+widthSZ]),:))))')';
            
            [Sin,Sout]=RadonSeparation_filt(A);
            
            
            Breakstd=FiltreMean(FiltreMean(nanstd(Sout(round(size(Sout,1)/20:19*size(Sout,1)/20),:)),round(size(Sout,2)/50)),round(size(Sout,2)/30));
            %Computation of threshold (thresh)  %Normalisation de std
            Breakstd=(Breakstd-nanmin(Breakstd))./nanmax(Breakstd-nanmin(Breakstd));
            Sout=Sout.*repmat(Breakstd,size(Sout,1),1);
            
            
            Sin1=Sin;
            Sout1=Sout;
            tot=detrend(std(double(Sout1+Sin1)'));tot=tot-smooth(tot,60)';
            
            
            
            In=[In std(double(Sout1))];
            Ou=[Ou std(double(Sin1))];
            To=[To tot];
            
            [frd mxOut]=max(Sout');
            SwashTIn=[SwashTIn mxOut];
            [frd mxIn]=max(Sin');
            SwashTOut=[SwashTOut mxIn];
            
            filt=0;% : filtrage (1) /pas de filtre (0)
            meth=1;% mean zero crossing
            [hs,htiers,hrms,trms,Hmax,h,Tp,t]=Wave_Char(SwashTOut,1/freq,filt,meth);
            htiersOut=[htiersOut htiers];hrmsOut=[hrmsOut hrms];trmsOut=[trmsOut trms];
            [hs,htiers,hrms,trms,Hmax,h,Tp,t]=Wave_Char(SwashTIn,1/freq,filt,meth);
            htiersIn=[htiersIn htiers];hrmsIn=[hrmsIn hrms];trmsIn=[trmsIn trms];
            
            
            Asym=[Asym asym(mxOut)];
            Skew=[Skew skew(mxOut)];
            ExtensionSWT1=[ExtensionSWT1 widthSZ];
            ExtensionSWT2=[ExtensionSWT2 4.*nanstd(mxOut)];
            tps=tims+(1:length(mxOut))./(24*3600*2);
            TimeSwash=[TimeSwash tps];
            Stacktot=[Stacktot imrotate(I(ind(1:1),:,:),90)];
            PosT=[PosT ShorelinePos];
            %end swash
            
        end
        end
    end
    
    cd(dirN2_p);
    File_name=strcat('GPP_Parameters_',ls_maj_p(j_maj,:))    
    save(File_name,'hsT','CT','depthT','ErT','rolLT','hmT','TmT',...
        'TpT','dateT','nbwaveT','X1T','Y1T','dxT','BreakstdT','Breakmean1T',...
        'Breakmean2T');
    
    

end


     catch
        
    cd(dirN2_p);
    File_name=strcat('GPP_Parameters_',ls_maj_p(j_maj,:))    
    save(File_name,'hsT','CT','depthT','ErT','rolLT','hmT','TmT',...
        'TpT','dateT','nbwaveT','X1T','Y1T','dxT','BreakstdT','Breakmean1T',...
        'Breakmean2T');

        
    end

end
end
end





