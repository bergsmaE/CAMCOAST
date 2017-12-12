% Aquisition_20141216.m

root='C:\Users\HP\Desktop\GrandPopo_Data\04-GrandPopo_videos_manual_processing\VidéosFilmsCorrompusAobserver - Copie/';%Input files from camera
rootimg='C:\Users\HP\Desktop\GrandPopo_Data\01-GrandPopo_Images\01.1-GrandPopo_Images\TESTS/';%Output files

duree=15;% Image duration (defaut is 15 min)
interval=15;% time between two images beginning (defaut is 15 min)

stackn=3;%Stack number

warning off
% Pour compiler: mcc -mv Aquisition_20130104.m
while 1==1
    nfiles=0;
    
    while nfiles<duree
        try
            L1=ls([root,'/20*']);
            file =[];
            for day=1:size(L1,1)
                Lday=ls([root,L1(day,:),'/*']);
                
                try
                    if size(Lday,1)<3; rmdir([root,L1(day,:)]) ; end
                end
                
                for hr=1:size(Lday,1)
                    Lhr=ls([root,L1(day,:),'/',Lday(hr,:),'/*.mp4']);
                    
                    try
                        if mean(size(Lhr))==0; rmdir([root,L1(day,:),'/',Lday(hr,:)]) ; end
                    end
                    
                    for mn=1:size(Lhr,1)
                        try
                            fileInfo = dir([root,L1(day,:),'/',Lday(hr,:),'/',Lhr(mn,1:6)]);
                            if strcmp(Lhr(mn,1),'t')==0 & strcmp(Lhr(mn,3),'_')==0 & fileInfo.bytes>0
                                file=[file [root,L1(day,:),'/',Lday(hr,:),'/',Lhr(mn,1:6)]'];
                            end
                        end
                    end
                end
            end
            disp(['          Waiting for camera video files '])
            pause(5)
            
            %Gestion des video corrompues
            if exist([root,'crash.mat'])~=0
                load([root,'crash'],'datelast')
                dates=datenum(file(size(file,1)-[17:-1:10 8:-1:7 5:-1:4],:)','yyyymmddHHMM');
                filesok=find(dates>datelast+3/(24*60));
                for fbad=1:size(file,2)-length(filesok)
                    delete(file(:,fbad)');%remove file
                    disp(['Img ',file(:,fbad)',' removed '])
                end
                try
                    file=file(:,filesok);
                catch
                    file=[];
                end
            end
            
            nfiles=size(file,2);
        end
    end%nfiles
    
    
    try;disp(['Number of videos to be processed :  ',num2str(nfiles)]);end
    dates=datenum(file(size(file,1)-[17:-1:10 8:-1:7 5:-1:4],:)','yyyymmddHHMM');
    
    
    for sn=1:stackn;
        if sn==1
            sk1=[150 300];
            sk2=[550 900];
        elseif sn==2
            sk1=[350 1];
            sk2=[160 1600];
        elseif sn==3
            sk1=[130 520];
            sk2=[360 1600];
            
        end
        %Def des coordonnées des stacks
        lmax=max(abs(sk2(2)-sk1(2)),abs(sk2(1)-sk1(1)));
        pas1=(sk2(1)-sk1(1))/lmax;
        % if sk1(1)>sk2(1) pas1=-pas1; end
        pas2=(sk2(2)-sk1(2))/lmax;
        % if sk1(2)>sk2(2) pas2=-pas2; end
        eval(['s1_',num2str(sn),'=round(sk1(1):pas1:sk2(1));s2_',num2str(sn),'=round(sk1(2):pas2:sk2(2));']);
    end
    datesmin=round(min(dates)/ (interval/(24*60)))*(interval/(24*60)):interval/(24*60):round(max(dates)/ (interval/(24*60)))*(interval/(24*60));
    for t=1:length(datesmin)
        try
            tic
            ind=find(dates-datesmin(t)>=0&dates-datesmin(t)<(duree-0.01)/(24*60));
            if length(ind)==duree
                for sn=1:stackn;
                    eval(['Sr_',num2str(sn),'=NaN*ones(10000,length(s1_',num2str(sn),'));Sg_',num2str(sn),'=NaN*ones(10000,length(s1_',num2str(sn),'));Sb_',num2str(sn),'=NaN*ones(10000,length(s1_',num2str(sn),')); clear S']);
                end
                cc=0;
                for i=1:length(ind)
                    %      vp = VideoPlayer(file(:,ind(i))', 'Verbose', false, 'ShowTime', false);
                    disp(['          Generating Img ',datestr(datesmin(t)),'   File n°', num2str(i),'/',num2str(length(ind)),' (',datestr(dates(ind(i))),' )'])
                    
                    datelast=dates(ind(end));%dernière vidéo ouverte
                    save([root,'crash'],'datelast');
                    
                    frameid=0;
                    while true
                        try
                            %           vp.nextFrame;
                            %           A=255*(vp.Frame);
                            frameid=frameid+1;
                            video = mmread(file(:,ind(i))',frameid);
                            A=video.frames.cdata;
                            cc=cc+1;
                            for sn=1:stackn;
                                eval(['StempR=([diag(A(s1_',num2str(sn),',s2_',num2str(sn),',1))]'')'';StempG=([diag(A(s1_',num2str(sn),',s2_',num2str(sn),',2))]'')'';StempB=([diag(A(s1_',num2str(sn),',s2_',num2str(sn),',3))]'')'';'])
                                eval(['Sr_',num2str(sn),'(cc,:)=StempR;Sg_',num2str(sn),'(cc,:)=StempG;Sb_',num2str(sn),'(cc,:)=StempB;']);
                            end
                            if cc==1
                                [n, m,c]=size(A);Amoy=zeros(n,m,c);
                            else
                                Amoy = Amoy+double(A);
                            end
                        catch
                            break
                        end
                    end
                    clear vp
                end
                
                
                
                %Concatenation des couleurs
                for sn=1:stackn;
                    eval(['Sr_',num2str(sn),'(find(Sr_',num2str(sn),' > 255))=255;Sr_',num2str(sn),'(find(Sr_',num2str(sn),' < 0))=0;Sg_',num2str(sn),'(find(Sg_',num2str(sn),' > 255))=255;Sg_',num2str(sn),'(find(Sg_',num2str(sn),' < 0))=0;Sb_',num2str(sn),'(find(Sb_',num2str(sn),' > 255))=255;Sb_',num2str(sn),'(find(Sb_',num2str(sn),' < 0))=0;']);
                    eval(['S_',num2str(sn),'=uint8(ones(cc,length(s1_',num2str(sn),'),3));S_',num2str(sn),'(:,:,1)=uint8(Sr_',num2str(sn),'(1:cc,:));S_',num2str(sn),'(:,:,2)=uint8(Sg_',num2str(sn),'(1:cc,:));S_',num2str(sn),'(:,:,3)=uint8(Sb_',num2str(sn),'(1:cc,:));']);
                end
                Amoy=Amoy./cc;
                Amoy(find(Amoy>255))=255;
                Amoy(find(Amoy<0))=0;
                Amoy = uint8(Amoy);
                mkdir([rootimg,datestr(datesmin(t),'yyyymmdd')])
                for sn=1:stackn;
                    eval(['imwrite(S_',num2str(sn),',[[rootimg,datestr(datesmin(t),''yyyymmdd'')],''/S_',num2str(sn),'_'',datestr(datesmin(t),''yyyymmddHHMM''),''.jpg''],''jpg'',''Quality'',100);']);
                end
                imwrite(Amoy,[[rootimg,datestr(datesmin(t),'yyyymmdd')],'/A_',datestr(datesmin(t),'yyyymmddHHMM'),'.jpg'],'jpg','Quality',100);
                imwrite(uint8(A),[[rootimg,datestr(datesmin(t),'yyyymmdd')],'/I_',datestr(datesmin(t),'yyyymmddHHMM'),'.jpg'],'jpg','Quality',100);
                [[rootimg,datestr(datesmin(t),'yyyymmdd')],'/I_',datestr(datesmin(t),'yyyymmddHHMM'),'.jpg']
                beep
                toc
                disp(['Img ',datestr(datesmin(t)),' Generated '])
                
                delete([root,'crash.mat']);
                
                for i=[1:min(ind) ind']
                    try
                        % mkdir([rootvid,file(length(root)+1:length(root)+11,i)'])
                        % copyfile(file(:,i)',[rootvid,file(length(root)+1:length(file(:,i)),i)'])
                        delete(file(:,i)');%remove file
                        disp(['Img ',file(:,i)',' copied and removed '])
                    end
                end
                
                
            end%if
            
            
        end%try
    end%number of file==duree
    
    if str2num(datestr(datenum(clock),'HH'))>20|str2num(datestr(datenum(clock),'HH'))<4
        reste=ls([root,'normal*']);
        for re=1:size(reste,1)
            delete([root,reste(re,:)])
        end
    end
    
end