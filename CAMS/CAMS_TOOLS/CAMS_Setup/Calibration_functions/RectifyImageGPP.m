%  RectifyImageGPP.m

function [Arect,X,Y]=RectifyImageGPP(ImageGPP,Rectfile,Zonefile,GridRes)





method = 'linear'     % Method for griddata


      A=imread(ImageGPP); 

    % Recupération de la resolution : resx et resy et des indices des points de la zone
    % : indo et des points de la zone : x1, y1
    load(Rectfile);
    load(Zonefile);

    % Application de la zone et passage du tableau en vecteur (1 colonne) pour les 
    % 3 composantes rouge, vert, bleu 

    nx=size(A,1);
    ny=size(A,2);
    rr=double(A(x1(:)+(y1(:)-1)*nx)');
    gg=double(A(x1(:)+(y1(:)-1)*nx+nx*ny)');
    bb=double(A(x1(:)+(y1(:)-1)*nx+2*nx*ny)');                  


    % Rectification des coordonnées
    disp('Rectification de l''image...')
    load(Rectfile);

    [Vcoord1,Vcoord2]=PixtoCoordGPP(y1,x1,ones(1,length(y1)));
X_I(1,:)=Vcoord1;
X_I(2,:)=Vcoord2;
    yr1=X_I(1,:);
    xr1=X_I(2,:);



  % Merging des images
  clear Arect
  disp(['Elaboration de l''image complète...'])  

    load(Rectfile);
  
  % Résolution finale de la figure et nouvelle grille
%       GridRes=1;
      xr=min(X_kk(1,:)):GridRes:max(X_kk(1,:));
      yr=min(X_kk(2,:))-10:GridRes:max(X_kk(2,:));
  [X,Y] = meshgrid(xr,yr);

      xx=yr1;
      yy=xr1;
      red   = rr;
      green = gg;
      blue  = bb;
  
  % Pack
  Arect(:,:,1) = griddata(xx,yy,red,X,Y,method);
  Arect(:,:,2) = griddata(xx,yy,green,X,Y,method);
  Arect(:,:,3) = griddata(xx,yy,blue,X,Y,method);

  % Contrôle des limites de bandes
  Arect(find(Arect>255))=255;
  Arect(find(Arect<0))=0;

  % Passage en uint8
  Arect = uint8(Arect);

  fig=figure(1);
  clf
  set(fig,'doubleBuffer','on','Position',[50 25 600 650]);
  image(X(1,:)-max(X(1,:)),Y(:,1)-max(Y(:,1)),Arect)
  axis image
  axis xy
  set(gcf,'Color','w')
  set(gca,'FontSize',12)
  xlabel ('Cross-shore position (m)')
  ylabel ('Longshore position (m)')
  grid on

  
end
%   Ecriture de l'image
%  hgsave(fig,nom);
%  AA= frame2im(getframe(fig));
%  imwrite(AA,[chemin,'Rectified\GPP_Zone.png'],'png');

