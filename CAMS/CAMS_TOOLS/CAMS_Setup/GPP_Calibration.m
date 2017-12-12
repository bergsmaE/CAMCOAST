% % Test.m
% 
% %Buscar puntos de contro en google earth
% % cambiar coordenadas con "deg2utm.m"
% %ex
% 
% for i=1:19
% [x,y,utmzone] = deg2utm(A(i,1),A(i,2));
% vv(i,:)=[x y];
% end
% %Llenar el archivo de coordenadas reales > pixeles
% %CamGPP.txt
% %ex :
% % 251282	6305197	0	125	97
% % lat	lon	vert	pixvert	pixhor
% 
% %correr la siguiente funcion
% % [err]=Findpos_Stoyanov_GPP(focale,dist1,dist2,ImageGPP,PointsControl);
% clear E
% for i=1:1:3
% [err]=Findpos_Stoyanov_GPP(2.5,0,0,'ImgGPP.jpg','GPP_GPS_Points_20130225.txt');
% E(i)=err;
% i
% end
% %Buscar los parametros optimales: Focale, distorsion 1, distorsion2
% %Que permiten reducir el error "err".
% % a traves de algo iterativo
% %ex: for focale=... end
% 
% % Definicion de la zona de interes
% ZoneGPP(ImageGPP,res)
% %res: resolucion (1)
ZoneGPP('ImgGPP.jpg',1);

% Rectificacion de la imagen
% [Arect,X,Y]=RectifyImageGPP(ImageGPP,Rectfile,Zonefile,GridRes);
%ex:
[Arect,X,Y]=RectifyImageGPP('ImgGPP.jpg','RectGPP','ZoneGPP',0.1);
%Arect : imagen rectificada
%X,Y, matrices de coordenadas

% Resolucion sobre un perfil
%[X,Y,s1,s2, res]=ResStackGPP(ImageGPP,Rectfile,sk1,sk2);
%sk1, 2 son coordenadas en pixeles de puntos extremos del perfil
%Si null, se hace a la mano con ginput
%ex: sk1=[174 702];sk2=[15 84];
%ex:[X,Y,s1,s2,ResStack]=ResStackGPP('ImageGPP.jpg','RectGPP',[],[]);