clearvars; close all; clc; 

%%
dx =0.05; 
[xg,yg]=meshgrid(-2.0:dx:2.0,-2.0:dx:2.0);
x0y0(1) = min(xg(:)); 
x0y0(2) = min(yg(:)); 

nrows = size(xg,1) ; 
ncols = size(xg,2) ; 

PSLG = [min(xg(:)),max(yg(:))
        max(xg(:)),max(yg(:)) 
        max(xg(:)),min(yg(:))
        min(xg(:)),min(yg(:))
        min(xg(:)),max(yg(:))]; 
    
%hold on; plot(PSLG(:,1),PSLG(:,2),'k-','linewi',2) ;

fid = fopen('PSLG.txt','w') ;
fprintf(fid,'%d %d\n',length(PSLG),2) ; 
for i = 1 : length(PSLG) 
   fprintf(fid,'%f %f \n',PSLG(i,1),PSLG(i,2)) ;
end

%% MESH SIZE IN X-DIRECTION 
szx = 0.0005 + 1.5*abs(1-(xg.^2  + yg.^2 ).^0.5);  
szy = 0.1*(xg.^2 + yg.^2).^0.5 + 1.5*abs(1-(xg.^2  + yg.^2 ).^0.5);  
% szy = zeros(nrows,ncols)+dx ; 
% elong =10; 
% grade =3; 
% 
% x=-1:dx:1; 
% %x=y.^2 ;
% y=linspace(-1.0,1.0,length(x)); 
% line = [x',y'];
% testset = [xg(:),yg(:)];
% [~, dist] = ourKNNsearch(line',testset',1) ;
% dist = reshape(dist,nrows,ncols); 
% 
% szx = elong*dx - dist*grade ;
% szx(szx < dx) = dx ; 

%szy = dx + dist*grade; 
%szy(szy > elong*dx) = elong*dx; 

figure;
pcolor(xg,yg,szx);
shading interp
title('Mesh size in x-dir') ; 
colorbar

figure;
pcolor(xg,yg,szy); 
shading interp
title('Mesh size in y-dir') ; 
colorbar

fid = fopen('MeshSizeXdir.txt','w') ; 
fprintf(fid,'%d %d %f %f %f\n',nrows,ncols,dx,x0y0(1),x0y0(2)) ; 
for i = 1 : ncols
    for j = 1 : nrows
       fprintf(fid,'%f ', szx(j,i)) ; 
    end
    fprintf(fid,'\n') ; 
end
fclose(fid) ;

%% MESH SIZE IN Y-DIR
fid = fopen('MeshSizeYdir.txt','w') ;
fprintf(fid,'%d %d %f %f %f\n',nrows,ncols,dx,x0y0(1),x0y0(2)) ;
for i = 1 : ncols
    for j = 1 : nrows
        fprintf(fid,'%f ',szy(j,i));
    end
    fprintf(fid,'\n') ;
end
fclose(fid) ;
%% ANGLE OF ELONGATION (in radians) 
fid = fopen('Angles.txt','w') ;
fprintf(fid,'%d %d %f %f %f\n',nrows,ncols,dx,x0y0(1),x0y0(2)) ;

angles = zeros(nrows,ncols) ;

for i = 1 : ncols
    for j = 1 : nrows
        angles(j,i)  = atand(xg(j,i)/yg(j,i)) + 90  ;
    end
end
angles(isnan(angles))=0.0; 
for i = 1 : ncols
    for j = 1 :  nrows
        fprintf(fid,'%f ',angles(j,i)*(pi/180.0));
    end
    fprintf(fid,'\n') ;
end
fclose(fid) ;

figure;
pcolor(xg,yg,angles); 
shading interp
title('Mesh angles') ; 
colorbar
