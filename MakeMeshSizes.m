clearvars; close all; clc; 

%%
dx =0.05; 
[xg,yg]=meshgrid(-1.0:dx:1.0,-1.0:dx:1.0);
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
szy = zeros(nrows,ncols)+dx ; 

y=-1:dx:1; 
x=y.^2 ;
line = [x',y'];
testset = [xg(:),yg(:)];
[~, dist] = ourKNNsearch(line',testset',1) ;
dist = reshape(dist,nrows,ncols); 

szx = 5*dx - dist*0.35 ;


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
for i = 1 : nrows
    for j = 1 : ncols
       fprintf(fid,'%f ', szx(i,j)) ; 
    end
    fprintf(fid,'\n') ; 
end
fclose(fid) ;

%% MESH SIZE IN Y-DIR
fid = fopen('MeshSizeYdir.txt','w') ;
fprintf(fid,'%d %d %f %f %f\n',nrows,ncols,dx,x0y0(1),x0y0(2)) ;
for i = 1 : nrows
    for j = 1 : ncols
        fprintf(fid,'%f ',szy(i,j));
    end
    fprintf(fid,'\n') ;
end
fclose(fid) ;
%% ANGLE OF ELONGATION (in radians) 
fid = fopen('Angles.txt','w') ;
fprintf(fid,'%d %d %f %f %f\n',nrows,ncols,dx,x0y0(1),x0y0(2)) ;

angles = zeros(nrows,ncols) ;

for i = 1 : nrows
    for j = 1 : ncols
        angles(i,j)  = atand(yg(i,j)/-xg(i,j)) ;
    end
end
angles(isnan(angles))=0.0; 
for i = 1 : nrows
    for j = 1 :  ncols
        fprintf(fid,'%f ', 0*angles(i,j)*(pi/180.0)); %(angles(j,i))*(pi/180)); % sz(i,j)) ;
    end
    fprintf(fid,'\n') ;
end
fclose(fid) ;

figure;
pcolor(xg,yg,angles); 
shading interp
title('Mesh angles') ; 
colorbar
