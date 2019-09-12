clearvars; close all; clc; 

%%
dx =0.01 ; 
[xg,yg]=meshgrid(-1.0:dx:1.0,-1.0:dx:1.0);
x0y0(1) = min(xg(:)); 
x0y0(2) = min(yg(:)); 

nrows = size(xg,1) ; 
ncols = size(xg,2) ; 

%% ISOTROPIC MESH SIZE 
% circular resolution
sz = xg.^2 + yg.^2 ;

% prevent size going to 0 
sz(sz < dx ) = dx ;

figure;
pcolor(xg,yg,sz); 
shading interp
title('Mesh sizes') ; 
colorbar

fid = fopen('MeshSizes.txt','w') ; 
fprintf(fid,'%d %d %f %f %f\n',nrows,ncols,dx,x0y0(1),x0y0(2)) ; 
for i = 1 : nrows
    for j = 1 : ncols
       fprintf(fid,'%f ',sz(i,j)) ; 
    end
    fprintf(fid,'\n') ; 
end
fclose(fid) ;

%%
sz = sin(pi*xg)*10 ;
% prevent elongation from going below 1
sz = sz - repmat(min(sz(:)),nrows,ncols) + repmat(1,nrows,ncols) ; 

figure ; 
pcolor(xg,yg,sz) ;
shading interp
title('Elongation') ; 
colorbar; 

fid = fopen('ElongSizes.txt','w') ; 
fprintf(fid,'%d %d %f %f %f\n',nrows,ncols,dx,x0y0(1),x0y0(2)) ; 
for i = 1 : nrows
    for j = 1 : ncols 
       fprintf(fid,'%f ', 0); %sz(i,j)); 
    end
    fprintf(fid,'\n') ; 
end
fclose(fid) ;

%% ANGLE OF ELONGATION

sz = sin(pi*xg)*40 ;

figure ;
pcolor(xg,yg,sz) ;
shading interp
title('Angles') ;

fid = fopen('Angles.txt','w') ;
fprintf(fid,'%d %d %f %f %f\n',nrows,ncols,dx,x0y0(1),x0y0(2)) ; 
for i = 1 : nrows
    for j = 1 : ncols 
       fprintf(fid,'%f ', 0); % sz(i,j)) ; 
    end
    fprintf(fid,'\n') ; 
end
fclose(fid) ;
