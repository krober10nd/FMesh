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
sz = xg.^2 / 2 ; 

% prevent size going to 0 
sz(sz < dx ) = dx ;
sz(sz > 5*dx) = 5*dx; 

% figure;
% pcolor(xg,yg,sz); 
% shading interp
% title('Mesh sizes') ; 
% colorbar

fid = fopen('MeshSizes.txt','w') ; 
fprintf(fid,'%d %d %f %f %f\n',nrows,ncols,dx,x0y0(1),x0y0(2)) ; 
for i = 1 : nrows
    for j = 1 : ncols
       fprintf(fid,'%f ', 2*dx) ; %sz(j,i)) ; 
    end
    fprintf(fid,'\n') ; 
end
fclose(fid) ;

%%
% Next create the circle in the image.
innerRadius = 0.25;
outerRadius = 0.50;
array2D = (xg - 0).^2 ...
    + (yg - 0).^2;

ringPixels = array2D >= innerRadius.^2 & array2D <= outerRadius.^2;
ringPixels = double(ringPixels) ;
ringPixels(ringPixels==1) = 20  ;
ringPixels(ringPixels==0) = 1 ;


dx2 = repmat(dx,[1,size(xg,2)]) ; 
hfun = zeros(nrows*ncols,1);
nn = 0;
for ipos = 1 : ncols
    for jpos = 1 : nrows
        nn = nn + 1;
        hfun(nn,1) = ringPixels(ipos,jpos);
    end
end
fdfdx=0.15 ; 
[hfun,flag] = limgradStruct(nrows,dx2,dx,hfun,...
    fdfdx,sqrt(length(hfun)));
% reshape it back
nn = 0;
for ipos = 1 : ncols
    for jpos = 1 : nrows
        nn = nn+1;
        ringPixels2(ipos,jpos) = hfun(nn);
    end
end

figure ;
pcolor(xg,yg,ringPixels2) ;
shading interp
title('Elong')

fid = fopen('ElongSizes.txt','w') ; 
fprintf(fid,'%d %d %f %f %f\n',nrows,ncols,dx,x0y0(1),x0y0(2)) ; 
for i = 1 : nrows
    for j = 1 : ncols 
       fprintf(fid,'%f ', ringPixels(j,i)); 
    end
    fprintf(fid,'\n') ; 
end
fclose(fid) ;

%% ANGLE OF ELONGATION

fid = fopen('Angles.txt','w') ;
fprintf(fid,'%d %d %f %f %f\n',nrows,ncols,dx,x0y0(1),x0y0(2)) ; 
for i = 1 : nrows
    for j = 1 : ncols 
       fprintf(fid,'%f ', 90*(pi/180)); % sz(i,j)) ; 
    end
    fprintf(fid,'\n') ; 
end
fclose(fid) ;


