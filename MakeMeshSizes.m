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
       fprintf(fid,'%f ', 5*dx) ; %sz(j,i)) ; 
    end
    fprintf(fid,'\n') ; 
end
fclose(fid) ;

%%
% Next create the circle in the image.
innerRadius = 0.40;
outerRadius = 0.50;
array2D = (xg - 0).^2 ...
    + (yg - 0).^2;
elong = 8 ; 
ringPixels = array2D >= innerRadius.^2 & array2D <= outerRadius.^2;
ringPixels = double(ringPixels) ;
ringPixels(ringPixels==1) = elong  ;
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
fdfdx=20.0 ; 
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
colorbar ;
shading interp
title('Elong')

fid = fopen('ElongSizes.txt','w') ; 
fprintf(fid,'%d %d %f %f %f\n',nrows,ncols,dx,x0y0(1),x0y0(2)) ; 
for i = 1 : ncols
    for j = 1 : nrows 
       fprintf(fid,'%f ',  ringPixels(j,i)); 
    end
    fprintf(fid,'\n') ; 
end
fclose(fid) ;

%% ANGLE OF ELONGATION
angles=zeros(nrows,ncols) ;
for i = 1 : ncols
    for j = 1 : nrows
        % if in ring
        if(ringPixels(i,j)==elong)
            if(xg(i,j) < 0)
                angles(i,j) = 270 - atand(yg(i,j)/-xg(i,j)) ;
            else
                angles(i,j) = 90 + atand(yg(i,j)/xg(i,j)) ;
            end
        end
    end
end

fid = fopen('Angles.txt','w') ;
fprintf(fid,'%d %d %f %f %f\n',nrows,ncols,dx,x0y0(1),x0y0(2)) ;
for i = 1 : ncols
    for j = 1 : nrows 
       fprintf(fid,'%f ', (angles(i,j))*(pi/180)); % sz(i,j)) ; 
    end
    fprintf(fid,'\n') ; 
end
fclose(fid) ;

figure ;
pcolor(xg,yg,angles) ;
colorbar ;
shading interp
title('Elong')

