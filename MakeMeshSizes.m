clearvars; close all; clc; 

dx =0.1 ; 
[xg,yg]=meshgrid(-5:dx:5,-5:dx:5);

sz = xg.^2 + yg.^2 ;
sz(sz<0.01) = 0.01 ; 
sz(sz>0.15) = 0.15 ; 
figure; pcolor(xg,yg,sz); shading interp
min(sz(:))
max(sz(:))

%% HEADER LOOKS LIKE THE BELOW LINE 
% 62   READ(1,*) SzFx%Ni,SzFx%Nj,SzFx%delta,SzFx%x0y0(1),SzFx%x0y0(2)
fid = fopen('MeshSizes.txt','w') ; 
fprintf(fid,'%d %d %f %f %f\n',101,101,0.1,-5.0,-5.0) ; 
for i = 1 : 101
    for j = 1 : 101 
       fprintf(fid,'%f ',sz(i,j)) ; 
    end
    fprintf(fid,'\n') ; 
end
fclose(fid) ;

fid = fopen('ElongSizes.txt','w') ; 
fprintf(fid,'%d %d %f %f %f\n',101,101,0.1,-5.0,-5.0) ; 
for i = 1 : 101
    for j = 1 : 101 
       fprintf(fid,'%f ', 10.0); 
    end
    fprintf(fid,'\n') ; 
end
fclose(fid) ;

fid = fopen('Angles.txt','w') ;
fprintf(fid,'%d %d %f %f %f\n',101,101,0.1,-5.0,-5.0) ; 
for i = 1 : 101
    for j = 1 : 101 
       fprintf(fid,'%f ',45.0) ; 
    end
    fprintf(fid,'\n') ; 
end
fclose(fid) ;
