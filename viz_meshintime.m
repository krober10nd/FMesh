close all; clc; 

t = []; p = []; 

nscreen=5; maxiter=500 ; 

figure; 
for i = 5:5:maxiter
    n_strPadded = sprintf( '%04d', i) ; 
    Pfname = ['POINTS',n_strPadded,'.TXT'] ;
    Tfname = ['FACETS',n_strPadded,'.TXT'] ;
    
    p = dlmread(Pfname) ;
    t = dlmread(Tfname) ; 
    length(p)
    
    cla, triplot(t,p(:,1),p(:,2));
    axis([-1 1 -1 1])
    axis equal
    %axis([   -0.4888   -0.0931    0.3473    0.7430])
    title(['ITERATION=',num2str(i)]) ;
    pause(0.01)
end