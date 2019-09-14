close all; clc; 

t = []; p = []; 

nscreen=10 ; maxiter=100 ; 

figure; 
for i = nscreen:nscreen:maxiter
    n_strPadded = sprintf( '%04d', i) ; 
    Pfname = ['POINTS',n_strPadded,'.TXT'] ;
    Tfname = ['FACETS',n_strPadded,'.TXT'] ;
    
    p = dlmread(Pfname) ;
    t = dlmread(Tfname) ;
    length(p)
    
    cla, triplot(t,p(:,1),p(:,2));
    title(['ITERATION=',num2str(i)]) ;
    axis([-1 1 -1 1])
    pause(0.01)
end