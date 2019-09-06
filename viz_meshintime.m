t = []; p = []; 

figure; 
for i = 0 :19
    n_strPadded = sprintf( '%04d', i) ; 
    Pfname = ['POINTS',n_strPadded,'.TXT'] ;
    Tfname = ['FACETS',n_strPadded,'.TXT'] ;
    
    p = dlmread(Pfname) ;
    t = dlmread(Tfname) ;
    length(p)
    
    cla, triplot(t,p(:,1),p(:,2));
    axis([-1 1 -1 1])
    title(['ITERATION=',num2str(i)]) ;
    pause
end