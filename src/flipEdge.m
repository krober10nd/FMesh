function [t,t2t,t2n] = flipEdge(p,t,Me)
% Flips edges of an arbitrary triangulation to achieve Delaunay triangulation.
% Uses the notion of a metric tensor 2x2 (Me) to create anisotropy. 
% Me can be calculated using the the CalcIdealMetric function
% 
% kjr, 2019
% 
% Requies mkt2t from DistMesh2D package to determine the node not shared in
% each triangle in a patch of two elements that are connected. 
%
%
% Shared edge in triangle 1 is edge index points
%shared1 = newt(1,[tix11,tix12]) ;
%shared1V = p(shared1(2),:)' - p(shared1(1),:)'  ;

% shared edge in triangle 2 is edge index points
%shared2 = newt(2,[tix21,tix22]) ;
%shared2V = p(shared2(2),:)' - p(shared2(1),:)';

nt = length(t) ;

if nargin < 4
    [t2t,t2n]=mkt2t(t);
end

mod3x1 =[2,3,1];
mod3x2 =[3,1,2];
mod3x3 =[1,2,3] ;

for t1 = 1 : nt % for each tri
    for n1 = 1 : 3 % for each conn tri
        t2 = t2t(t1,n1) ; % connected triangle (shares edge)
        if t2 >= 1
            n2=t2n(t1,n1); % local node not shared by the two triangles that share an edges
            
            % once we know the node not shared by the two triangles, we can
            % figure out the shared edges easily
            tix11=mod3x1(n1);
            tix12=mod3x2(n1);
            tix13=mod3x3(n1);
            
            tix21=mod3x1(n2);
            tix22=mod3x2(n2);
            tix23=mod3x3(n2);
            
            % two connected triangles that share an edge
            newt=[t(t1,1) t(t1,2) t(t1,3);
                  t(t2,1) t(t2,2) t(t2,3) ];
            
            tpmid = squeeze(mean(reshape(p(newt,:),[],3,2),2));      % Compute centroid of patch
            tpmid = mean(tpmid,1) ; 
                        
            edgeBV = p(newt(1,tix13),:)' - p(newt(1,tix11),:)' ;
            
            edgeDV = p(newt(2,tix23),:)' - p(newt(2,tix21),:)' ;
            
            edgeCV = p(newt(1,tix13),:)' - p(newt(2,tix21),:)' ;
            
            edgeAV = p(newt(2,tix23),:)' - p(newt(1,tix11),:)' ;
            
            % cross-->(v1.X*v2.Y) - (v1.Y*v2.X);
            cp1 = ( (edgeAV(1)*edgeBV(2)) - (edgeAV(2)*edgeBV(1)) ) ;
            cp2 = ( (edgeCV(1)*edgeDV(2)) - (edgeCV(2)*edgeDV(1)) ) ;
            
            % del. criterion
            flip = cp1*(edgeCV'*Me*edgeDV)+ ...
                   cp2*(edgeAV'*Me*edgeBV) > 0 ;
            
            if flip
                % swap edge
                newt(1,tix12)=newt(2,n2);
                newt(2,tix22)=newt(1,n1);
                
                % Insert new triangles (with flipped edges)
                t(t1,:) = newt(1,:) ;
                t(t2,:) = newt(2,:) ;
                
                % Update t2t and t2n
                nbt=t2t(t2,tix21);
                nbn=t2n(t2,tix21);
                
                t2t(t1,n1)=nbt;
                t2n(t1,n1)=nbn;
                
                if nbt>=1
                    t2t(nbt,nbn)=t1;
                    t2n(nbt,nbn)=n1;
                end
                
                nbt=t2t(t1,tix11) ;
                nbn=t2n(t1,tix11) ;
                t2t(t2,n2)=nbt;
                t2n(t2,n2)=nbn;
                
                if nbt>=1
                    t2t(nbt,nbn)=t2;
                    t2n(nbt,nbn)=n2;
                end
                
                t2t(t1,tix11)=t2;
                t2n(t1,tix11)=tix21;
                t2t(t2,tix21)=t1;
                t2n(t2,tix21)=tix11;
            end % if flip
        end % if connected tria
    end % for each connected tria
end % for each tria
end

