function [maxinlier,P] = RANSAC_Essential(p2,p1,pixtol)

ind = find(isfinite(p1(1,:)) & isfinite(p2(1,:)));

mmax = [];
maxinlier = 0;
for ii = 1:50
    %slumpa 5 pkt & bestäm geometri.
    randind = randperm(length(ind));
    
    pp1 = p1(:,randind(1:5));
    pp2 = p2(:,randind(1:5));
    Evec = calibrated_fivepoint(pp2,pp1);

    for iiii = 1:size(Evec,2);
        E = reshape(Evec(:,iiii),3,3);
        [U,S,V] = svd(E);
        if det(U*V') < 0
            V = -V;
        end
        W = [0 1 0; -1 0 0 ; 0 0 1];
        t = U(:,3);
        R = U*W*V';
        P2 = [];       
        P{1} = [eye(3) zeros(3,1)];
        P{2} = [R t];
        
        UU = intsec2views_midpoint(P{1},P{2},pp1,pp2); %Triangulation using the midpoint method, not recommended but seems to work
       
        posDepth = sum([P{1}(3,:)*UU P{2}(3,:)*UU] > 0);
        P2 = P{2};
        
        P{2} = [R -t];
        UU = intsec2views_midpoint(P{1},P{2},pp1,pp2);
        if sum([P{1}(3,:)*UU P{2}(3,:)*UU] > 0) > posDepth
            P2 = P{2}; 
            posDepth = sum([P{1}(3,:)*UU P{2}(3,:)*UU] > 0);
        end
        
        R = U*W'*V';
        P{2} = [R t];
        UU = intsec2views_midpoint(P{1},P{2},pp1,pp2);
        if sum([P{1}(3,:)*UU P{2}(3,:)*UU] > 0) > posDepth
            P2 = P{2}; 
            posDepth = sum([P{1}(3,:)*UU P{2}(3,:)*UU] > 0);
        end
        
        P{2} = [R -t];
        UU = intsec2views_midpoint(P{1},P{2},pp1,pp2);
        if sum([P{1}(3,:)*UU P{2}(3,:)*UU] > 0) > posDepth
            P2 = P{2}; 
            posDepth = sum([P{1}(3,:)*UU P{2}(3,:)*UU] > 0);
        end
        
        if det(P{2}(:,1:3)) < 0
            disp('Ej Rotation');
        end
        
        if ~isempty(P2)
            %Triangulate all points
            P{1} = [eye(3) zeros(3,1)];
            P{2} = P2;
            U = intsec2views_midpoint(P{1},P{2},p1,p2);
            %check number of inliers
            err = sqrt(sum((p1-pflat(P{1}*U)).^2)+sum((p2-pflat(P{2}*U)).^2));
            mindepth = min(P{1}(3,:)*U,P{2}(3,:)*U);
            inlier = err < pixtol & mindepth > 0;
           
            if sum(maxinlier) < sum(inlier);
                Pmax = P;                
                maxinlier = inlier;
            end
        end
    end
end
%If we find more than 5 inliers
if sum(maxinlier) > 5    
    P = Pmax;
    U = intsec2views(P{1},P{2},p1,p2); %Optimal two view triangulation
    
    err = sqrt(sum((p1-pflat(P{1}*U)).^2)+sum((p2-pflat(P{2}*U)).^2));
    mindepth = min(P{1}(3,:)*U,P{2}(3,:)*U);
    maxinlier = err < pixtol & mindepth > 0; 
   
    u2.pointnr = size(p1(:,maxinlier),2);
    u2.points{1} = p1(:,maxinlier);
    u2.index{1} = 1:u2.pointnr;
    u2.points{2} = p2(:,maxinlier);
    u2.index{2} = 1:u2.pointnr;
    
    [U,P] = modbundle_sparse(U(:,maxinlier),P,u2,20,0.001);
end


function y = pflat(x)
y = x./repmat(x(end,:),[size(x,1) 1]);

function y = pextend(x)
y = [x; ones(1,size(x,2))];

