function imnames = RANSAC_orientations(settings)
save_path = settings.save_path;
mincorrnr = settings.mincorrnr;
roterrtol = settings.roterrtol;

load(fullfile(save_path,'pairwise_geom.mat'));
load(fullfile(save_path,'impoints.mat'));

G = zeros(size(pairwise_geom));
for i = 1:size(pairwise_geom,1)
    for j = 1:size(pairwise_geom,1)
        if ~isempty(pairwise_geom{i,j})
            G(i,j) = length(pairwise_geom{i,j}.inliers);
            L(i,j) = -1;
            if G(i,j) < mincorrnr;
                G(i,j) = 0;
                L(i,j) = 0;
            end
        end
    end
end
G = (G+G');

maxinl = 0;
for jj = 1:100
    A = cell(1,size(pairwise_geom,1));
    S = zeros(size(pairwise_geom,1),1);
    %Randomly select first camera
    i = ceil(rand*length(S));
    S(i) = 1;
    A{i} = eye(3);
    campairs = [];
    while sum(S) < length(S);
        if 0
            %Choose the camera with most matches
            [i,j,W] = find(S*(1-S)'.*G);
            if isempty(W)
                break;
            end
            ind = find(max(W) == W);
            ind = ind(1);
            i = i(ind);
            j = j(ind);
        end
        if 1
            %Choose camera with many matches
            [i,j,W] = find(S*(1-S)'.*G);
            if isempty(W)
                break;
            end
            prob = rand;
            ind = find(cumsum(W.^5)./sum(W.^5) >= prob);
            ind = ind(1);
            i = i(ind);
            j = j(ind);
        end
        
        campairs = [campairs [i;j;W(ind)]];
        %Relative orientaiton (from epipolar geometry)
        if i < j
            P1 = pairwise_geom{i,j}.P{1};
            P2 = pairwise_geom{i,j}.P{2};
        else
            P1 = pairwise_geom{j,i}.P{2};
            P2 = pairwise_geom{j,i}.P{1};
        end
        R1 = inv(P1(:,1:3));
        R2 = A{i};
        A{j} = P2(:,1:3)*R1*R2;
        S(j) = 1;
    end
    
    res = zeros(size(pairwise_geom));
    weightres = zeros(size(pairwise_geom));
    inl = zeros(size(pairwise_geom));
    for i = 1:size(pairwise_geom,1);
        for j = i+1:size(pairwise_geom,2);
            if ~isempty(pairwise_geom{i,j}) & G(i,j) >= mincorrnr;
                P2 = pairwise_geom{i,j}.P{2};
                if ~isempty(A{i}) && ~isempty(A{j});
                    T = inv(A{i}(:,1:3));
                    A1 = A{j}*T;
                    
                    if 0 %Correct metric, but somthing wrong...
                        relrot = A1'*P2(:,1:3);
                        [~,~,V] = svd(relrot-eye(3));
                        rotaxis = V(:,end);
                        if abs(rotaxis(1)) > 0.1;
                            perpvec = [rotaxis(2); -rotaxis(1); 0];
                        else
                            perpvec = [0; rotaxis(3); -rotaxis(2)];
                        end
                        perpvec = perpvec./norm(perpvec);
                        res(i,j) = acos(perpvec'*relrot*perpvec);
                    else %Easy metric
                        res(i,j) = norm(A1-P2(:,1:3));
                    end
                    
                    weightres(i,j) = length(pairwise_geom{i,j}.inliers)*norm(A1-P2(:,1:3));
                    if res(i,j) < roterrtol
                        inl(i,j) = 1;
                    end
                end
            end
        end
    end
    if sum(sum(inl)) > sum(sum(maxinl));
        Amax = A;
        maxinl = inl;
        maxres = res;
        maxweightres = weightres;
        maxcampairs = campairs;
    end
    if settings.debug_match
        fprintf('RANSAC_ROT1, Iteration: %3d. Max Inliers: %5d.\n',jj,sum(sum(maxinl)));
    end
    %[jj sum(sum(maxinl))]
end

if settings.debug_match
    fprintf('RANSAC_ROT1, Max Res: %f.\n',sum(maxres(:)));
end
if settings.debug_match
    fprintf('RANSAC_ROT1, Max Inl: %4d.\n',sum(maxinl(:)));
end
%sum(maxres(:))
%sum(maxinl(:))
%keyboard;
A = Amax;
save(fullfile(save_path,'rotations_inliers.mat'),'maxinl');
save(fullfile(save_path,'rotations.mat'),'A');
