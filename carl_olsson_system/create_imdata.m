function create_imdata(settings)
save_path = settings.save_path;
merge_tracks = settings.merge_tracks;

load(strcat(save_path,'pairwise_matchings.mat'));

G = zeros(size(pairwiseEst));
for i = 1:size(pairwiseEst,1);
    for j = i+1:size(pairwiseEst,2);
        if ~isempty(pairwiseEst{i,j})
            G(i,j) = length(pairwiseEst{i,j}.ind1);
        end
    end
end
G = G+G';

%Start with a random image
S = false(size(G,1),1);
S(ceil(size(G,1)*rand),1) = true;

siftind = cell(size(imnames));

impoints.pointnr = 0;
impoints.points = cell(size(imnames));
if settings.storesift == 1,
    impoints.sift = cell(size(imnames));
end
impoints.index = cell(size(imnames));

while max(G(:)) ~= 0

    %Choose the pair with the most matches.
    Gtmp = repmat(S,[1 length(S)]).*G;
    [i,j] = find(Gtmp == max(Gtmp(:)));
    if i(1) > j(1);
        tmp = j(1);
        j = i(1);
        i = tmp;
        S(i) = true;
    else
        j = j(1);
        i = i(1);
        S(j) = true;
    end
    
    %FInd index to corresponding points
    ind1 = pairwiseEst{i,j}.ind1;
    ind2 = pairwiseEst{i,j}.ind2;
    %Remove non-unique points
    [uniqind2,I,J] = unique(ind2);
    ind2 = ind2(I);
    ind1 = ind1(I);
    [uniqind1,I,J] = unique(ind1);
    ind2 = ind2(I);
    ind1 = ind1(I);
        
    %Check which points have been added to tracks eariler.
    [members1,positions1] = ismember(ind1,siftind{i});
    [members2,positions2] = ismember(ind2,siftind{j});
    
    %Add match that have not been seen earlier.
    unseen = (~members1) & (~members2);
    %camera i
    impoints.index{i} = [impoints.index{i} impoints.pointnr+[1:sum(unseen)]];
    siftind{i} = [siftind{i} ind1(unseen)];
    impoints.points{i} = [impoints.points{i} SIFT{i}.locs(ind1(unseen),[2 1])'];
    if settings.storesift==1,
        impoints.sift{i} = [impoints.sift{i} SIFT{i}.desc(ind1(unseen),:)'];
    end
    %camera j
    impoints.index{j} = [impoints.index{j} impoints.pointnr+[1:sum(unseen)]];
    siftind{j} = [siftind{j} ind2(unseen)];
    impoints.points{j} = [impoints.points{j} SIFT{j}.locs(ind2(unseen),[2 1])'];
    if settings.storesift==1,
        impoints.sift{j} = [impoints.sift{j} SIFT{j}.desc(ind2(unseen),:)'];
    end
    %Update the numer of tracks.
    impoints.pointnr = impoints.pointnr + sum(unseen);
    
    
    %Add points that has been seen in camera i previously
    unseen_cam2 = members1 & (~members2);
    %Find corresponding track number
    newindex_cam2 = impoints.index{i}(positions1(unseen_cam2));
    %Find index that is already taken
    occupied_newind2 = ismember(newindex_cam2,impoints.index{j});
    new_siftind2 = ind2(unseen_cam2);
    new_siftind2 = new_siftind2(~occupied_newind2);
    
    impoints.index{j} = [impoints.index{j} newindex_cam2(~occupied_newind2)];
    impoints.points{j} = [impoints.points{j} SIFT{j}.locs(new_siftind2,[2 1])'];
    if settings.storesift==1,
        impoints.sift{j} = [impoints.sift{j} SIFT{j}.desc(new_siftind2,:)'];
    end
    siftind{j} = [siftind{j} new_siftind2];
           
    
    %Add points that has been seen in camera j previously
    unseen_cam1 = (~members1) & members2;
    %Find corresponding track number
    newind_cam1 = impoints.index{j}(positions2(unseen_cam1));
    %Find index that is already taken
    occupied_newind1 = ismember(newind_cam1,impoints.index{i});
    new_siftind1 = ind1(unseen_cam1);
    new_siftind1 = new_siftind1(~occupied_newind1);
    
    impoints.index{i} = [impoints.index{i} newind_cam1(~occupied_newind1)];
    impoints.points{i} = [impoints.points{i} SIFT{i}.locs(new_siftind1,[2 1])'];
    if settings.storesift==1,
        impoints.sift{i} = [impoints.sift{i} SIFT{i}.desc(new_siftind1,:)'];
    end
    siftind{i} = [siftind{i} new_siftind1];
    
    %Points that have been seen in both cameras previously
    seen = members1 & members2;
    if merge_tracks
        %Find the track index
        seen_ind1 = impoints.index{i}(positions1(seen));
        seen_ind2 = impoints.index{j}(positions2(seen));
        %Find those that don't have the same index
        merge1 = seen_ind1(seen_ind1 ~= seen_ind2);
        merge2 = seen_ind2(seen_ind1 ~= seen_ind2);
        
        %Check if they can be merged
        tracks1 = false(length(merge1),length(impoints.index));
        tracks2 = false(length(merge1),length(impoints.index));
        if size(tracks1,1) ~= 0
            for ii = 1:length(impoints.index)
                tracks1(:,ii) = ismember(merge1,impoints.index{ii});
                tracks2(:,ii) = ismember(merge2,impoints.index{ii});
            end
            non_conflict = (sum(tracks1 & tracks2,2) == 0);
            merge1 = merge1(non_conflict);
            merge2 = merge2(non_conflict);
            
            %Merge
            for ii = 1:length(impoints.index);
                [mem,pos] = ismember(merge2,impoints.index{ii});
                impoints.index{ii}(pos(mem)) = merge1(mem);
            end
            
            %Clean up the index
            index = 1:impoints.pointnr;
            newindex = cumsum(~ismember(index,merge2));
            for ii = 1:length(impoints.index);
                impoints.index{ii} = newindex(impoints.index{ii});
            end
            impoints.pointnr = impoints.pointnr - sum(non_conflict);
        end
    end
    [sum(unseen) sum(unseen_cam2) sum(unseen_cam1) sum(seen)]
    
    G(i,j) = 0;
    G(j,i) = 0;
    
    for iii = 1:length(impoints.index)
        if length(unique(impoints.index{iii})) ~= length(impoints.index{iii})
            disp('non unique index'); %Not good if this message is displayed
        end
    end
end

save(strcat(save_path,'impoints.mat'),'impoints','imnames');

if 0 %Debug code
    i = 1;
    j = 2;
    p1 = NaN*ones(2,impoints.pointnr);
    p1(:,impoints.index{i}) = impoints.points{i};
    p2 = NaN*ones(2,impoints.pointnr);
    p2(:,impoints.index{j}) = impoints.points{j};
    vis = isfinite(p1(1,:)) & isfinite(p2(1,:));
    figure(1);
    plot(imagedata(strcat(img_path,imnames(i).name),p1(:,vis)));
    figure(2);
    plot(imagedata(strcat(img_path,imnames(j).name),p2(:,vis)));
end

