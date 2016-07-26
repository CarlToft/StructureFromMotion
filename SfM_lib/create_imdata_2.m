function create_imdata_2(settings)
mininlnr = settings.mininlnr;
save_path = settings.save_path;
merge_tracks = settings.merge_tracks;
KK = settings.KK;
kc = settings.kc;

load(strcat(save_path,'impoints2.mat'));
load(strcat(save_path,'pairwise_geom2.mat'));
load(strcat(save_path,'rotations_inliers2.mat'));


G = zeros(size(pairwise_geom));
for i = 1:size(pairwise_geom,1)
    for j = 1:size(pairwise_geom,1)
        if ~isempty(pairwise_geom{i,j})
            if length(pairwise_geom{i,j}.inliers) > mininlnr &&  maxinl(i,j) == 1;
                G(i,j) = length(pairwise_geom{i,j}.inliers);
            end
        end
    end
end
G = G+G';

%Starta med slumpässig bild
S = false(size(G,1),1);
S(ceil(size(G,1)*rand),1) = true;

pointind = cell(size(imnames));
u.pointnr = 0;
u.points = cell(size(imnames));
if settings.storesift == 1,
    u.sift = cell(size(imnames));
end
u.index = cell(size(imnames));

Gtmp = repmat(S,[1 length(S)]).*G;
while max(Gtmp(:)) ~= 0

    %Välj par med mest matchningar där minst en finns i S
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
    
    %Ta bara med inliers från parvisa_geometrin
    ind1 = pairwise_geom{i,j}.inliers;
    ind2 = pairwise_geom{i,j}.inliers;

    %Plocka bort icke unika siftar.
    [uniqind2,I,J] = unique(ind2);
    ind2 = ind2(I);
    ind1 = ind1(I);
    
    %Kolla vilka som blivit matchade tidigare
    [members1,positions1] = ismember(ind1,pointind{i});
    [members2,positions2] = ismember(ind2,pointind{j});
    
    p1 = NaN*ones(2,impoints.pointnr);
    p1(:,impoints.index{i}) = impoints.points{i};
    p2 = NaN*ones(2,impoints.pointnr);
    p2(:,impoints.index{j}) = impoints.points{j};
    if settings.storesift == 1,
        p1sift = zeros(128,impoints.pointnr,'uint8');
        p1sift(:,impoints.index{i}) = impoints.sift{i};
        p2sift = zeros(128,impoints.pointnr,'uint8');
        p2sift(:,impoints.index{j}) = impoints.sift{j};
    end


    %Lägg till punkter som inte detekterats i någon av kamrorna (i paret) tidigare
    unseen = (~members1) & (~members2);
    %kamera i
    u.index{i} = [u.index{i} u.pointnr+[1:sum(unseen)]];
    pointind{i} = [pointind{i} ind1(unseen)];
    u.points{i} = [u.points{i} p1(:,ind1(unseen))];
    if settings.storesift==1,
        u.sift{i} = [u.sift{i} p1sift(:,ind1(unseen))];
    end
    %kamera j
    u.index{j} = [u.index{j} u.pointnr+[1:sum(unseen)]];
    pointind{j} = [pointind{j} ind2(unseen)];
    u.points{j} = [u.points{j} p2(:,ind2(unseen))];
    if settings.storesift==1,
        u.sift{j} = [u.sift{j} p2sift(:,ind2(unseen))];
    end
    %uppdatera antalet 3D-punkter
    u.pointnr = u.pointnr + sum(unseen);
    
    
    %Punkter som bara detekterats i kamera i tidigare
    unseen_cam2 = members1 & (~members2);
    %hitta motsvarande index för 3D-punkter
    newindex_cam2 = u.index{i}(positions1(unseen_cam2));
    %hitta index som redan är upptagna
    occupied_newind2 = ismember(newindex_cam2,u.index{j});
    new_pointind2 = ind2(unseen_cam2);
    new_pointind2 = new_pointind2(~occupied_newind2);
    
    u.index{j} = [u.index{j} newindex_cam2(~occupied_newind2)];
    u.points{j} = [u.points{j} p2(:,new_pointind2)];
    if settings.storesift==1,
        u.sift{j} = [u.sift{j} p2sift(:,new_pointind2)];
    end
    pointind{j} = [pointind{j} new_pointind2];
           
    
    %Punkter som bara detekterats i kamera j tidigare
    unseen_cam1 = (~members1) & members2;
    %hitta motsvarande index för 3D-punkter
    newind_cam1 = u.index{j}(positions2(unseen_cam1));
    %hitta index som redan är upptagna
    occupied_newind1 = ismember(newind_cam1,u.index{i});
    new_pointind1 = ind1(unseen_cam1);
    new_pointind1 = new_pointind1(~occupied_newind1);
    
    u.index{i} = [u.index{i} newind_cam1(~occupied_newind1)];
    u.points{i} = [u.points{i} p1(:,new_pointind1)];
    if settings.storesift==1,
        u.sift{i} = [u.sift{i} p1sift(:,new_pointind1)];
    end
    pointind{i} = [pointind{i} new_pointind1];
    
    %Punkter som detekterats i båda kamrorna
    seen = members1 & members2;
    if merge_tracks
        %Hitta indexen
        seen_ind1 = u.index{i}(positions1(seen));
        seen_ind2 = u.index{j}(positions2(seen));
        %Hitta dom som inte har samma index
        merge1 = seen_ind1(seen_ind1 ~= seen_ind2);
        merge2 = seen_ind2(seen_ind1 ~= seen_ind2);
        
        %Kolla om det går att slå ihop
        tracks1 = false(length(merge1),length(u.index));
        tracks2 = false(length(merge1),length(u.index));
        if size(tracks1,1) ~= 0
            for ii = 1:length(u.index)
                tracks1(:,ii) = ismember(merge1,u.index{ii});
                tracks2(:,ii) = ismember(merge2,u.index{ii});
            end
            non_conflict = (sum(tracks1 & tracks2,2) == 0);
            merge1 = merge1(non_conflict);
            merge2 = merge2(non_conflict);
            
            %Slå ihop
            for ii = 1:length(u.index);
                %Byt ut index merge2 mot merge1
                [mem,pos] = ismember(merge2,u.index{ii});
                u.index{ii}(pos(mem)) = merge1(mem);
            end
            
            %Rensa upp bland indexen
            index = 1:u.pointnr;
            newindex = cumsum(~ismember(index,merge2));
            for ii = 1:length(u.index);
                u.index{ii} = newindex(u.index{ii});
            end
            u.pointnr = u.pointnr - sum(non_conflict);
        end
    end
    [sum(unseen) sum(unseen_cam2) sum(unseen_cam1) sum(seen)]
    
    G(i,j) = 0;
    G(j,i) = 0;
    
    for iii = 1:length(u.index)
        if length(unique(u.index{iii})) ~= length(u.index{iii})
            disp('non unique index');
        end
    end
    Gtmp = repmat(S,[1 length(S)]).*G;
end

fc = KK([1 5]);
cc = KK(1:2,3);
alpha_c = KK(1,2)/fc(1);

%kalibrera och plocka bort
for i=1:length(u.points);
    if ~isempty(u.index{i})
        u.points{i} = pextend(normalize(u.points{i}([1,2],:),fc,cc,kc,alpha_c));
    end
end
u.points = u.points(S);
if settings.storesift==1
    u.sift = u.sift(S);
end
u.index = u.index(S);
imnames = imnames(S);
save(strcat(save_path,'impoints3.mat'),'u','imnames');
load(strcat(save_path,'rotations2.mat'),'A');
A = A(S);
save(strcat(save_path,'rotations2.mat'),'A');


if 0
    i = ceil(rand*length(imnames));
    j = ceil(rand*length(imnames));
    while isempty(pairwise_geom{i,j})
        i = ceil(rand*length(imnames));
        
        j = ceil(rand*length(imnames));
    end
    p1 = NaN*ones(3,impoints.pointnr);
    p2 = NaN*ones(3,impoints.pointnr);
    p1(:,u.index{i}) = u.points{i};
    p2(:,u.index{j}) = u.points{j};
    
    vis = isfinite(p1(1,:)) & isfinite(p2(1,:));
    figure(1);
    plot(imagedata(strcat(img_path,imnames(i).name),KK*p1(:,vis)));
    figure(2);
    plot(imagedata(strcat(img_path,imnames(j).name),KK*p2(:,vis)));
    
    ind1 = pairwise_geom{i,j}.inliers;
    ind2 = pairwise_geom{i,j}.inliers;
    p1 = NaN*ones(2,impoints.pointnr);
    p1(:,impoints.index{i}) = impoints.points{i};
    p2 = NaN*ones(2,impoints.pointnr);
    p2(:,impoints.index{j}) = impoints.points{j};
    figure(3);
    plot(imagedata(strcat(img_path,imnames(i).name),pextend(p1(:,ind1))));
    figure(4);
    plot(imagedata(strcat(img_path,imnames(j).name),pextend(p2(:,ind2))));
end

if 0
    outliertracks = false(length(imnames),u.pointnr);
    for i = 1:length(P);
        uu = NaN*ones(3,u.pointnr);
        uu(:,u.index{i}) = u.points{i};
        vis = find(isfinite(uu(1,:)));
        res = length(vis);
        outl = s(1:res) > 1e-5;
        outliertracks(i,vis) = outl';
        s = s(res+1:end);
    end
    
    
    track = ceil(rand*u.pointnr);
    while sum(outliertracks(:,track)) == 0
        track = ceil(rand*u.pointnr);
    end
    for i = 1:length(imnames);
        [mem,pos]=ismember(track,u.index{i});
        if mem
            figure(1);
            if outliertracks(i,track) == 1
                plot(imagedata(strcat(img_path,imnames(i).name),KK*u.points{i}(:,pos)),'r*');
            else
                plot(imagedata(strcat(img_path,imnames(i).name),KK*u.points{i}(:,pos)),'*');
            end
            p = KK*P{i}*U(:,track);
            hold on;
            plot(imagedata([],p),'go');
            hold off;
            pause
        end
    end
end

function y = pflat(x)
y = x./repmat(x(end,:),[size(x,1) 1]);

function y = pextend(x)
y = [x; ones(1,size(x,2))];
