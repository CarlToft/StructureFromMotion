function rel_or_5points(settings)
save_path = settings.save_path;
RANSAC_pixtol = settings.RANSAC_pixtol;
KK = settings.KK;
kc = settings.kc;
camera_graph = settings.camera_graph;



load(fullfile(save_path,'impoints.mat'));

%beräkna alla parvisa geometrier

pairwise_geom = cell(length(imnames),length(imnames));
for i = 1:length(imnames);
    for j = i+1:length(imnames);
        
        if isempty(camera_graph) || camera_graph(i,j) == 1,
            
            
            %[i j length(imnames)]
            p1 = NaN*ones(2,impoints.pointnr);
            p2 = NaN*ones(2,impoints.pointnr);
            p1(:,impoints.index{i}) = impoints.points{i};
            p2(:,impoints.index{j}) = impoints.points{j};
            
            vis = isfinite(p1(1,:)) & isfinite(p2(1,:));
            
            fc = KK([1 5]);
            cc = KK(1:2,3);
            alpha_c = KK(1,2)/fc(1);
            pp1 = pextend(normalize(p1([1,2],vis),fc,cc,kc,alpha_c));
            pp2 = pextend(normalize(p2([1,2],vis),fc,cc,kc,alpha_c));
            
            %        pp1 = pflat(inv(KK)*pextend(p1(:,vis)));
            %        pp2 = pflat(inv(KK)*pextend(p2(:,vis)));
            
            if sum(vis) > settings.mininlnr,
                %5-punkt RANSAC
                [maxinliers,Pmax] = RANSAC_Essential(pp2,pp1,RANSAC_pixtol/KK(1,1));
                if ~isempty(Pmax);
                    if settings.debug_match
                        fprintf('FivePointRansac: (%3d,%3d) of %3d. Inliers %5d of %5d.\n',i,j,length(imnames),sum(maxinliers),size(pp2,2));
                    end
                    %[i j length(imnames)]
                    %sum(maxinliers)
                    vis = find(vis);
                    pairwise_geom{i,j}.P = Pmax;
                    pairwise_geom{i,j}.inliers = vis(maxinliers);
                end
            end
        end
    end
end

save(fullfile(save_path,'pairwise_geom.mat'),'pairwise_geom');

%keyboard;

if 0
    load(fullfile(save_path,'pairwise_geom.mat'));
    
    i = 167;
    j = 168;
    P = cell(1,2);
    if ~isempty(pairwise_geom{i,j});
        m = motion(pairwise_geom{i,j}.P);
        inl = pairwise_geom{i,j}.inliers;
        p1 = NaN(3,impoints.pointnr);
        p1(:,impoints.index{i}) = pextend(impoints.points{i});
        imseq2{1} = imagedata([],inv(KK)*p1(:,inl));
        p2 = NaN(3,impoints.pointnr);
        p2(:,impoints.index{j}) = pextend(impoints.points{j});
        imseq2{2} = imagedata([],inv(KK)*p2(:,inl));
        curr=pwd;
        cd ../visionary
        s = intsec2views(m,imseq2);
        cd (curr);
        figure(1);
        plot(m,s);
        figure(2);
        plot(imagedata(strcat(settings.img_path,imnames(i).name), p1(:,inl)));
        plot(imagedata([],KK*getcameras(m,1)*getpoints(s)),'numberedro');
        figure(3);
        plot(imagedata(strcat(settings.img_path,imnames(j).name), p2(:,inl)));
        plot(imagedata([],KK*getcameras(m,2)*getpoints(s)),'numberedro');
    else
        disp('Empty!');
    end
end

function y = pflat(x)
y = x./repmat(x(end,:),[size(x,1) 1]);

function y = pextend(x)
y = [x; ones(1,size(x,2))];

