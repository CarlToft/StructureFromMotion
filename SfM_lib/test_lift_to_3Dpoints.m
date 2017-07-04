
%load ../data/oxford/result_shortseq1-4.mat
load ../data/oxford/result_seq1-6_block1.mat

%shortseq1
result_name = 'shortseq1';
models = 20;nbr_images = 892;

result_name = 'seq1_block1';
models = 3;nbr_images = 149;%seq1-6_block1

KK = settings.KK;
kc = settings.kc;
fc = KK([1 5]);
cc = KK(1:2,3);
alpha_c = KK(1,2)/fc(1);
pixtol = 5;

settings.LIFT = 1;
settings.LIFT_thr = 5;
settings.epipoledistance = 20;
settings.storesift=1;

minvis = 3; %minimum number of views  
settings.visviews = minvis;
seqlength = 10;

seq_main = 1;


start_index = [round(1:nbr_images/models:nbr_images),nbr_images];

for run=1:models,
    run
    index = [start_index(run):start_index(run+1)];
    %add overlap
    if run==models,
%        index = [index,1:5]; %only for circular
    else
        index = [index,start_index(run+1)+[1:5]];
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

seqlist = {};imname = settings.imnames(1).name(1:19);
nn = 1;
for ii=2:length(settings.imnames),
    if strcmp(imname,settings.imnames(ii).name(1:19))==0,
        imname = settings.imnames(ii).name(1:19);
        seqlist{end+1}=nn;
        nn = ii;
    else
        nn(end+1)=ii;
    end
end
seqlist{end+1}=nn;

%fix camera graph if images have been removed
nbr = size(settings.imnames,2);
%index = seqlist{seq_main};
tmp = zeros(nbr,nbr);

for ii=1:length(index)-1,
    index2 = index(ii+1:min(ii+seqlength-1,length(index)));
    tmp(index(ii),index2)=1;
end
tmp=max(tmp,tmp');
settings.camera_graph = tmp;

seq=sort(index);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute [str,mot,imseqnew]

pairwise_matching(settings, seq);
create_imdata(settings);

save_path = settings.save_path;
load(strcat(save_path,'pairwise_matchings.mat'), 'SIFT', 'pairwiseEst', 'imnames');
load(strcat(save_path,'impoints.mat'),'impoints','imnames');

mot = motion(P_uncalib(seq));
imseq = cell(1,length(seq));
imseqnew = cell(1,length(seq));
lift_descriptors = cell(1,length(seq));
clear u;
u.pointnr = 0;

for ii=1:length(seq),
    x_norm = normalize(impoints.points{ii},fc,cc,kc,alpha_c);
    tmp=NaN*ones(2,impoints.pointnr);
    tmp(:,impoints.index{ii})=x_norm;
    imseq{ii}=imagedata([],tmp);
    imseqnew{ii} = imagedata;
    lift_descriptors{ii}=zeros(128,0);
    u.points{ii} = zeros(3,0);
    u.index{ii} = zeros(1,0);
end
str = structure;
for ii=1:impoints.pointnr,
    tmpimseq={};
    tmpmot = motion;
    vis = 0;
    for jj=1:length(seq);
        if ~isnan(getpoints(imseq{jj},ii)),
            vis = vis +1;
        end
    end
    if vis>=minvis; %DEMAND VISIBLE IN ALL?
        xlist=zeros(2,0);Plist={};
        for jj=1:length(seq);
            if ~isnan(getpoints(imseq{jj},ii)),
                tmp=getpoints(imseq{jj},ii);
                tmpimseq{end+1}=imagedata([],tmp);
                xlist(:,end+1)=tmp(1:2,1);
                Plist{end+1}=getcameras(mot,jj);
                tmpmot = addcameras(tmpmot,getcameras(mot,jj));
            end
        end
        tmpstr=intsecpoints(tmpmot,tmpimseq);
        tmpstr = bundleplc(tmpstr,tmpmot,tmpimseq,{'structure','iteration=10','outputoff'});
        rms=rmspoints(tmpstr,tmpmot,tmpimseq);

        if rms>pixtol/KK(1,1),
            Utmp=linf_triangulation(xlist,Plist,20,0.01);
            tmpstr = bundleplc(structure(Utmp),tmpmot,tmpimseq,{'structure','iteration=10','outputoff'});
            rms = rmspoints(tmpstr,tmpmot,tmpimseq);
        end
        
        if rms<pixtol/KK(1,1),
            for jj=1:length(seq);
                pt = getpoints(imseq{jj},ii);
                imseqnew{jj}=addpoints(imseqnew{jj},pt);
                if ~isnan(pt(1)),
                    ind = find(impoints.index{jj}==ii);
                    lift_descriptors{jj}(:,end+1)=impoints.sift{jj}(:,ind);
                    u.points{jj}(:,end+1)=pt;
                    u.index{jj}(:,end+1)=size(imseqnew{jj},1);
                end
            end
            str = str + tmpstr;
        end
    end
end
U = pflat(str);
u.pointnr = size(U,2);
P = getcameras(mot);
disp(['Points: ',num2str(size(str,1)),' of ',num2str(impoints.pointnr)]);

eval(['save ../data/oxford/lift_models/',result_name,'_model',num2str(run),'.mat settings str mot imseqnew u U P lift_descriptors']);
% Compute [str,mot,imseqnew]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%
end %run
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    


