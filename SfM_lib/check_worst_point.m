function [the_index,neg_depths]=check_worst_point(settings,P,U,u, ignoreindex)

if nargin<5,
    ignoreindex=[];
end

KK = settings.KK;
kc = settings.kc;
fc = KK([1 5]);
cc = KK(1:2,3);
alpha_c = KK(1,2)/fc(1);

the_rms = -1;

neg_depths = [];

%    'hej'
for i=1:length(P);
    ind=u.index{i};
    [ind,ind2] = setdiff(ind,[]);
    [ind,ind2] = setdiff(ind,ignoreindex);
    
    tmp = KK*pextend(apply_distortion(u.points{i}(1:2,ind2),kc));
    unflat = P{i}*U(:,ind);
    depths = unflat(3,:);
    neg_depths = union(neg_depths,ind(find(depths<0)));
    
    PU = pflat(unflat);
    PU = KK*pextend(apply_distortion(PU(1:2,:),kc));
    
    rms = sum((tmp-PU).^2);
    if max(rms)>the_rms,
        [aa,bb]=max(rms);
        the_rms = aa;
        the_index = ind(bb);
        the_i = i;
    end
end

[the_index, the_i]
plot_1point(settings,P,U,u,the_index);

neg_depths


function y = pflat(x)
y = x./repmat(x(end,:),[size(x,1) 1]);
