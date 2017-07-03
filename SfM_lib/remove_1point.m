function [U,u]=remove_1point(P,U,u, theUindex)

theUindex = sort(theUindex);

for kk=length(theUindex):-1:1,
    Uindex = theUindex(kk);

    UU = [1:Uindex-1,0,Uindex:size(U,2)-1];
    U(:,Uindex)=[];
    u.pointnr = u.pointnr-1;

    for i=1:length(P);
        ind=find(u.index{i}==Uindex);
        if ~isempty(ind),
            u.points{i}(:,ind)=[];
            u.index{i}(ind)=[];
        end
        u.index{i}=UU(u.index{i});
    end
end

