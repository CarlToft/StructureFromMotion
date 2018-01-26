function [tribottle, Ubottle] = annotate_get_bottle_model()
    alpha = 0:0.4:2*pi;rr = 0.13;height=0.8;nn=length(alpha);
    Ubottle=[rr*cos(alpha),rr*cos(alpha),0;rr*sin(alpha),rr*sin(alpha),0;zeros(size(alpha)),height*ones(size(alpha)),1];
    tribottle = [1:nn,      nn+1:2*nn, 2*nn+1*ones(1,nn) ; ...
                 [2:nn,1],  [nn+2:2*nn,nn+1]  , nn+1:2*nn, ;...
                 nn+1:2*nn, [2:nn,1]       , [nn+2:2*nn,nn+1]];
end
