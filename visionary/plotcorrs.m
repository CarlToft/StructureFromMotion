function fnbr=plotmatches(im1,im2,ps1,ps2,pmode,plotargs)
if nargin<6
    plotargs='-';
end
if nargin<5
    pmode=0;
end


[m1,n1,~]=size(im1);
[m2,n2,~]=size(im1);

imtot=[im1 im2];
fnbr=figure();

if pmode==0
    imagesc(imtot);hold on;
    plot([ps1(1,:);ps2(1,:)+n1],[ps1(2,:);ps2(2,:)],plotargs);
else
    for ii=1:pmode:length(ps1)-pmode
        imagesc(imtot);hold on;
        plot([ps1(1,ii:ii+pmode-1);ps2(1,ii:ii-1+pmode)+n1],...
            [ps1(2,ii:ii-1+pmode);ps2(2,ii:ii-1+pmode)],plotargs);
        hold off;
        pause
    end
end
hold off;