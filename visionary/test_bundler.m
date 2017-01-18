

%test-script
bibdir='C:\bundler\bundler-v0.3-binary\examples\ET';

[str,mot,Kmot,f,pp,kk, imseq0, imundist]=bundler_import(bibdir);

figure(1);plot(imseq0{1},'numbered')
figure(2);plot(imseq0{2},'numbered')

figure(3);
plot(mot,str);axis equal;


figure(4);
rmspoints(str,mot,imundist);
figure(5);
reproject(str,mot,imundist);

% Do calibrated bundle adjustment
%[ss,mm]=bundleplc(str,mot, imundist, {'autocalib=11111'});
