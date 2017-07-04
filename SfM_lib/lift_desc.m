function kp = lift_desc(im);


imwrite(im,'C:/cygwin64/home/fredrik/tmp.tif');
system('C:\cygwin64\bin\bash.exe --login -c "./extractFeatures.sh tmp.tif"');
system('C:\cygwin64\bin\bash.exe --login -c "tar xvf results_tmp.tif.tar.gz"');

kp = readKPandDesc('C:/cygwin64/home/fredrik/results');





