



if 0,
    load C:\Users\fredrik\Documents\MATLAB\tunnels\experiment3_data\data6\result_data

    u0 = u;
    P0 = P;
    U0 = U;
    
    %add a couple of manual correspondences
    ii=1;jj=101;
    img_path = settings.img_path;
    imnames = settings.imnames;
    
    figure(1);
    im = imread(strcat(img_path,imnames(ii).name));
    imagesc(im); hold on;
    figure(2);
    im = imread(strcat(img_path,imnames(jj).name));
    imagesc(im); hold on;
    
   x1 = 1000*[ 
   0.474693548387097   0.115232824427481;...
   0.517274193548387   0.122103053435115;...
   0.536629032258064   0.122103053435115;...
   0.575338709677419   0.124851145038168;...
   0.617919354838710   0.128973282442748;...
   0.350822580645161   0.895690839694657;...
   0.537919354838710   0.734927480916031;...
   1.241145161290323   0.600270992366412;...
   1.541790322580645   1.011110687022901]';
   x2 = 1000*[
   0.795983870967742   0.189431297709924;...
   0.846306451612903   0.199049618320611;...
   0.864370967741936   0.201797709923664;...
   0.905661290322581   0.211416030534351;...
   0.945661290322581   0.218286259541985;...
   1.254048387096774   0.956148854961832;...
   1.060500000000000   0.800881679389313;...
   1.500500000000000   0.648362595419848;...
   1.828241935483871   0.908057251908397]';

   KK = settings.KK;
   kc = settings.kc;
   fc = KK([1 5]);
   cc = KK(1:2,3);
   alpha_c = KK(1,2)/fc(1);
   p1 = pextend(normalize(x1,fc,cc,kc,alpha_c));
   p2 = pextend(normalize(x2,fc,cc,kc,alpha_c));
   U = intsec2views(P{ii},P{jj},p1,p2);
   
   %reproject(structure(U),motion(P([ii,jj])),{imagedata([],p1),imagedata([],p2)},[],'numbered');
   
   
   %add to u0 and U0
   nbr0 = size(U0,2);
   nbr = size(U,2);
   U0 = [U0,U];
   u0.pointnr = nbr0 + nbr;
   u0.points{ii}=[u0.points{ii},p1];
   u0.index{ii}=[u0.index{ii},nbr0+1:nbr0+nbr];
   u0.points{jj}=[u0.points{jj},p2];
   u0.index{jj}=[u0.index{jj},nbr0+1:nbr0+nbr];
   
   [U0,P0,lambda] = modbundle_sparse(U0,P0,u0,20,0.01);
   

    %add a couple of manual correspondences
    ii=95;jj=200;
    img_path = settings.img_path;
    imnames = settings.imnames;
    
    figure(1);
    im = imread(strcat(img_path,imnames(ii).name));
    imagesc(im); hold on;
    figure(2);
    im = imread(strcat(img_path,imnames(jj).name));
    imagesc(im); hold on;
    
   x1 = 1000*[
   0.958564516129032   0.509583969465649;...
   0.986951612903226   0.495843511450382;...
   1.026951612903226   0.479354961832061;...
   1.097919354838710   0.453248091603054;...
   1.234693548387097   0.394164122137405;...
   1.563725806451613   0.270500000000000;...
   0.345661290322581   0.315843511450382;...
   0.492758064516129   0.396912213740458;...
   0.570177419354839   0.439507633587786;...
   0.620500000000000   0.464240458015267;...
   0.654048387096774   0.482103053435115]';
   x2 = 1000*[
   1.055338709677420   0.528820610687023;...
   1.092758064516129   0.524698473282443;...
   1.154693548387097   0.515080152671756;...
   1.248887096774194   0.501339694656489;...
   1.402435483870968   0.479354961832061;...
   1.683725806451613   0.438133587786260;...
   0.688887096774193   0.405156488549619;...
   0.723725806451613   0.450500000000000;...
   0.744370967741936   0.477980916030534;...
   0.759854838709678   0.497217557251908;...
   0.768887096774194   0.506835877862596]';

   KK = settings.KK;
   kc = settings.kc;
   fc = KK([1 5]);
   cc = KK(1:2,3);
   alpha_c = KK(1,2)/fc(1);
   p1 = pextend(normalize(x1,fc,cc,kc,alpha_c));
   p2 = pextend(normalize(x2,fc,cc,kc,alpha_c));
   U = intsec2views(P{ii},P{jj},p1,p2);
   
   %reproject(structure(U),motion(P([ii,jj])),{imagedata([],p1),imagedata([],p2)},[],'numbered');
   
   
   %add to u0 and U0
   nbr0 = size(U0,2);
   nbr = size(U,2);
   U0 = [U0,U];
   u0.pointnr = nbr0 + nbr;
   u0.points{ii}=[u0.points{ii},p1];
   u0.index{ii}=[u0.index{ii},nbr0+1:nbr0+nbr];
   u0.points{jj}=[u0.points{jj},p2];
   u0.index{jj}=[u0.index{jj},nbr0+1:nbr0+nbr];
   
   [U0,P0,lambda] = modbundle_sparse(U0,P0,u0,20,0.01);
   
   save C:\Users\fredrik\Documents\MATLAB\tunnels\experiment3_data\final\result_seq1 U0 P0 u0 settings
   
end

%SEQUENCE 2

%start with getting coordinate system of sequence 1
load C:\Users\fredrik\Documents\MATLAB\tunnels\experiment3_data\final\result_seq1

T1=[P0{100};[0,0,0,1]]; %im12730
T2=[P0{200};[0,0,0,1]]; %im10296

clear P0 u0 U0 settings;

load C:\Users\fredrik\Documents\MATLAB\tunnels\experiment3_data\data2\result_data

settings1 = settings;
u1 = u;
P1 = P;
U1 = U;
tmp=[P1{11};[0,0,0,1]]; %im12730
T = inv(tmp)*T1;
for ii=1:length(P1);
    P1{ii} = P1{ii}*T;
end
U1 = inv(T)*U1;

load C:\Users\fredrik\Documents\MATLAB\tunnels\experiment3_data\data7\result_data
settings2 = settings;
u2 = u;
P2 = P;
U2 = U;

clear P u U settings;

tmp=[P2{10};[0,0,0,1]]; %im10296
T = inv(tmp)*T2;
for ii=1:length(P2);
    P2{ii} = P2{ii}*T;
end
U2 = inv(T)*U2;

%merge
nbr1 = size(U1,2);
nbr2 = size(U2,2);

P0 = [P1,P2];
U0 = [U1,U2];

u0 = u1;
u0.pointnr = nbr1+nbr2;
u0.points = [u1.points,u2.points];
u0.sift = [u1.sift,u2.sift];
tmp = u2.index;
for ii=1:length(tmp);
    tmp{ii}=tmp{ii}+nbr1;
end
u0.index = [u1.index,tmp];

settings0 = settings1;
settings0.img_path = 'C:\Users\fredrik\Dropbox\tunnels\';
for ii=1:100,
    settings0.imnames(ii).name=['experiment3\',settings1.imnames(ii).name];
    settings0.imnames(ii+100).name=['experiment4\',settings2.imnames(ii).name];
end

%add manual points
   x1 = 1000*[
   0.958564516129032   0.509583969465649;...
   0.986951612903226   0.495843511450382;...
   1.026951612903226   0.479354961832061;...
   1.097919354838710   0.453248091603054;...
   1.234693548387097   0.394164122137405;...
   1.563725806451613   0.270500000000000;...
   0.345661290322581   0.315843511450382;...
   0.492758064516129   0.396912213740458;...
   0.570177419354839   0.439507633587786;...
   0.620500000000000   0.464240458015267;...
   0.654048387096774   0.482103053435115]';
   x2 = 1000*[
   1.055338709677420   0.528820610687023;...
   1.092758064516129   0.524698473282443;...
   1.154693548387097   0.515080152671756;...
   1.248887096774194   0.501339694656489;...
   1.402435483870968   0.479354961832061;...
   1.683725806451613   0.438133587786260;...
   0.688887096774193   0.405156488549619;...
   0.723725806451613   0.450500000000000;...
   0.744370967741936   0.477980916030534;...
   0.759854838709678   0.497217557251908;...
   0.768887096774194   0.506835877862596]';

   KK = settings1.KK;
   kc = settings1.kc;
   fc = KK([1 5]);
   cc = KK(1:2,3);
   alpha_c = KK(1,2)/fc(1);
   p1 = pextend(normalize(x1,fc,cc,kc,alpha_c));
   p2 = pextend(normalize(x2,fc,cc,kc,alpha_c));
   ii=6;%im12710
   jj=10+100; %im10296
   U = intsec2views(P0{ii},P0{jj},p1,p2);


   %reproject(structure(U),motion(P0([ii,jj])),{imagedata([],p1),imagedata([],p2)},[],'numbered');

   
   %add to u0 and U0
   nbr0 = size(U0,2);
   nbr = size(U,2);
   U0 = [U0,U];
   u0.pointnr = nbr0 + nbr;
   u0.points{ii}=[u0.points{ii},p1];
   u0.index{ii}=[u0.index{ii},nbr0+1:nbr0+nbr];
   u0.points{jj}=[u0.points{jj},p2];
   u0.index{jj}=[u0.index{jj},nbr0+1:nbr0+nbr];

    %add a couple of manual correspondences
    ii=91;jj=200;kk=100;
    
    figure(1);
    im = imread(strcat(settings0.img_path,settings0.imnames(ii).name));
    imagesc(im); hold on;
    figure(2);
    im = imread(strcat(settings0.img_path,settings0.imnames(jj).name));
    imagesc(im); hold on;
    figure(3);
    im = imread(strcat(settings0.img_path,settings0.imnames(kk).name));
    imagesc(im); hold on;
    
   x1 = 1000*[
   0.4366222366157701   0.4813542359051075;...
   0.5031345855720019   0.5379204625623677;...
   0.4450318439550638   0.5944866892196278;...
   0.6089376393637581   0.7651859726814506;...
   0.6727079466542082   0.6931792407925268;...
   0.7378645649727114   0.6272862880262474;...
   0.8989620612032534   0.5958223571811065;...
   0.9115204081632651   0.5987966281056361;...
   0.9169025568604132   0.5675667833980739;...
   1.121009554140128   0.524588397790055]';
   x2 = 1000*[
   0.6823611565333729   0.5257801128781376;...
   0.7144705663743124   0.5740293967649534;...
   0.694745723172628   0.636893442622951;...
   1.088898133748056   0.804106557377049;...
   1.005290046656299   0.729352459016393;...
   0.909737947122862   0.652631147540983;...
   0.993346034214619   0.609352459016393;...
   1.026192068429238   0.615254098360656;...
   1.026192068429238   0.579844262295082;...
   1.306876360808709   0.573942622950819]';
   x3 = 1000*[
   0.135430875576037   0.408255102040816;...
   0.303541474654378   0.521607871720117;...
   0.139854838709677   0.619217201166181;...
   0.387596774193548   1.009654518950438;...
   0.701698156682028   0.669596209912536;...
   0.905200460829493   0.606622448979592;...
   0.936168202764977   0.619217201166181;...
   0.940592165898617   0.575135568513120;...
   1.458195852534562   0.493269679300291]';
       

   KK = settings1.KK;
   kc = settings1.kc;
   fc = KK([1 5]);
   cc = KK(1:2,3);
   alpha_c = KK(1,2)/fc(1);
   p1 = pextend(normalize(x1,fc,cc,kc,alpha_c));
   p2 = pextend(normalize(x2,fc,cc,kc,alpha_c));
   p3 = pextend(normalize(x3,fc,cc,kc,alpha_c));
   U = intsec2views(P0{ii},P0{jj},p1,p2);
   
   %reproject(structure(U),motion(P0([ii,jj])),{imagedata([],p1),imagedata([],p2)},[],'numbered');
   
   
   %add to u0 and U0
   nbr0 = size(U0,2);
   nbr = size(U,2);
   U0 = [U0,U];
   u0.pointnr = nbr0 + nbr;
   u0.points{ii}=[u0.points{ii},p1];
   u0.index{ii}=[u0.index{ii},nbr0+1:nbr0+nbr];
   u0.points{jj}=[u0.points{jj},p2];
   u0.index{jj}=[u0.index{jj},nbr0+1:nbr0+nbr];
   u0.points{kk}=[u0.points{kk},p3];
   u0.index{kk}=[u0.index{kk},[nbr0+[1,2,3,5,6,7,8,9,10]]];
   
   lambda = 0.01;
   [U0,P0,lambda] = modbundle_sparse(U0,P0,u0,20,lambda);
   
 
   

   %save C:\Users\fredrik\Documents\MATLAB\tunnels\experiment3_data\final\result_seq2 U0 P0 u0 settings0
   
   

   
   
   
   
   
   
   
   
   
  
   
   
   
   
   
%SEQUENCE 3

%start with getting coordinate system of sequence 2
load C:\Users\fredrik\Documents\MATLAB\tunnels\experiment3_data\final\result_seq2

T1=[P0{100};[0,0,0,1]]; %im13086
T2=[P0{200};[0,0,0,1]]; %im10656

clear P0 u0 U0 settings0;

load C:\Users\fredrik\Documents\MATLAB\tunnels\experiment3_data\data3\result_data

settings1 = settings;
u1 = u;
P1 = P;
U1 = U;

clear P u U settings;

tmp=[P1{10};[0,0,0,1]]; %im13086
T = inv(tmp)*T1;
for ii=1:length(P1);
    P1{ii} = P1{ii}*T;
end
U1 = inv(T)*U1;

load C:\Users\fredrik\Documents\MATLAB\tunnels\experiment3_data\data8\result_data
settings2 = settings;
u2 = u;
P2 = P;
U2 = U;

clear P u U settings;

tmp=[P2{10};[0,0,0,1]]; %im10656
T = inv(tmp)*T2;
for ii=1:length(P2);
    P2{ii} = P2{ii}*T;
end
U2 = inv(T)*U2;

%merge
nbr1 = size(U1,2);
nbr2 = size(U2,2);

P0 = [P1,P2];
U0 = [U1,U2];

u0 = u1;
u0.pointnr = nbr1+nbr2;
u0.points = [u1.points,u2.points];
u0.sift = [u1.sift,u2.sift];
tmp = u2.index;
for ii=1:length(tmp);
    tmp{ii}=tmp{ii}+nbr1;
end
u0.index = [u1.index,tmp];

settings0 = settings1;
settings0.img_path = 'C:\Users\fredrik\Dropbox\tunnels\';
for ii=1:100,
    settings0.imnames(ii).name=['experiment3\',settings1.imnames(ii).name];
    settings0.imnames(ii+100).name=['experiment4\',settings2.imnames(ii).name];
end

    %add a couple of manual correspondences
    ii=1;jj=110;kk=10;
    
    figure(1);
    im = imread(strcat(settings0.img_path,settings0.imnames(ii).name));
    imagesc(im); hold on;
    figure(2);
    im = imread(strcat(settings0.img_path,settings0.imnames(jj).name));
    imagesc(im); hold on;
    figure(3);
    im = imread(strcat(settings0.img_path,settings0.imnames(kk).name));
    imagesc(im); hold on;
    
   x1 = 1000*[
   0.4366222366157701   0.4813542359051075;...
   0.5031345855720019   0.5379204625623677;...
   0.4450318439550638   0.5944866892196278;...
   0.6089376393637581   0.7651859726814506;...
   0.6727079466542082   0.6931792407925268;...
   0.7378645649727114   0.6272862880262474;...
   0.8989620612032534   0.5958223571811065;...
   0.9115204081632651   0.5987966281056361;...
   0.9169025568604132   0.5675667833980739;...
   1.121009554140128   0.524588397790055]';
   x2 = 1000*[
   0.6823611565333729   0.5257801128781376;...
   0.7144705663743124   0.5740293967649534;...
   0.694745723172628   0.636893442622951;...
   1.088898133748056   0.804106557377049;...
   1.005290046656299   0.729352459016393;...
   0.909737947122862   0.652631147540983;...
   0.993346034214619   0.609352459016393;...
   1.026192068429238   0.615254098360656;...
   1.026192068429238   0.579844262295082;...
   1.306876360808709   0.573942622950819]';
   x3 = 1000*[
   0.135430875576037   0.408255102040816;...
   0.303541474654378   0.521607871720117;...
   0.139854838709677   0.619217201166181;...
   0.387596774193548   1.009654518950438;...
   0.701698156682028   0.669596209912536;...
   0.905200460829493   0.606622448979592;...
   0.936168202764977   0.619217201166181;...
   0.940592165898617   0.575135568513120;...
   1.458195852534562   0.493269679300291]';
%point 4 missing in x3
figure(1);plot(imagedata([],x1),'numbered');
figure(2);plot(imagedata([],x2),'numbered');
figure(3);plot(imagedata([],x3),'numbered');
       

   KK = settings1.KK;
   kc = settings1.kc;
   fc = KK([1 5]);
   cc = KK(1:2,3);
   alpha_c = KK(1,2)/fc(1);
   p1 = pextend(normalize(x1,fc,cc,kc,alpha_c));
   p2 = pextend(normalize(x2,fc,cc,kc,alpha_c));
   p3 = pextend(normalize(x3,fc,cc,kc,alpha_c));
   U = intsec2views(P0{ii},P0{jj},p1,p2);
   
   %reproject(structure(U),motion(P0([ii,jj])),{imagedata([],p1),imagedata([],p2)},[],'numbered');
   
   %add to u0 and U0
   nbr0 = size(U0,2);
   nbr = size(U,2);
   U0 = [U0,U];
   u0.pointnr = nbr0 + nbr;
   u0.points{ii}=[u0.points{ii},p1];
   u0.index{ii}=[u0.index{ii},nbr0+1:nbr0+nbr];
   u0.points{jj}=[u0.points{jj},p2];
   u0.index{jj}=[u0.index{jj},nbr0+1:nbr0+nbr];
   u0.points{kk}=[u0.points{kk},p3];
   u0.index{kk}=[u0.index{kk},[nbr0+[1,2,3,5,6,7,8,9,10]]];
   
   ii=85;jj=200;
   
   figure(1);
   im = imread(strcat(settings0.img_path,settings0.imnames(ii).name));
   imagesc(im); hold on;
   figure(2);
   im = imread(strcat(settings0.img_path,settings0.imnames(jj).name));
   imagesc(im); hold on;
   
   x1 = 1000*[
   0.147642857142857   0.397778106508876;...
   0.153357142857143   0.623576923076923;...
   0.250500000000000   0.531979289940828;...
   0.473357142857143   0.536239644970414;...
   0.667642857142857   0.555411242603550;...
   0.724785714285714   0.525588757396450;...
   0.636214285714286   0.736476331360947;...
   0.687642857142857   0.681091715976331;...
   0.744785714285715   0.623576923076923;...
   1.530500000000000   0.698133136094675]';
   x2 = 1000*[
   0.515854107648725   0.479669960474308;...
   0.540330028328612   0.729393280632411;...
   0.586562322946176   0.612001976284585;...
   0.684466005665722   0.579986166007905;...
   0.760613314447592   0.571448616600791;...
   0.782369688385269   0.535164031620553;...
   1.078800283286119   0.808365612648221;...
   0.991774787535411   0.733662055335968;...
   0.899310198300283   0.654689723320158;...
   1.815797450424930   0.735796442687747]';


   KK = settings1.KK;
   kc = settings1.kc;
   fc = KK([1 5]);
   cc = KK(1:2,3);
   alpha_c = KK(1,2)/fc(1);
   p1 = pextend(normalize(x1,fc,cc,kc,alpha_c));
   p2 = pextend(normalize(x2,fc,cc,kc,alpha_c));
   U = intsec2views(P0{ii},P0{jj},p1,p2);


   %reproject(structure(U),motion(P0([ii,jj])),{imagedata([],p1),imagedata([],p2)},[],'numbered');

   
   %add to u0 and U0
   nbr0 = size(U0,2);
   nbr = size(U,2);
   U0 = [U0,U];
   u0.pointnr = nbr0 + nbr;
   u0.points{ii}=[u0.points{ii},p1];
   u0.index{ii}=[u0.index{ii},nbr0+1:nbr0+nbr];
   u0.points{jj}=[u0.points{jj},p2];
   u0.index{jj}=[u0.index{jj},nbr0+1:nbr0+nbr];
   
   
   A={};u={};
   for ii=1:length(P0);
       A{ii} = P0{ii}(:,1:3);
       tmp=NaN*ones(3,u0.pointnr);
       tmp(:,u0.index{ii})=u0.points{ii};
       u{ii}=tmp(1:2,:);
   end
    [U,P]=linf_knownrotation(u,A,10,0.01);       
   
    U0 = U; P0 = P;
    
   lambda = 0.01;
   [U0,P0,lambda] = modbundle_sparse(U0,P0,u0,20,lambda);
   
 
   

   %save C:\Users\fredrik\Documents\MATLAB\tunnels\experiment3_data\final\result_seq3 U0 P0 u0 settings0
   
   
   
   
   
   
   
   
   
   