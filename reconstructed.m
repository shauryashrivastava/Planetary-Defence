%% code to get reconstruction of lightcurves by incorporating 'n' dominant frequencies

clear variables
close all

load('PSDold.mat')
load('PSDnew.mat') % loading the python power spectra

n = 8;% specify the number of frequencies to be reconstructed here

fr=linspace(1e-6*pi, 1e-4*pi, 1000000);
fr=fr/(2*pi);
freq=zeros(300,10);

% maxf=zeros(300,1);
for i=1:300
    %importing the lightcurve data
    pathold = 'lightcurves/lcvold';
    pathnew = 'lightcurves/lcvnew';
    file=sprintf('%s%.3d.dat',pathold,i);
    filen=sprintf('%s%.3d.dat',pathnew,i);
    data=importdata(file);
    datan=importdata(filen);
    data(:,1)=data(:,1)*60;
    datan(:,1)=datan(:,1)*60;

    %% computing the power spectra

    %% pre collision frequency analysis
    % data10=data; data10(:,2)=data10(:,2)-mean(data10(:,2));
    % s1=[data10(:,2)' data10(:,2)' data10(:,2)' data10(:,2)' data10(:,2)' data10(:,2)' data10(:,2)' data10(:,2)' data10(:,2)' ];
    % ss1=[s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1];
    % t=zeros(size(ss1));
    % for j=1:length(data10):length(t)
    %     t(j:j+length(data10)-1)=data10(:,1)'+(data10(end,1)+1)*((j-1)/length(data10));
    % end
    % [pxx,fro] = plomb(ss1,t);
    % clear s1 ss1 loc

    % %% post collision frequency analysis
    % data10=datan; data10(:,2)=data10(:,2)-mean(data10(:,2));
    % s1=[data10(:,2)' data10(:,2)' data10(:,2)' data10(:,2)' data10(:,2)' data10(:,2)' data10(:,2)' data10(:,2)' data10(:,2)' ];
    % ss1=[s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1];
    % t=zeros(size(ss1));
    % for j=1:length(data10):length(t)
    %     t(j:j+length(data10)-1)=data10(:,1)'+(data10(end,1)+1)*((j-1)/length(data10));
    % end
    % [pxx2,frn] = plomb(ss1,t);

    % finding the first 'n' dominant frequencies 
    [~,loc]=findpeaks(powerold(i,:),'NPeaks',n,'Sortstr','descend');
    fo=fr(loc);
    [~,loc2]=findpeaks(powernew(i,:),'NPeaks',n,'Sortstr','descend');
    fn=fr(loc2);
    % generating time signal with no missing points
    time_old=data(1,1):600:data(end,2);
    time_new=datan(1,1):600:datan(end,2);
    % maxf(i)=fr(loc(1));

    %% fitting the reconstruction of first 'n' frequencies to the data
    % old data
    ft = 2*pi*data(:,1)*fo;
    ABC = [ones(size(ft(:,1))) cos(ft) sin(ft)] \ data(:,2);
    to=data(1,1):600:data(end,1);
    fx=2*pi*to'*fo;
    recon_old=[ones(size(fx(:,1))) cos(fx) sin(fx)] * ABC;
    clear ft ABC fx

    % new data
    ft = 2*pi*datan(:,1)*fn;
    ABC = [ones(size(ft(:,1))) cos(ft) sin(ft)] \ datan(:,2);
    tn=datan(1,1):600:datan(end,1);
    fx=2*pi*tn'*fn;
    recon_new=[ones(size(fx(:,1))) cos(fx) sin(fx)] * ABC;

    %% saving the reconstructed files
    pathre_old='reconstructed/rlcvold';
    pathre_new='reconstructed/rlcvnew';
    file_old=sprintf('%s%.3d.csv',pathre_old,i);
    file_new=sprintf('%s%.3d.csv',pathre_new,i);
    do=[to' recon_old]; dn=[tn' recon_new];
    writematrix(do,file_old);
    writematrix(do,file_new);

    freq(i,:)=[fo(1:5) fn(1:5)];
end


D=importdata('newfeaturelist_1903.csv');
d=D.data;
n=D.textdata;
dd=[d(:,1:3) freq d(:,14:end)];
featuretable=array2table(dd);
featuretable.Properties.VariableNames=n;
writetable(featuretable,'newfeaturelist_0109.csv');


% mfn=zeros(300,1);
% for i=1:300
%     [~,ind]=max(powernew(i,:));
%     mfn(i)=fr(ind);
% end

% r=((1./maxf)/(2*pi/sqrt(mu))).^(2/3);
% a=((1./mfn)/(2*pi/sqrt(mu))).^(2/3);
% msp=700;vsp=5000;ms=1.3e9;
% betac=(sqrt(2-r./a)-1).*sqrt(mu./r)*ms/(msp*vsp);
% rd=855.19;
% betacc=(sqrt(2-rd./a)-1).*sqrt(mu/rd)*ms/(msp*vsp);
% given  = importdata('lightcurves/parameters.csv');
% beta=given(1,:);