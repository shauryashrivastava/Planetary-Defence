%% code to find rise and  interval times for each
% also adds other relevant features based on statistics of the light curve data
%
% Please refer to allfeatures.m for more information on the feature extraction algorithm

clear variables 
close all

sd=zeros(300,3);
max_slope=zeros(300,3);
MAD=zeros(300,3);
mbrp=zeros(300,3);
skew=zeros(300,3);
scatter_res=zeros(300,3);
varabove=zeros(300,3);
for i=1:300

    pathold = 'lightcurves/lcvold';
    pathnew = 'lightcurves/lcvnew';
    file=sprintf('%s%.3d.dat',pathold,i);
    filen=sprintf('%s%.3d.dat',pathnew,i);
    data=importdata(file);
    datan=importdata(filen);

    sd(i,1)=std(data(:,2));
    slope=(data(2:end,2)-data(1:end-1,2))./(data(2:end,1)-data(1:end-1,1));
    max_slope(i,1)=max(slope);
    MAD(i,1)=mad(data(:,2));
    fracmed=( (data(:,2)-median(data(:,2)))/median(data(:,2))) <0.10;
    mbrp(i,1)=length(find(fracmed==1));
    skew(i,1)=skewness(data(:,2));
    scatter_res(i,1)=mad(data(:,2)-mean(data(:,2)))/mean(data(:,2));
    abovedata=data(data(:,2)>mean(data(:,2)),2);
    varabove(i,1)=var(abovedata);

    data=datan;
    sd(i,2)=std(data(:,2));
    slope=(data(2:end,2)-data(1:end-1,2))./(data(2:end,1)-data(1:end-1,1));
    max_slope(i,2)=max(slope);
    MAD(i,2)=mad(data(:,2));
    fracmed=( (data(:,2)-median(data(:,2)))/median(data(:,2))) <0.10;
    mbrp(i,2)=length(find(fracmed==1));
    skew(i,2)=skewness(data(:,2));
    scatter_res(i,2)=mad(data(:,2)-mean(data(:,2)))/mean(data(:,2));
    abovedata=data(data(:,2)>mean(data(:,2)),2);
    varabove(i,2)=var(abovedata);


    sd(i,3)=sd(i,2)-sd(i,1);
    max_slope(i,3)=max_slope(i,2)-max_slope(i,1);
    MAD(i,3)=MAD(i,2)-MAD(i,1);
    mbrp(i,3)=mbrp(i,2)-mbrp(i,1);
    skew(i,3)=skew(i,2)-skew(i,1);
    scatter_res(i,3)=scatter_res(i,2)-scatter_res(i,1);
    varabove(i,3)=(varabove(i,2)-varabove(i,1))/varabove(i,1);

end


D=importdata('newfeaturelist_2peaks_1708.csv');
d=D.data;
n=D.textdata;
features=[d sd max_slope MAD mbrp skew scatter_res varabove]; 
featurename={'beta' 'J2' 'a/c' 'f1' 'f2' 'f3' 'f4' 'f5' 'fn1' 'fn2' 'fn3' 'fn4' 'fn5' 'r' 'a' 'beta_estimate' 'delv' 'delp' 'io' 'in' 'dto' 'dtn' 'deltadt' 'ato' 'atn' 'dat' 'mho' 'mhn' 'dmh' 'p2po' 'p2pn' 'dp2p' 'ilo' 'iln' 'iho' 'ihn'  'sdo' 'sdn' 'dsd' 'mso' 'msn' 'dms' 'mado' 'madn' 'dmad' 'mbrpo' 'mbrpn' 'dmbrp'  'skewo' 'skewn' 'dskew' 'sro' 'srn' 'dsrn' 'vabvo' 'vabvn' 'dfvabv'};
featuretable=array2table(features);
featuretable.Properties.VariableNames=featurename;
writetable(featuretable,'newfeaturelist_1903.csv');
