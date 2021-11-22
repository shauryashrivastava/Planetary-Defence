%%  MATLAB code to extract meaningful features for a binary asteroid using lightcurve data 
% Input: light curve data samples 
% Output: meaningful features using fourier transform and signal processing 

% author: Shaurya Shrivastava (shaurya039@gmail.com)

clear variables
close all

%% Loading constants, data and other parameters
% In our implementation, we found the Python implementation of the Lomb-Scargle Periodogram to be 
% more accurate and precise compared to the MATLAB implementation. You may use the MATLAB plomb function instead.

load('PSDold.mat')
load('PSDnew.mat') % loading the python power spectra

fr=linspace(1e-6*pi, 1e-4*pi, 1000000);
fr=fr/(2*pi); % frequency vector
mu=6.67*3; msp=700;vsp=5000;ms=1.3e9; %spacecraft weight and other parameters


f=zeros(300,75);
wierd=[217; 218; 239; 254; 258; 262; 266; 278; 293; 294];
 % for the given data, these samples were outliers, according to some metrices (eg: variance> some constant or no of peaks < 6, in our case)
past=importdata('newfeaturelist_1903.csv');
abovetime=past.data(:,21:26);

harmonics=zeros(300,36); % primary harmonics in the light curve data
sd=zeros(300,3);         % standard deviation
max_slope=zeros(300,3);  % max slope
MAD=zeros(300,3);        % mean absolute deviation
mbrp=zeros(300,3);       % median buffer range percentage
skew=zeros(300,3);       % skewness
scatter_res=zeros(300,3);% MAD(y-y_m)/MAD(y)
varabove=zeros(300,3);   % variance of the fraction of data which is above the mean
diptime=zeros(300,6);    % time taken to go from mean to minima

% saving the features data (with training variables as well)
features(1:200,1:3)=given;
features(201:300,1:3)=NaN;

for i=1:300
    %importing the lightcurve data (available on http://kelvins.esa.int/planetary-defence/)
    pathold = 'lightcurves/lcvold';
    pathnew = 'lightcurves/lcvnew';
    file=sprintf('%s%.3d.dat',pathold,i);
    filen=sprintf('%s%.3d.dat',pathnew,i);
    data=importdata(file);
    datan=importdata(filen);
    data(:,1)=data(:,1)*60; % converting to seconds
    datan(:,1)=datan(:,1)*60;

    % finding the first 'n' dominant frequencies 
    [~,loc]=findpeaks(powerold(i,:),'NPeaks',8,'Sortstr','descend');
    fo=fr(loc);
    [~,loc2]=findpeaks(powernew(i,:),'NPeaks',8,'Sortstr','descend');
    fn=fr(loc2);
    % generating time signal with no missing points
    time_old=data(1,1):600:data(end,2);
    time_new=datan(1,1):600:datan(end,2);

    if i==278           % strange case! % can be modified to handle outliers
        fn=[fn fn(end)];
    end

    %% fitting the reconstruction of first 'n' frequencies to the data to get high frequencies
    % old data
    ft = 2*pi*data(:,1)*fo;
    ABC = [ones(size(ft(:,1))) cos(ft) sin(ft)] \ data(:,2);
    recon_old=[ones(size(ft(:,1))) cos(ft) sin(ft)] * ABC;
    clear ft ABC

    % new data
    ft = 2*pi*datan(:,1)*fn;
    ABC = [ones(size(ft(:,1))) cos(ft) sin(ft)] \ datan(:,2);
    recon_new=[ones(size(ft(:,1))) cos(ft) sin(ft)] * ABC;
    
    %% pre collision data analysis
    timediff=data(2:end,1)-data(1:end-1,1);
    ind=find(abs(timediff-600)==0);
    data10=data(ind+1,:);
    s1=[data10(:,2)' data10(:,2)' data10(:,2)' data10(:,2)' data10(:,2)' data10(:,2)' data10(:,2)' data10(:,2)' data10(:,2)' ];
    ss1=[s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1];
    t=zeros(size(ss1));
    for j=1:length(data10):length(t)
        t(j:j+length(data10)-1)=data10(:,1)'+(data10(end,1)+1)*((j-1)/length(data10));
    end
    [pxx,fr] = plomb(ss1,t);
    [~,indo]=max(pxx);
    [~,loc]=findpeaks(pxx,'NPeaks',8,'Sortstr','descend');
    features(i,4:11)=fr(loc)';
    r(i)=((1/fr(indo))/(2*pi/sqrt(mu)))^(2/3);
    clear s1 ss1 loc




    %% calculating the non-dominant frequencies


    %% pre collision frequency analysis
    data10=data; data10(:,2)=data10(:,2)-recon_old;
    data10(:,2)=data10(:,2)-mean(data10(:,2));
    s1=[data10(:,2)' data10(:,2)' data10(:,2)' data10(:,2)' data10(:,2)' data10(:,2)' data10(:,2)' data10(:,2)' data10(:,2)' ];
    ss1=[s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1];
    t=zeros(size(ss1));
    for j=1:length(data10):length(t)
        t(j:j+length(data10)-1)=data10(:,1)'+(data10(end,1)+1)*((j-1)/length(data10));
    end
    [pxx,fro] = plomb(ss1,t);
    clear s1 ss1 loc



    % %% post collision frequency analysis
    data10=datan; data10(:,2)=data10(:,2)-recon_new;
    data10(:,2)=data10(:,2)-mean(data10(:,2));
    s1=[data10(:,2)' data10(:,2)' data10(:,2)' data10(:,2)' data10(:,2)' data10(:,2)' data10(:,2)' data10(:,2)' data10(:,2)' ];
    ss1=[s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1];
    t=zeros(size(ss1));
    for j=1:length(data10):length(t)
        t(j:j+length(data10)-1)=data10(:,1)'+(data10(end,1)+1)*((j-1)/length(data10));
    end
    [pxx2,frn] = plomb(ss1,t);


    % finding the first 'n' dominant frequencies 
    [~,locho]=findpeaks(pxx,'NPeaks',10,'Sortstr','descend');
    
    [~,lochn]=findpeaks(pxx2,'NPeaks',10,'Sortstr','descend');


    ftotalo=[fo fro(locho)'];
    ftotaln=[fn frn(lochn)'];


    %% fitting the reconstruction with higher frequencies
    % old data
    ft = 2*pi*data(:,1)*ftotalo;
    ABC = [ones(size(ft(:,1))) cos(ft) sin(ft)] \ data(:,2);
    to=data(1,1):600:data(end,1);
    fx=2*pi*to'*ftotalo;
    totalrecon_old=[ones(size(fx(:,1))) cos(fx) sin(fx)] * ABC;

    nf=(length(ABC)-1)/2;
    Amp=zeros(nf,1);

    for j=2:2:length(ABC)
        Amp(j/2)=sqrt(ABC(j)^2+ABC(j+1)^2);
    end

    clear ft ABC fx

    % new data
    ft = 2*pi*datan(:,1)*ftotaln;
    ABC = [ones(size(ft(:,1))) cos(ft) sin(ft)] \ datan(:,2);
    tn=datan(1,1):600:datan(end,1);
    fx=2*pi*tn'*ftotaln;
    totalrecon_new=[ones(size(fx(:,1))) cos(fx) sin(fx)] * ABC;

    nf=(length(ABC)-1)/2;
    Ampn=zeros(nf,1);

    for j=2:2:length(ABC)
        Ampn(j/2)=sqrt(ABC(j)^2+ABC(j+1)^2);
    end

    erro=zeros(length(data),1);
    errn=zeros(length(datan),1);

    % calculating error in lightcurve after fourier fit into light curve data
    for k=1:length(data)
        ixo=find(data(k,1)==to);
        erro(k)=data(k,2)-totalrecon_old(ixo);
    end

    for k=1:length(datan)
        ixn=find(datan(k,1)==tn);
        errn(k)=datan(k,2)-totalrecon_new(ixn);
    end

    %% assigning frequency according to frequency rules
    % fo1: 2.5-2.9e-5  fn1: 2.1-2.7e-5
    % fo2: 1.8-2.5e-5  fn2: 1-2.1e-5
    % fo3: 3-3.5e-5    fn3: 3-4e-5
    % fo4: >3.5e-5     fn4: >4e-5
    % fo5: <1.8e-5     fn5: <1e-5
    % hfo and hfn: high frequency after subtracting the reconstruction

    % defining rules for pre-collision data
    rule1=find((fo>2.5e-5 & fo<2.9e-5));
    rule2=find((fo>1.8e-5 & fo<2.5e-5));
    rule3=find((fo>3e-5 & fo<3.5e-5));
    rule4=find(fo>3.5e-5);
    rule5=find(fo<1.8e-5);
    % defining rules for post-collision data
    rulen1=find((fn>2.1e-5 & fn<2.7e-5));
    rulen2=find((fn>1e-5 & fn<2.1e-5));
    rulen3=find((fn>3e-5 & fn<4e-5));
    rulen4=find(fn>4e-5);
    rulen5=find(fn<1e-5);

    % assigning frequencies for pre-collision data
    if ~isempty(rule1)
        f(i,1)=fo(rule1(1));
        fo(rule1(1))=0;
    end
    if ~isempty(rule2)
        f(i,2)=fo(rule2(1));
        fo(rule2(1))=0;
    end
    if ~isempty(rule3)
        f(i,3)=fo(rule3(1));
        fo(rule3(1))=0;
    end
    if ~isempty(rule4)
        f(i,4)=fo(rule4(1));
        fo(rule4(1))=0;
    end
    if ~isempty(rule5)
        f(i,5)=fo(rule5(1));
        fo(rule5(1))=0;
    end
    % checking if there are any other freq which don't obey the rules
    for j=1:5
        if f(i,j)==0
            ind=find(fo>0);
            f(i,j)=fo(ind(1));
            fo(ind(1))=0;
        end
    end

    % frequency assigning for post collision 
    if ~isempty(rulen1)
        f(i,6)=fn(rulen1(1));
        fn(rulen1)=0;
    end
    if ~isempty(rulen2)
        f(i,7)=fn(rulen2(1));
        fn(rulen2)=0;
    end
    if ~isempty(rulen3)
        f(i,8)=fn(rulen3(1));
        fn(rulen3)=0;
    end
    if ~isempty(rulen4)
        f(i,9)=fn(rulen4(1));
        fn(rulen4)=0;
    end
    if ~isempty(rulen5)
        f(i,10)=fn(rulen5(1));
        fn(rulen5)=0;
    end
    for j=6:10
        if f(i,j)==0 
            if j==6
                ind=find(fn<f(i,1));
                [~,idx]=min(abs(fn(ind)-2.1e-5));
                if ~isempty(idx)
                    f(i,j)=fn(ind(idx));
                else
                    f(i,6)=2.31e-5;
                end
            else
                ind=find(fn>0);
                if ~isempty(ind)
                    f(i,j)=fn(ind(1));
                    fn(ind(1))=0;
                else
                    f(i,j)=f(i,j-5)-0.5e-5;
                end
            end
        end 
    end

    % for hfo and hfn
    f(i,11)=fro(locho(1));
    f(i,12)=frn(lochn(1));


    if ~isempty(find(i==wierd))
         % for wierd cases exception
        f(i,6)=2.31e-5;
    end
    


    % adding the subtraction to the feature list
    f(i,13:17)= (1./f(i,6:10)-1./f(i,1:5))/60;
    f(i,18)=(1./f(i,12)-1./f(i,11))/60;

    % straight line regression 
    f(i,19:20)=polyfit(data(:,1),data(:,2),1);
    f(i,21:22)=polyfit(datan(:,1),datan(:,2),1);
    f(i,23)=f(i,21)-f(i,18);
    f(i,24)=f(i,22)-f(i,20);

    % amplitude of all frequencies
    f(i,25:27)=Amp(1:3);
    f(i,28:29)=Amp(2:3)/Amp(1);
    f(i,30:32)=Ampn(1:3);
    f(i,33:34)=Ampn(2:3)/Ampn(1);
    f(i,35:39)=f(30:34)-f(i,25:29);

    % statistics of the data left after fourier fit
    f(i,40)=var(erro);
    f(i,41)=rms(erro);
    f(i,42)=var(errn);
    f(i,43)=rms(errn);
    f(i,44:45)=f(i,42:43)./f(i,40:41)-1;

    % features about statistics not considered before
    f(i,46)=min(data(:,2)); % amplitude
    f(i,47)=min(datan(:,2));
    f(i,48)=(f(i,47)-f(i,46))/f(i,46);
    f(i,49)=kurtosis(data(:,2));
    f(i,50)=kurtosis(datan(:,2));
    f(i,51)=f(i,51)-f(i,50);
    f(i,52)=(prctile(data(:,2),95)-prctile(data(:,2),5))/median(data(:,2));
    f(i,53)=(prctile(datan(:,2),95)-prctile(datan(:,2),5))/median(datan(:,2));
    f(i,54)=f(i,53)/f(i,52)-1; % fractional change
    f(i,55)=max(data(:,2)-median(data(:,2)))/median(data(:,2));
    f(i,56)=max(datan(:,2)-median(datan(:,2)))/median(datan(:,2));
    f(i,57)=f(i,57)-f(i,56);
    slope=(totalrecon_old(2:end)-totalrecon_old(1:end-1))./((to(2:end)-to(1:end-1))');
    f(i,58)=max(slope);
    slope=(totalrecon_new(2:end)-totalrecon_new(1:end-1))./((tn(2:end)-tn(1:end-1))');
    f(i,59)=max(slope);
    f(i,60)=f(i,59)-f(i,58);

    % pre-collision dip analysis
    [dips,pos]=findpeaks(-totalrecon_old,to,'NPeaks',8,'Sortstr','descend','MinPeakHeight',-mean(totalrecon_old),'MinPeakDistance',12000);
    dips=-dips;
    dips=rmoutliers(dips);
    f(i,61)=mean(dips(dips>mean(dips)));
    timediff=(sort(pos));
    timediff=timediff(2:end)-timediff(1:end-1);    
    timediff=rmoutliers(timediff);
    f(i,62)=mean(timediff(timediff>=mean(timediff)));
    f(i,63)=mean(timediff(timediff<mean(timediff)));
    if isnan(f(i,63))
        f(i,63)=f(i,62);
    end
    f(i,64)=(f(i,62)+f(i,63))/2;
    [~,ppos]=findpeaks(totalrecon_old,to,'NPeaks',8,'Sortstr','descend','MinPeakHeight',-mean(totalrecon_old),'MinPeakDistance',12000);
    cc=0;dto=0;
    for j=1:length(ppos)
        fff=min(pos(pos>ppos(j)));
        if ~isempty(fff)
            dto=dto+fff-ppos(j);
            cc=cc+1;
        end
    end
    diptime(i,1)=dto;

    %% post-collision dip analysis
    [dips,pos]=findpeaks(-totalrecon_new,tn,'NPeaks',8,'Sortstr','descend','MinPeakHeight',-mean(totalrecon_new),'MinPeakDistance',12000);
    dips=-dips;
    dips=rmoutliers(dips);
    f(i,65)=mean(dips(dips>mean(dips)));
    timediff=(sort(pos));
    timediff=timediff(2:end)-timediff(1:end-1);    
    timediff=rmoutliers(timediff);
    f(i,66)=mean(timediff(timediff>=mean(timediff)));
    f(i,67)=mean(timediff(timediff<mean(timediff)));
    if isnan(f(i,67))
        f(i,67)=f(i,66)-2400;
    end
    f(i,68)=(f(i,66)+f(i,67))/2;
    [~,ppos]=findpeaks(totalrecon_new,tn,'NPeaks',8,'Sortstr','descend','MinPeakHeight',-mean(totalrecon_new),'MinPeakDistance',12000);
    cc=0;dto=0;
    for j=1:length(ppos)
        fff=min(pos(pos>ppos(j)));
        if ~isempty(fff)
            dto=dto+fff-ppos(j);
            cc=cc+1;
        end
    end
    dto=dto/cc;
    diptime(i,2)=dto;
    % change in orbital parameters
    f(i,69:72)=f(i,65:68)-f(i,61:64);


    % radius of orbit and beta estimate
    f(i,73)= ((1/f(i,1))/(2*pi/sqrt(mu)))^(2/3);   % radius of precollision orbit
    f(i,74)= ((1/f(i,6))/(2*pi/sqrt(mu)))^(2/3);   % semi-major axis of post-collision orbit
    f(i,75)= (sqrt(2-f(i,73)/f(i,74))-1)*sqrt(mu/f(i,73))*ms/(msp*vsp);


    %% other statistics of light curves from intervaltime.m code

    % pre-collision data
    sd(i,1)=std(data(:,2));
    slope=(data(2:end,2)-data(1:end-1,2))./(data(2:end,1)-data(1:end-1,1));
    timediff=data(2:end,1)-data(1:end-1,1);
    slope=slope(timediff<(50*60));
    max_slope(i,1)=max(slope);
    MAD(i,1)=mad(data(:,2));
    fracmed=( (data(:,2)-median(data(:,2)))/median(data(:,2))) <0.10;
    mbrp(i,1)=length(find(fracmed==1));
    skew(i,1)=skewness(data(:,2));
    scatter_res(i,1)=mad(data(:,2)-mean(data(:,2)))/mean(data(:,2));
    abovedata=data(data(:,2)>mean(data(:,2)),2);
    varabove(i,1)=var(abovedata);

    % post-collision data
    sd(i,2)=std(datan(:,2));
    slope=(datan(2:end,2)-datan(1:end-1,2))./(datan(2:end,1)-datan(1:end-1,1));
    timediff=datan(2:end,1)-datan(1:end-1,1);
    slope=slope(timediff<(50*60));
    max_slope(i,2)=max(slope);
    MAD(i,2)=mad(datan(:,2));
    fracmed=( (datan(:,2)-median(datan(:,2)))/median(datan(:,2))) <0.10;
    mbrp(i,2)=length(find(fracmed==1));
    skew(i,2)=skewness(datan(:,2));
    scatter_res(i,2)=mad(datan(:,2)-mean(datan(:,2)))/mean(datan(:,2));
    abovedatan=datan(datan(:,2)>mean(datan(:,2)),2);
    varabove(i,2)=var(abovedatan);


    sd(i,3)=sd(i,2)-sd(i,1);
    max_slope(i,3)=max_slope(i,2)-max_slope(i,1);
    MAD(i,3)=MAD(i,2)-MAD(i,1);
    mbrp(i,3)=mbrp(i,2)-mbrp(i,1);
    skew(i,3)=skew(i,2)-skew(i,1);
    scatter_res(i,3)=scatter_res(i,2)-scatter_res(i,1);
    varabove(i,3)=(varabove(i,2)-varabove(i,1))/varabove(i,1);
    



    % harmonics reconstruction and features
    clear ABC
    fo=f(i,1:3);
    fn=f(i,6:8);
    ho=[fo 2*fo 3*fo];
    hn=[fn 2*fn 3*fn];
    
    % old data
    fho = 2*pi*datan(:,1)*ho;
    ABC = [ones(size(fho(:,1))) cos(fho) sin(fho)] \ datan(:,2);
    basephase=atan(ABC(2)/ABC(1));

    for j=2:2:length(ABC)
        harmonics(i,j/2)=sqrt(ABC(j)^2+ABC(j+1)^2);
        harmonics(i,9+j/2)= atan(ABC(j+1)/ABC(j)) - (ceil(j/6))*ho(mod(j/2,3)+1)*basephase/ho(1);
    end
    clear ABC

    % new data
    fhn = 2*pi*datan(:,1)*hn;
    ABC = [ones(size(fhn(:,1))) cos(fhn) sin(fhn)] \ datan(:,2);


    for j=2:2:length(ABC)
        harmonics(i,18+j/2)=sqrt(ABC(j)^2+ABC(j+1)^2);
        harmonics(i,27+j/2)=atan(ABC(j+1)/ABC(j)) - (ceil(j/6))*hn(mod(j/2,3)+1)*basephase/hn(1);
    end


    % above times


    [~,pos]=findpeaks(-totalrecon_old,'NPeaks',8,'Sortstr','descend','MinPeakHeight',-mean(totalrecon_old),'MinPeakDistance',10);

    pre_pts=zeros(12,1); post_pts=zeros(12,1);
    for j=1:length(pos)
        tdiff = (to-to(pos(j)))';
        idg = find(tdiff>0 &  totalrecon_old>mean(totalrecon_old));
        idl = find(tdiff<0 &  totalrecon_old>mean(totalrecon_old));
        if ~isempty(idg)
            % disp(j);
            % disp(data(idg(greater),1));
            post_pts(j)=idg(1);
        end
        if ~isempty(idl)
            pre_pts(j)=idl(end);
        end
    end

    pre=unique(sort(pre_pts));
    post=unique(sort(post_pts));


    ipre=find(pre>0);
    ipost=find(post>0);

    pre=pre(ipre(1):end);
    post=post(ipost(1):end);


    if length(pre)>length(post)
        pre=pre(1:end-1);
    elseif length(pre)<length(post)
        post=post(2:end);
    end

    timediff=to(pre(2:end))-to(post(1:end-1));
    timediff=timediff(timediff>0);
    diptime(i,3)=mean(timediff);



    %% post collision data
    [~,pos]=findpeaks(-totalrecon_new,'NPeaks',8,'Sortstr','descend','MinPeakHeight',-mean(totalrecon_new),'MinPeakDistance',10);

    pre_pts=zeros(12,1); post_pts=zeros(12,1);
    for j=1:length(pos)
        tdiff = (tn-tn(pos(j)))';
        idg = find(tdiff>0 &  totalrecon_new>mean(totalrecon_new));
        idl = find(tdiff<0 &  totalrecon_new>mean(totalrecon_new));
        if ~isempty(idg)
            % disp(j);
            % disp(data(idg(greater),1));
            post_pts(j)=idg(1);
        end
        if ~isempty(idl)
            pre_pts(j)=idl(end);
        end
    end

    pre=unique(sort(pre_pts));
    post=unique(sort(post_pts));


    ipre=find(pre>0);
    ipost=find(post>0);

    pre=pre(ipre(1):end);
    post=post(ipost(1):end);


    if length(pre)>length(post)
        pre=pre(1:end-1);
    elseif length(pre)<length(post)
        post=post(2:end);
    end

    timediff=tn(pre(2:end))-tn(post(1:end-1));
    timediff=timediff(timediff>0);
    diptime(i,4)=mean(timediff);
    diptime(i,5)=diptime(i,2)-diptime(i,1);
    diptime(i,6)=diptime(i,4)-diptime(i,3);

end

% saving all data into a new file
D=importdata('newfeaturelist_1903.csv');
d=D.data(:,1:3);
features=[d f sd max_slope MAD mbrp skew scatter_res varabove harmonics diptime abovetime];
featurename={'beta' 'J2' 'a/c' 'f1' 'f2' 'f3' 'f4' 'f5' 'fn1' 'fn2' 'fn3' 'fn4' 'fn5' 'hfo' 'hfn' 'df1' 'df2' 'df3' 'df4' 'df5' 'dhf' 'mo' 'co' 'mn' 'cn' 'dm' 'dc' 'afo1' 'afo2' 'afo3' 'afo21' 'afo31' 'afn1' 'afn2' 'afn3' 'afn21' 'afn31' 'daf1' 'daf2' 'daf3' 'daf21' 'daf31' 'rvaro' 'rmseo' 'rvarn' 'rmsen' 'drvar' 'drmse' 'ao' 'an' 'dfa' 'ko' 'kn' 'dk' 'pdfpo' 'pdfpn' 'dfpdfp' 'pao' 'pan' 'dpa' 'rmso' 'rmsn' 'drms' 'bo' 'to1' 'to2' 'mto' 'bn' 'tn1' 'tn2' 'mtn' 'db' 'dto1' 'dto2' 'dmt' 'r' 'a' 'betac' 'sdo' 'sdn' 'dsd' 'mso' 'msn' 'dms' 'mado' 'madn' 'dmad' 'mbrpo' 'mbrpn' 'dmbrp'  'skewo' 'skewn' 'dskew' 'sro' 'srn' 'dsrn' 'vabvo' 'vabvn' 'dfvabv' 'a11o' 'a21o' 'a31o' 'a12o' 'a22o' 'a32o' 'a13o' 'a23o' 'a33o' 'p11o' 'p21o' 'p31o' 'p12o' 'p22o' 'p32o' 'p13o' 'p23o' 'p33o' 'a11n' 'a21n' 'a31n' 'a12n' 'a22n' 'a32n' 'a13n' 'a23n' 'a33n' 'p11n' 'p21n' 'p31n' 'p12n' 'p22n' 'p32n' 'p13n' 'p23n' 'p33n' 'rdto' 'rdtn' 'rato' 'ratn' 'rddt' 'rdat' 'dto' 'dtn' 'ddt' 'ato' 'atn' 'dat'};
featuretable=array2table(features);
featuretable.Properties.VariableNames=featurename;
writetable(featuretable,'allfeaturelist.csv');