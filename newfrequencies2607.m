% code for finding dominant frequencies in old and new data
% run this first, followed by newfeatures_1903.m and then allfeatures.m

clear variables
close all
%% declaring file path
tic
pathold = 'lightcurves/lcvold';
pathnew = 'lightcurves/lcvnew';
given  = importdata('lightcurves/parameters.csv');
features=zeros(300,26);
a=zeros(300,1); r=zeros(300,1); mu=6.67*3;
betac=zeros(300,1);
msp=700;vsp=5000;ms=1.3e9; %spacecraft weight and other parameters
features(1:200,1:3)=given;
features(201:300,1:3)=NaN;
radius=zeros(300,2);
%% finding all the dominant frequencies
for i=1:300

    %importing the lightcurve data
    file=sprintf('%s%.3d.dat',pathold,i);
    filen=sprintf('%s%.3d.dat',pathnew,i);
    data=importdata(file);
    datan=importdata(filen);
    
    %data simplification
    intensity_old=min(data(:,2));
    intensity_new=min(datan(:,2));
    data(:,1)=data(:,1)*60; % converting minutes to seconds
    datan(:,1)=datan(:,1)*60;
    data(:,2)=data(:,2)-mean(data(:,2)); % subtracting the mean from the data for fft
    datan(:,2)=datan(:,2)-mean(datan(:,2));

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

    %% post collision data analysis
    data=datan;
    % finding datapoints where the successive difference is 10 minutes
    timediff=data(2:end,1)-data(1:end-1,1);
    ind=find(abs(timediff-600)==0);
    data10=data(ind+1,:);
    s1=[data10(:,2)' data10(:,2)' data10(:,2)' data10(:,2)' data10(:,2)' data10(:,2)' data10(:,2)' data10(:,2)' data10(:,2)' ];
    ss1=[s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1 s1];
    t=zeros(size(ss1));
    for j=1:length(data10):length(t)
        t(j:j+length(data10)-1)=data10(:,1)'+(data10(end,1)+1)*((j-1)/length(data10));
    end
    [pxx,fr2] = plomb(ss1,t);
    [~,indn]=max(pxx);
    a(i)=( (1/fr2(indn)) /( 2*pi/sqrt(mu))  )^(2/3);
    [~,loc]=findpeaks(pxx,'NPeaks',8,'Sortstr','descend');
    features(i,12:19)=fr2(loc)';

    %other features
    betac(i)=(sqrt(2-r(i)/a(i))-1)*sqrt(mu/r(i))*ms/(msp*vsp);
    delv=betac(i)*msp*vsp/ms;

    if fr2(indn)<fr(indo)
        delp=((1/fr2(indn))-(1/fr(indo)));
    else
        indn=loc(fr2(loc)<fr(indo));
        if isempty(indn)
            [~,loc]=findpeaks(pxx,'NPeaks',25,'Sortstr','descend'); % considering more no of peaks
            indn=loc(fr2(loc)<fr(indo));
            delp=((1/fr2(indn(1)))-(1/fr(indo)));
        else
            delp=((1/fr2(indn(1)))-(1/fr(indo)));
        end
    end
    clear s1 ss1 loc

    %finding orbit radius from hohmann transfer equations using beta values
    % if i<=200 %checking if data includes beta
    %     options = optimoptions('fsolve','Display','iter','TolFun',1e-30,'TolX',1e-30);
    %     delv=features(i,1)*msp*vsp/ms;
    %     guess=370:50:50*370;
    %     for j=1:length(guess)
    %         x = fsolve(@(x)root2d(x,delp,delv),[guess(j) 2*guess(j)],options);
    %         if(imag(x)<1e-6)
    %             hohmann_radius=x;
    %             zzzz=guess(j);
    %         end
    %     end
    %     radius(i,1)=real(hohmann_radius(1));
    %     radius(i,2)=real(hohmann_radius(2));
    % else  %if beta not given use estimate of beta to calculate radius and delv
    %     delv=betac(i)*msp*vsp/ms;
    %     radius(i,1)=r(i);
    %     radius(i,2)=a(i);
    % end
    
    
    %saving other relevant features
    features(i,20)=betac(i);
    features(i,21)=delv;
    features(i,22)=intensity_old;
    features(i,23)=intensity_new;
    features(i,24)=delp;
    features(i,25)=r(i);
    features(i,26)=a(i);
end

featuretable=array2table(features);
featurelist={'beta' 'J2' 'a/c' 'f1' 'f2' 'f3' 'f4' 'f5' 'f6' 'f7' 'f8' 'fn1' 'fn2' 'fn3' 'fn4' 'fn5' 'fn6' 'fn7' 'fn8' 'beta_estimate' 'delv' 'io' 'in' 'delp' 'r' 'a'};
featuretable.Properties.VariableNames=featurelist;
%saving the data as a CSV file 

%% function to calculate roots of the hohmann transfer equations
function F =root2d(x,delp,delv)
    mu=6.67*3;
    k=(delv*sqrt(x(1)/mu)+1)^2;
    F(1)=k*(x(1)+x(2))-2*x(2);
    F(2)=2*pi*( ((x(1)+x(2))/2)^(3/2) - x(1)^(3/2) )/sqrt(mu) - delp;
end