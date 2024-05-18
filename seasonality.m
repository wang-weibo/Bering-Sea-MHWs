%% An example analysing seasonality of MHWs
% Here we provide an example about seasonality of MHWs
clc
clear
close all
sstinpath='..\SST\';
%% 1. Loading data

% Load NOAA OI SST V2 data
sst_full=NaN(200,60,datenum(2022,12,31)-datenum(1982,1,1)+1);
for i=1982:2022
    file_here=[sstinpath,'sst_' num2str(i)];
    load(file_here);
    eval(['sst_' num2str(i) '=sst;'])
    clear sst
    eval(['data_here=sst_' num2str(i) ';'])
    sst_full(:,:,(datenum(i,1,1):datenum(i,12,31))-datenum(1982,1,1)+1)=data_here;
end

load([sstinpath,'lon_and_lat.mat']);

%% 2. Detecting MHWs and MCSs

% Here we detect marine heatwaves in the Bering Sea based on the
% traditional definition of MHWs (Hobday et al. 2016). We detected MHWs
% during 1993 to 2022 for climatologies and thresholds in 1993 to 2005.

[MHW,mclim,m90,mhw_ts]=detect(sst_full,datenum(1982,1,1):datenum(2022,12,31),datenum(1993,1,1),datenum(2005,12,31),datenum(1993,1,1),datenum(2022,12,31)); %take about 30 seconds.


%% 3. Generating monthly and seasonal MHW metrics
% Here we calculate monthly and seasonal MHW metrics including numbers of
% MHW days and mean MHW intensity

% Generating date matrix
date_used=datevec(datenum(1993,1,1):datenum(2022,12,31));

% Determining land index
land_index=isnan(nanmean(mhw_ts,3));

% Monthly
mhwday_month=NaN(size(mhw_ts,1),size(mhw_ts,2),12); % lon-lat-month
mhwint_month=NaN(size(mhw_ts,1),size(mhw_ts,2),12); % lon-lat-month
for i=1:12
    index_used=date_used(:,2)==i;
    mhwday_month(:,:,i)=sum(~isnan(mhw_ts(:,:,index_used)),3,'omitnan')./(2022-1993+1);
    mhwint_month(:,:,i)=mean(mhw_ts(:,:,index_used),3,'omitnan');
end
mhwday_month(repmat(land_index,1,1,12))=nan;
% mhwday_month is the average number of MHW days in each month during
% 1993-2016
% mhwint_month is the average intensity of MHW days in each month during
% 1993-2016

% Seasonal
% Determining austral seasons
% AUT-SON WIN-DJF SPR-MAM SUM-JJA
seas=[9 10 11;...
    12 1 2;...
    3 4 5;...
    6 7 8];
mhwday_seas=NaN(size(mhw_ts,1),size(mhw_ts,2),4); % lon-lat-seasons
mhwint_seas=NaN(size(mhw_ts,1),size(mhw_ts,2),4); % lon-lat-seasons
for i=1:4
    index_used=ismember(date_used(:,2),seas(i,:));
    mhwday_seas(:,:,i)=sum(~isnan(mhw_ts(:,:,index_used)),3,'omitnan')./(3*(2022-1993+1));
    mhwint_seas(:,:,i)=mean(mhw_ts(:,:,index_used),3,'omitnan');
end
mhwday_seas(repmat(land_index,1,1,4))=nan;
% mhwday_seas (days/month) is the average number of MHW days in each season during
% 1993-2016
% mhwint_seas (^{o}C) is the average intensity of MHW days in each season during
% 1993-2016

%% seasonal MHW matrix
% Generating date matrix
date_used=datevec(datenum(1993,1,1):datenum(2022,12,31));

% Determining land index
land_index=isnan(nanmean(mhw_ts,3));

% Monthly
mhwday_monthly=NaN(size(mhw_ts,1),size(mhw_ts,2),12*(2022-1993+1)); % lon-lat-month
mhwint_monthly=NaN(size(mhw_ts,1),size(mhw_ts,2),12*(2022-1993+1)); % lon-lat-month
flag=1;
for year=1993:2022
    for month=1:12
        index_used=date_used(:,1)==year & date_used(:,2)==month;
        mhwday_monthly(:,:,flag)=sum(~isnan(mhw_ts(:,:,index_used)),3,'omitnan')./sum(index_used);
        mhwint_monthly(:,:,flag)=mean(mhw_ts(:,:,index_used),3,'omitnan');
        flag=flag+1;
    end
end
mhwday_monthly(repmat(land_index,1,1,12*(2022-1993+1)))=nan;

