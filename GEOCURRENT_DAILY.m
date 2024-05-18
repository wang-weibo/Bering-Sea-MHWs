% READ GEO current
% 2023/8/25

clc
clear
close all

inpath='..\DAILY\';
lim=[50,70,160,205];
period=datetime(1993,1,1):datetime(2021,12,31);
period=datevec(period);

% READ lat lon
year=1993;
month=1;
day=1;
fpath=[inpath,sprintf('%04d',year),'\',sprintf('%02d',month),'\'];
fname=['dt_global_allsat_phy_l4_',sprintf('%04d',year),sprintf('%02d',month),sprintf('%02d',day),'*.nc'];
files=ls([fpath,fname]);
lat=ncread([fpath,files],'latitude');
lon=ncread([fpath,files],'longitude');
lon(lon<0)=lon(lon<0)+360;
[lon,lid]=sort(lon);

latin=lat<=lim(2) & lat>=lim(1);
lonin=lon<=lim(4) & lon>=lim(3);
lat=lat(latin);
lon=lon(lonin);
[sshlon,sshlat]=meshgrid(lon,lat);
clear lon lat

% load ref geo
load('ECMWFgeo.mat');
[lon,lat]=meshgrid(lon,lat);
[m,n]=size(lon);

% read geo-current
gu=NaN(n,m,size(period,1));
gv=NaN(n,m,size(period,1));
for i=1:size(period,1)
    year=period(i,1);
    month=period(i,2);
    day=period(i,3);
    fpath=[inpath,sprintf('%04d',year),'\',sprintf('%02d',month),'\'];
    fname=['dt_global_allsat_phy_l4_',sprintf('%04d',year),sprintf('%02d',month),sprintf('%02d',day),'*.nc'];
    files=ls([fpath,fname]);
    if(~isempty(files))
        u=ncread([fpath,files],'ugos');
        v=ncread([fpath,files],'vgos');
        u=u(lid,:);
        v=v(lid,:);
        u=u(lonin,latin);
        v=v(lonin,latin);
        u=interp2(sshlon,sshlat,u',double(lon),double(lat),'linear');
        v=interp2(sshlon,sshlat,v',double(lon),double(lat),'linear');
        gu(:,:,i)=u';
        gv(:,:,i)=v';
        clear u v
    else
        disp([fname,' has not founded! please check.'])
    end
end
save('geocurrent.mat','lon',"gv","gu","lat")

u1=squeeze(gu(:,:,10));
v1=squeeze(gv(:,:,10));
quiversc(lon,lat,u1',v1')
hold on
borders('countries','red','center',180)