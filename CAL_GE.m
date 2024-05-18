%calculate Ekman drift
% 2024/4/4

clc
clear
close all

sicpath='..\SIC\';
deppath='.\';
windpath='..\ECMWF\WIND\';
opath='..\ECMWF\GE\';
load([deppath,'Etop.mat']);
h=h';
glat=lat';
glon=lon';
lat=lat(:,1);
lon=lon(1,:);
lon=lon';
lid=island(glat,glon);

for iyr=1979:2023
    for imo=1:12
        prefix=[sprintf('%04d',iyr),'-',sprintf('%02d',imo)];
        load([sicpath,'sic-',prefix,'.mat'],'sic');
        sic=sic/100;
        sic(isnan(sic))=0;
        load([windpath,'u10-',prefix,'.mat'],'u10');
        load([windpath,'v10-',prefix,'.mat'],'v10');
        [Ue,Ve]=ekman_finite(glat,glon,u10,v10,h,sic);
        %[Ue,Ve]=ekman(glat,glon,u10,v10,'ci',sic);
        Ue(lid)=NaN;
        Ve(lid)=NaN;
        save([opath,'ge-',prefix,'.mat'],'lon','lat',"Ve","Ue");
        clear sic u10 v10 Ue Ve
    end
end
