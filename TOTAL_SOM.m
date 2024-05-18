% SOM
% 2023/8/23

clc
clear 
close all

hfinpath='..\ECMWF\NHF\';
windpath='..\ECMWF\WIND\';
satpath='..\ECMWF\SAT\';
sstpath='..\SST\';
geopath='..\GEC\';
inpath='..\SHUM\';

% read net heat flux
load([hfinpath,'nhf.mat'],'nhf');

% load wind
load([windpath,'u10.mat'],'u10');
load([windpath,'v10.mat'],'v10');

% load sat
load([satpath,'t2m.mat'],'t2m','lat','lon');

% load geo current
load([geopath,'GEC.mat'],'u','v');

% load sst
load([sstpath,'sst.mat'],'sst');
load([sstpath,'lon_and_lat.mat']);

% load shum
load([inpath,'shum.mat'],'shum');

% merge all dataset
m1=length(lon);
n1=length(lat);
el_num=m1*n1;
m=length(lon_used);
n=length(lat_used);
sub=NaN(size(sst,3),el_num*6+m*n);
for i=1:size(sst,3)
    %sst
    a=squeeze(sst(:,:,i));
    sub(i,1:m*n)=reshape(a,m*n,1);
    
    %sat
    a=squeeze(t2m(:,:,i));
    sub(i,m*n+1:m*n+el_num)=reshape(a,el_num,1);

    %wind
    a=squeeze(u10(:,:,i));
    sub(i,m*n+el_num+1:m*n+2*el_num)=reshape(a,el_num,1);
    a=squeeze(v10(:,:,i));
    sub(i,m*n+2*el_num+1:m*n+3*el_num)=reshape(a,el_num,1);    

    %nhf
    a=squeeze(nhf(:,:,i));
    sub(i,m*n+3*el_num+1:m*n+4*el_num)=reshape(a,el_num,1);    

    % geo current
    a=squeeze(u(:,:,i));
    sub(i,m*n+4*el_num+1:m*n+5*el_num)=reshape(a,el_num,1); 
    a=squeeze(v(:,:,i));
    sub(i,m*n+5*el_num+1:m*n+6*el_num)=reshape(a,el_num,1);        
end

sData = som_data_struct(sub,'name','Net Heat Flux');
sDiris = som_normalize(sData,'var');
msize=[4,3];
sMap=som_randinit(sDiris,'msize',msize,'lattice','rect');
sMap=som_seqtrain(sMap,sDiris,'msize',msize,'lattice','rect','neigh','ep','radius_ini',4,'radius_fin',0.1,'alpha_type','linear','trainlen',100);
%[q,t] = som_quality(sMap,sDiris);

g_node=NaN(size(sst,3),1);
for i=1:size(sst,3)    
    dist_node=NaN(msize(1)*msize(2),1);
    b=sub(i,:);
    for num=1:msize(1)*msize(2)
        a=sMap.codebook(num,:);
        dist_node(num) = dist_E(a,b);
    end
    [~,g_node(i)]=min(dist_node);
end

per=NaN(msize(1)*msize(2),1);
for num=1:msize(1)*msize(2)
    gid=g_node==num;
    per(num)=sum(gid)/length(g_node);
end

%%
node=12;
a=sMap.codebook(node,:);
a=som_denormalize(a,sMap);
r_sst=reshape(a(1:m*n),[m,n]);
r_sat=reshape(a(m*n+1:m*n+el_num),[m1,n1]);
r_u10=reshape(a(m*n+el_num+1:m*n+2*el_num),[m1,n1]);
r_v10=reshape(a(m*n+2*el_num+1:m*n+3*el_num),[m1,n1]);
r_nhf=reshape(a(m*n+3*el_num+1:m*n+4*el_num),[m1,n1]);
r_u=reshape(a(m*n+4*el_num+1:m*n+5*el_num),[m1,n1]);
r_v=reshape(a(m*n+5*el_num+1:m*n+6*el_num),[m1,n1]);
save(['node_',num2str(node),'.mat'],"r_v","r_u","r_nhf","r_v10","r_u10","r_sat","r_sst",'lat','lon','lon_used','lat_used');
 
subplot(1,2,1)
contourf(lon_used,lat_used,r_sst')
clim([-2,2]);
cmocean('balance','pivot',0);
colorbar;
hold on
quiversc(lon,lat,r_u10',r_v10');
borders('countries','red','center',180)

subplot(1,2,2)
contourf(lon,lat,r_nhf')
cmocean('balance','pivot',0);
clim([-50,50])
colorbar;
hold on
streamslice(double(lon),double(lat),r_u',r_v');
borders('countries','red','center',180)

