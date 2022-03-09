# MatLab-NetCDF
MatLab-NetCDF: Chlorophyll, Temperature and Salinity projects source code!

### Chlorophyll Visualization
``` matlab
close;clear all;clc;
%% ��ȡҶ����Chlorophyll NetCDF����
datadir='D:\Documents\MATLAB\20220307\CHL-MOD-SR4.6-DAY-20160101-20171231\'; 
filelist=dir([datadir,'*.nc']);
n=length(filelist);
% �鿴NetCDF�ļ���ϸ����
file1=[datadir,filelist(1).name];
ncdisp(file1);

% �趨Ŀ�꾭γ�ȷ�Χ
lon_min=105;lon_max=110;
lat_min=17;lat_max=23;

for s=1:n
    % ƴ�Ӿ���·�����ļ���Ϊ�ļ�ȫ·��
    filename=[datadir,filelist(s).name];
    % ���ļ���ȡNC ID
    nc_id=netcdf.open(filename,'NC_NOWRITE');
    % ��ȡNetCDF�ļ��о��������sΪ����ά��¼
    CHM(:,:,s)= ncread(filename,'CHL1_mean');
    % ��ȡNetCDF�ļ��о��Ⱥ�ά�ȸ����±���
    lon_nc  = ncread(filename,'lon');
    lat_nc  = ncread(filename,'lat');
    % ����NC ID�ر�NetCDF�ļ�
    netcdf.close(nc_id);
end;

lon_reg_index=find(lon_nc>= lon_min & lon_nc<=lon_max);
lat_reg_index=find(lat_nc>=lat_min& lat_nc<=lat_max);
lon_reg=lon_nc(lon_reg_index);
lat_reg=lat_nc(lat_reg_index);
CHMY(:,:,:)=CHM(lon_reg_index,lat_reg_index,:);
%% ���ݵ�����ѡȡ����
band_file='D:\Documents\MATLAB\20220307\etopo1_bedrock.nc'
ncdisp(band_file);

lon2_nc=(ncread(band_file,'lon'));
lat2_nc=flipud(ncread(band_file,'lat'));
band_nc=fliplr(ncread(band_file,'Band1'));

lon2_reg_index=find(lon2_nc>lon_min&lon2_nc<=lon_max);
lat2_reg_index=find(lat2_nc>lat_min&lat2_nc<=lat_max);

lon2_reg=lon2_nc(lon2_reg_index);
lat2_reg=lat2_nc(lat2_reg_index);
band_reg=band_nc(lon2_reg_index,lat2_reg_index);
%contourf(lon2_nc,lat2_nc,band_reg')
band_new=band_reg(1:2.5:end,1:2.5:end);
%contourf(x,y,ban')

st=ones(size(band_new,1),size(band_new,2));
for i=1:size(band_new,1);
    for j=1:size(band_new,2);
        if band_new(i,j)>=-100;
            st(i,j)=nan;
        end;
    end;
end;
for i=1:n;
    CHLL(:,:,i)=CHMY(:,:,i).*st;
end
chlaa=nanmean(CHLL(:,:,:),3);

% m_proj��Ҫ����M_Map�����ڻ��Ƶ�ͼ
lon=lon_reg;lat=lat_reg;
m_proj('Equidistant Cylindrical','lat',[lat_min lat_max],'lon',[lon_min lon_max]);
hold on;
m_contourf(lon,lat,chlaa',30,'linestyle','none');
caxis([0 0.3]);
colorbar;
m_coast('linewidth',1,'color','k');
m_coast('patch',[0.783,0.741,0.721]);
m_grid('linestyle','none','box','fancy');
gtext('(mg/m^{3})');

% ʱ�����л�
abc(:,:,:)=nanmean(CHLL(:,:,:),2);
chll(:,:)=squeeze(nanmean(abc(:,:,:),1));
%% �����ڷֲ�ͼ
l=1
for k=1:11:14
    chl1(:,:,l)=nanmean(CHLL(:,:,k:(k+1)),3);%����
    chl2(:,:,l)=nanmean(CHLL(:,:,(k+2):(k+4)),3);%����
    chl3(:,:,l)=nanmean(CHLL(:,:,(k+5):(k+7)),3);%�＾
    chl4(:,:,l)=nanmean(CHLL(:,:,(k+8):(k+10)),3);%�ļ�
    l=l+1;
end
%% ���� Winter
chl_win(:,:,:)=nanmean(chl1(:,:,:),3);
m_proj('Equidistant Cylindrical','lat',[lat_min lat_max],'lon',[lon_min lon_max]);
hold on;
m_contourf(lon,lat,chl_win',30,'linestyle','none');
caxis([0 0.3]);
colorbar;
hold on;
m_coast('linewidth',1,'color','k');
m_coast('patch',[0.783,0.741,0.721]);
m_grid('linestyle','none','box','fancy');
gtext('(mg/m^{3})');
%% ���� Spring
chl_spr(:,:,:)=nanmean(chl2(:,:,:),3);
m_proj('Equidistant Cylindrical','lat',[lat_min lat_max],'lon',[lon_min lon_max]);
hold on;
m_contourf(lon,lat,chl_spr',30,'linestyle','none');
caxis([0 0.3]);
colorbar;
hold on;
m_coast('linewidth',1,'color','k');
m_coast('patch',[0.783,0.741,0.721]);
m_grid('linestyle','none','box','fancy');
gtext('(mg/m^{3})');
%% �＾ Autumn
chl_aut(:,:,:)=nanmean(chl3(:,:,:),3);
m_proj('Equidistant Cylindrical','lat',[lat_min lat_max],'lon',[lon_min lon_max]);
hold on;
m_contourf(lon,lat,chl_aut',30,'linestyle','none');
caxis([0 0.3]);
colorbar;
hold on;
m_coast('linewidth',1,'color','k');
m_coast('patch',[0.783,0.741,0.721]);
m_grid('linestyle','none','box','fancy');
gtext('(mg/m^{3})');
%% �ļ� Summer
chl_sum(:,:,:)=nanmean(chl4(:,:,:),3);
m_proj('Equidistant Cylindrical','lat',[lat_min lat_max],'lon',[lon_min lon_max]);
hold on;
m_contourf(lon,lat,chl_sum',30,'linestyle','none');
caxis([0 0.3]);
colorbar;
hold on;
m_coast('linewidth',1,'color','k');
m_coast('patch',[0.783,0.741,0.721]);
m_grid('linestyle','none','box','fancy');
gtext('(mg/m^{3})');
```