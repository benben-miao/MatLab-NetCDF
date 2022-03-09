# MatLab-NetCDF
MatLab-NetCDF: Chlorophyll, Temperature and Salinity projects source code!

### Chlorophyll Visualization
``` matlab
close;clear all;clc;
%% 读取叶绿素Chlorophyll NetCDF数据
datadir='D:\Documents\MATLAB\20220307\CHL-MOD-SR4.6-DAY-20160101-20171231\'; 
filelist=dir([datadir,'*.nc']);
n=length(filelist);
% 查看NetCDF文件详细描述
file1=[datadir,filelist(1).name];
ncdisp(file1);

% 设定目标经纬度范围
lon_min=105;lon_max=110;
lat_min=17;lat_max=23;

for s=1:n
    % 拼接绝对路径和文件名为文件全路径
    filename=[datadir,filelist(s).name];
    % 打开文件获取NC ID
    nc_id=netcdf.open(filename,'NC_NOWRITE');
    % 读取NetCDF文件中具体变量，s为第三维记录
    CHM(:,:,s)= ncread(filename,'CHL1_mean');
    % 读取NetCDF文件中经度和维度赋予新变量
    lon_nc  = ncread(filename,'lon');
    lat_nc  = ncread(filename,'lat');
    % 根据NC ID关闭NetCDF文件
    netcdf.close(nc_id);
end;

lon_reg_index=find(lon_nc>= lon_min & lon_nc<=lon_max);
lat_reg_index=find(lat_nc>=lat_min& lat_nc<=lat_max);
lon_reg=lon_nc(lon_reg_index);
lat_reg=lat_nc(lat_reg_index);
CHMY(:,:,:)=CHM(lon_reg_index,lat_reg_index,:);
%% 根据等深线选取区域
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

% m_proj需要下载M_Map来用于绘制地图
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

% 时间序列化
abc(:,:,:)=nanmean(CHLL(:,:,:),2);
chll(:,:)=squeeze(nanmean(abc(:,:,:),1));
%% 画季节分布图
l=1
for k=1:11:14
    chl1(:,:,l)=nanmean(CHLL(:,:,k:(k+1)),3);%冬季
    chl2(:,:,l)=nanmean(CHLL(:,:,(k+2):(k+4)),3);%春季
    chl3(:,:,l)=nanmean(CHLL(:,:,(k+5):(k+7)),3);%秋季
    chl4(:,:,l)=nanmean(CHLL(:,:,(k+8):(k+10)),3);%夏季
    l=l+1;
end
%% 冬季 Winter
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
%% 春季 Spring
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
%% 秋季 Autumn
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
%% 夏季 Summer
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