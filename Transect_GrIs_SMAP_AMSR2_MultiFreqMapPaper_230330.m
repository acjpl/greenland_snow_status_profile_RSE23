close all; clear; 

global HOMEDIR;

%% This scripts plots the TB and air temperature for the K-transect in Greenland.
%
% Author: Andreas Colliander, Mohammad Mousavi, March 2023, NASA-JPL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
savemovie       = 0;
savefig         = 0;

driveletter     = 'Y';

SAVEFIG_PATH    = [driveletter, ':\IceSheets\matlab\Figs\GrIs\Transects230123\'];

%% Load SMAP data
%local directory for faster reading
datadir = [HOMEDIR, '_Temp\icesheetdata\']%'Y:\IceSheets\data\Greenland Data\'%Greenland_combined_2015_2020'

load([datadir,'Greenland_combined_2015_2020.mat'])

datadir1 = [driveletter, ':\IceSheets\data\'];
% datadir = [datadir1, '\Greenland Data\'];
% load([datadir,'SMAP_DATA_Greenland_2015_2019.mat']);

%
%% load AMSR2 data

% datadir_amsr2='Y:\AMSR2\DATA\MAT09KM_GrIs\'
datadir_amsr2 = [HOMEDIR, '_Temp\icesheetdata\'];
x_amsr2       = load([datadir_amsr2, 'AMSR2_L1B_GrIs_20150401_20191231.mat']);

case_id='pm';

if strcmpi('PM',case_id)

    case_select = 'PM';
    tbh_06      = x_amsr2.TBH_06_A;
    tbv_06      = x_amsr2.TBV_06_A;
    tbh_10      = x_amsr2.TBH_10_A;
    tbv_10      = x_amsr2.TBV_10_A;
    tbh_18      = x_amsr2.TBH_18_A;
    tbv_18      = x_amsr2.TBV_18_A;
    tbh_36      = x_amsr2.TBH_36_A;
    tbv_36      = x_amsr2.TBV_36_A;
    dnum_amsr2  = x_amsr2.dnum_A;

elseif strcmpi('AM',case_id)

    case_select = 'AM';
    tbh_06      = x_amsr2.TBH_06_D;
    tbv_06      = x_amsr2.TBV_06_D;
    tbh_10      = x_amsr2.TBH_10_D;
    tbv_10      = x_amsr2.TBV_10_D;
    tbh_18      = x_amsr2.TBH_18_D;
    tbv_18      = x_amsr2.TBV_18_D;
    tbh_36      = x_amsr2.TBH_36_D;
    tbv_36      = x_amsr2.TBV_36_D;
    dnum_amsr2  = x_amsr2.dnum_D;
end

%% Load ice edge
promicedir = [datadir1, '\Promice_ice_mask\'];
land       = shaperead('landareas.shp', 'UseGeoCoords', true);

% A_GrIS = 1716555; %km2 -- https://www.promice.dk/RetreatingIce.html;jsessionid=0AB4246BAECD765D889EE6659AB92686?cid=3816

if exist([promicedir, 'iceedge.mat'], 'file')
    load([promicedir, 'iceedge.mat']);
else
    
    [iceedge] = shaperead([promicedir, 'ice_mask.shp']);
    
    addpath([driveletter, ':\IceSheets\matlab\EPSG']);
    
    for i = 1:length(iceedge)
        [iceedge(i).Lat, iceedge(i).Lon] = polarstereo_inv(iceedge(i).X, iceedge(i).Y, 6378137.0, 0.08181919, 70, -45);
    end
    
    % The last struct seems to have the edges
    % Remove trailing nan from shapefile
    ice_lat = iceedge(end).Lat(1:end-1);
    ice_lon = iceedge(end).Lon(1:end-1);
    
    ice_mask = inpolygon(lat_area, lon_area, ice_lat, ice_lon);
    
    save([promicedir, 'iceedge.mat'], 'iceedge', 'ice_mask', 'ice_lat', 'ice_lon');
    
end
disp('Ice edge loaded');

ice_mask2 = double(ice_mask);
ice_mask2(ice_mask2==0) = NaN;

%%

for n = 1:length(dnum_pm(1,1,:))
    
    dnum_pm_x = dnum_pm(:,:,n);
    
    dnum_pm_m(n,1) = nanmean_mathworks(dnum_pm_x(:));
    
    dvec_pm_m = datevec(dnum_pm_m);
end

for n = 1:length(dnum_amsr2(1,1,:))
    dnum_amsr2_x = dnum_amsr2(:,:,n);
    dnum_amsr2_m(n,1) =(datenum(1993,1,1,0,0,0)+ nanmean(dnum_amsr2_x(:))./86400);
    dvec_amsr2_m = datevec(dnum_amsr2_m);
end

%% Read/scan the selected met data file

addpath([driveletter, ':\IceSheets\Weather Stations Data\Greenland\PROMICE']);

type                    = 'hour';
% type='day';
ver                     = 3;
stat_name               = 'KAN_U';
[dnum_KAN_U,Tair_KAN_U] = read_AWS_station(stat_name,ver,type, 'Y');
stat_name               = 'KAN_L';
[dnum_KAN_L,Tair_KAN_L] = read_AWS_station(stat_name,ver,type, 'Y');
stat_name               = 'KAN_M';
[dnum_KAN_M,Tair_KAN_M] = read_AWS_station(stat_name,ver,type, 'Y');
stat_name               = 'KAN_B';
[dnum_KAN_B,Tair_KAN_B] = read_AWS_station(stat_name,ver,type, 'Y');

%% Calculate the T_air statistics (min/max/avg)

dnums_KAN_U = floor(dnum_KAN_U);
dnumu_KAN_U = unique(floor(dnum_KAN_U));
dnums_KAN_L = floor(dnum_KAN_L);
dnumu_KAN_L = unique(floor(dnum_KAN_L));
dnums_KAN_M = floor(dnum_KAN_M);
dnumu_KAN_M = unique(floor(dnum_KAN_M));
dnums_KAN_B = floor(dnum_KAN_B);
dnumu_KAN_B = unique(floor(dnum_KAN_B));

for n = 1:length(dnumu_KAN_U)
    
    b = dnumu_KAN_U(n) == dnums_KAN_U;
    
    Tair_ave_KAN_U(n,1) = nanmean(Tair_KAN_U(b));
    Tair_max_KAN_U(n,1) = nanmax(Tair_KAN_U(b));
    Tair_min_KAN_U(n,1) = nanmin(Tair_KAN_U(b));
    
end
Tair_max_KAN_U_diff = gradient(Tair_max_KAN_U);

for n = 1:length(dnumu_KAN_L)
    
    b = dnumu_KAN_L(n) == dnums_KAN_L;
    
    Tair_ave_KAN_L(n,1) = nanmean(Tair_KAN_L(b));
    Tair_max_KAN_L(n,1) = nanmax(Tair_KAN_L(b));
    Tair_min_KAN_L(n,1) = nanmin(Tair_KAN_L(b));
    
end
Tair_max_KAN_L_diff = gradient(Tair_max_KAN_L);

for n = 1:length(dnumu_KAN_M)
    
    b = dnumu_KAN_M(n) == dnums_KAN_M;
    
    Tair_ave_KAN_M(n,1) = nanmean(Tair_KAN_M(b));
    Tair_max_KAN_M(n,1) = nanmax(Tair_KAN_M(b));
    Tair_min_KAN_M(n,1) = nanmin(Tair_KAN_M(b));
    
end
Tair_max_KAN_M_diff = gradient(Tair_max_KAN_M);

for n = 1:length(dnumu_KAN_B)
    
    b = dnumu_KAN_B(n) == dnums_KAN_B;
    
    Tair_ave_KAN_B(n,1) = nanmean(Tair_KAN_B(b));
    Tair_max_KAN_B(n,1) = nanmax(Tair_KAN_B(b));
    Tair_min_KAN_B(n,1) = nanmin(Tair_KAN_B(b));
    
end
Tair_max_KAN_B_diff = gradient(Tair_max_KAN_B);

%%
yrs             = 2015:2020;
yr_select       = 2019;
m               = find(yrs==yr_select);

dnum_1          = datenum(yrs(m), 06, 15);

ind_time        = find(floor(dnum_pm_m)== dnum_1);
ind_time_amsr2  = find(floor(dnum_amsr2_m)== dnum_1);
ind_time_KAN_U  = find(dnumu_KAN_U == dnum_1);
ind_time_KAN_L  = find(dnumu_KAN_L == dnum_1);
ind_time_KAN_M  = find(dnumu_KAN_M == dnum_1);
ind_time_KAN_B  = find(dnumu_KAN_B == dnum_1);

%%

% L=KAN_L
lat_KAN_L = 67.0955;
lon_KAN_L = -49.9513;

% KAN_M
lat_KAN_M = 67.0670;
lon_KAN_M = -48.8355;

% KAN_U
lat_KAN_U = 67.0003;
lon_KAN_U = -47.0253;

% KAN_B
lat_KAN_B = 67.1252;
lon_KAN_B = -50.1832;

lat_avg   = mean([lat_KAN_L, lat_KAN_M, lat_KAN_U, lat_KAN_B]);
%Greenland min and max r/c
minr      = 1040;
maxr      = 1280;
minc      = 700;
maxc      = 980;

lon_min   = lon_KAN_B-0.1;
lon_max   = -43;
lon       = lon_min:0.01:lon_max;

row       = NaN(length(lon),1);
col       = NaN(length(lon),1);

for ii=1:length(lon)
    
    [col1, row1] = get_grid_col_row_for_lat_lon(lat_avg, lon(ii), 9, 'N');
    
    row(ii)      = row1-minr+1;
    col(ii)      = col1-minc+1;
    
end

coords              = [row,col];
coords              = unique(coords, 'rows', 'stable');

dnum1               = datenum(yrs(m), 01, 01);
dnum2               = datenum(yrs(m), 12, 31);

dnum_ref1           = datenum(yrs(m), 04,  1);
dnum_ref2           = datenum(yrs(m), 04,  6);

inds_ref            = find((dnum_pm_m    > dnum_ref1 & dnum_pm_m    < dnum_ref2));
inds_ref_amsr2      = find((dnum_amsr2_m > dnum_ref1 & dnum_amsr2_m < dnum_ref2));

inds_smap_period    = find(dnum_pm_m > dnum1 & dnum_pm_m < dnum2);
inds_amsr2_period   = find(dnum_amsr2_m > dnum1 & dnum_amsr2_m < dnum2);

npr_pm_masked       = npr_pm;
tbv_pm_masked       = tbv_pm;
tbh_pm_masked       = tbh_pm;


tbv_est_wref        = nanmean_mathworks(tbv_pm_masked(:,:,inds_ref),3);
tbv_06_wref         = nanmean_mathworks(tbv_06(:,:,inds_ref_amsr2),3);
tbv_10_wref         = nanmean_mathworks(tbv_10(:,:,inds_ref_amsr2),3);
tbv_18_wref         = nanmean_mathworks(tbv_18(:,:,inds_ref_amsr2),3);
tbv_36_wref         = nanmean_mathworks(tbv_36(:,:,inds_ref_amsr2),3);

tbh_est_wref        = nanmean_mathworks(tbh_pm_masked(:,:,inds_ref),3);
tbh_06_wref         = nanmean_mathworks(tbh_06(:,:,inds_ref_amsr2),3);
tbh_10_wref         = nanmean_mathworks(tbh_10(:,:,inds_ref_amsr2),3);
tbh_18_wref         = nanmean_mathworks(tbh_18(:,:,inds_ref_amsr2),3);
tbh_36_wref         = nanmean_mathworks(tbh_36(:,:,inds_ref_amsr2),3);

tbv_wref_diff       = (tbv_pm_masked)- tbv_est_wref;
tbv_06_wref_diff    = (tbv_06)- tbv_06_wref;
tbv_10_wref_diff    = (tbv_10)- tbv_10_wref;
tbv_18_wref_diff    = (tbv_18)- tbv_18_wref;
tbv_36_wref_diff    = (tbv_36)- tbv_36_wref;

tbh_wref_diff       = (tbh_pm_masked)- tbh_est_wref;
tbh_06_wref_diff    = (tbh_06)- tbh_06_wref;
tbh_10_wref_diff    = (tbh_10)- tbh_10_wref;
tbh_18_wref_diff    = (tbh_18)- tbh_18_wref;
tbh_36_wref_diff    = (tbh_36)- tbh_36_wref;

[tbv_smap_absdiff_max, tbv_smap_absdiff_max_ind] ...
                    = nanmax(abs(tbv_wref_diff(:,:,inds_smap_period)),[], 3);

tbv_wref_diff_per   = tbv_wref_diff(:,:,inds_smap_period);
tbv_smap_maxdiff    = NaN*tbv_pm_masked(:,:,1);

for n1 = 1:length(tbv_smap_maxdiff(:,1))
    for n2 = 1:length(tbv_smap_maxdiff(1,:))
        tbv_smap_maxdiff(n1,n2) = tbv_wref_diff_per(n1,n2,tbv_smap_absdiff_max_ind(n1,n2));
    end
end

tbv_06_maxdiff      = nanmax(tbv_06_wref_diff,[], 3);
tbv_10_maxdiff      = nanmax(tbv_10_wref_diff,[], 3);
tbv_18_maxdiff      = nanmax(tbv_18_wref_diff,[], 3);
tbv_36_maxdiff      = nanmax(tbv_36_wref_diff,[], 3);

tbv_smap_scaled     = tbv_wref_diff./tbv_smap_maxdiff;
tbv_06_scaled       = tbv_06_wref_diff./tbv_06_maxdiff;
tbv_10_scaled       = tbv_10_wref_diff./tbv_10_maxdiff;
tbv_18_scaled       = tbv_18_wref_diff./tbv_18_maxdiff;
tbv_36_scaled       = tbv_36_wref_diff./tbv_36_maxdiff;
    

for ii=1:size(coords,1)
    
    tbv_tran_diff(ii)    = tbv_wref_diff(coords(ii,1),coords(ii,2),ind_time);    
    tbv_06_tran_diff(ii) = tbv_06_wref_diff(coords(ii,1),coords(ii,2),ind_time_amsr2);
    tbv_10_tran_diff(ii) = tbv_10_wref_diff(coords(ii,1),coords(ii,2),ind_time_amsr2);
    tbv_18_tran_diff(ii) = tbv_18_wref_diff(coords(ii,1),coords(ii,2),ind_time_amsr2);
    tbv_36_tran_diff(ii) = tbv_36_wref_diff(coords(ii,1),coords(ii,2),ind_time_amsr2);
    
    tbh_tran_diff(ii)    = tbh_wref_diff(coords(ii,1),coords(ii,2),ind_time);    
    tbh_06_tran_diff(ii) = tbh_06_wref_diff(coords(ii,1),coords(ii,2),ind_time_amsr2);
    tbh_10_tran_diff(ii) = tbh_10_wref_diff(coords(ii,1),coords(ii,2),ind_time_amsr2);
    tbh_18_tran_diff(ii) = tbh_18_wref_diff(coords(ii,1),coords(ii,2),ind_time_amsr2);
    tbh_36_tran_diff(ii) = tbh_36_wref_diff(coords(ii,1),coords(ii,2),ind_time_amsr2);
    
    tbv_smap_tran_scaled(ii) = tbv_smap_scaled(coords(ii,1),coords(ii,2),ind_time);    
    tbv_06_tran_scaled(ii) = tbv_06_scaled(coords(ii,1),coords(ii,2),ind_time_amsr2);
    tbv_10_tran_scaled(ii) = tbv_10_scaled(coords(ii,1),coords(ii,2),ind_time_amsr2);
    tbv_18_tran_scaled(ii) = tbv_18_scaled(coords(ii,1),coords(ii,2),ind_time_amsr2);
    tbv_36_tran_scaled(ii) = tbv_36_scaled(coords(ii,1),coords(ii,2),ind_time_amsr2);    
    
    lons(ii)             = lon_area(coords(ii,1),coords(ii,2));
    lats(ii)             = lat_area(coords(ii,1),coords(ii,2));
    ice_mask_edge(ii)    = ice_mask2(coords(ii,1),coords(ii,2));
end

%finidng the lon of the edge of the ice mask for the transect (west side)
%as the transect is on the west side.
ind_tmp           = find(isnan(ice_mask_edge));
ice_mask_lon_edge = lons(ind_tmp(end));

dt_tmp            = (cellstr(datetime(datevec(dnum_1),'Format','MMMM d, y')));

%% 

left_color  = [0 0 0];
right_color = [1 0 0];

%%

datesel = [...
    2019, 4, 25;
    2019, 5, 28;
    2019, 6,  9;
    2019, 8,  8;
    2019, 8, 10;
    2019, 8, 17;
    2019, 9, 25];

dnums = datenum(datesel);

% Verify 2017 late October melt
% dnums = datenum(2017,10,1):datenum(2017,11,15);

for wn = 1:length(dnums)
    ind_time2(wn) = find(floor(dnum_amsr2_m) == dnums(wn));
end



%% TBV

VarNameMap   = '\DeltaTBV';
VarName      = '\DeltaTBV';% 'Scaled \DeltaTBV';

savefilename = 'TBV_SMAP_AMSR2';

Y_lim_map    = [-100 100];
Y_lim        = [-100 100];%[-1 1.2];

for ii=1:length(ind_time2)
    
    ind_smap                        = find(floor(dnum_amsr2_m(ind_time2(ii))) == floor(dnum_pm_m))

    dt_tmp                          = (cellstr(datetime(datevec(dnum_amsr2_m(ind_time2(ii))),'Format','MMMM d, y')));

    % Use SMAP measurement time stamps to find the closest to in other stations
    
    [close_check_U,closestIndex_U]  = min(abs(dnumu_KAN_U-floor(dnum_amsr2_m(ind_time2(ii)))'));
    [close_check_L,closestIndex_L]  = min(abs(dnumu_KAN_L-floor(dnum_amsr2_m(ind_time2(ii)))'));
    [close_check_M,closestIndex_M]  = min(abs(dnumu_KAN_M-floor(dnum_amsr2_m(ind_time2(ii)))'));
    [close_check_B,closestIndex_B]  = min(abs(dnumu_KAN_B-floor(dnum_amsr2_m(ind_time2(ii)))'));
   
    
    if close_check_U==0
        bar_Tair_U = (Tair_max_KAN_U(closestIndex_U));
    else
        bar_Tair_U = NaN;
    end
    
    if close_check_L==0
        bar_Tair_L = (Tair_max_KAN_L(closestIndex_L));
    else
        bar_Tair_L = NaN;
    end
    
    if close_check_M==0
        bar_Tair_M = (Tair_max_KAN_M(closestIndex_M));
    else
        bar_Tair_M = NaN;
    end
    
    if close_check_B==0
        bar_Tair_B = (Tair_max_KAN_B(closestIndex_B));
    else
        bar_Tair_B = NaN;
    end


    for jj=1:size(coords,1)
        
        if isempty(ind_smap)
            tbv_tran_diff(jj)    = NaN;
        else
            tbv_tran_diff(jj)    = tbv_wref_diff(   coords(jj,1), coords(jj,2), ind_smap); 
        end
        tbv_06_tran_diff(jj) = tbv_06_wref_diff(coords(jj,1), coords(jj,2), ind_time2(ii));
        tbv_10_tran_diff(jj) = tbv_10_wref_diff(coords(jj,1), coords(jj,2), ind_time2(ii));
        tbv_18_tran_diff(jj) = tbv_18_wref_diff(coords(jj,1), coords(jj,2), ind_time2(ii));
        tbv_36_tran_diff(jj) = tbv_36_wref_diff(coords(jj,1), coords(jj,2), ind_time2(ii));

    end
    
    if isempty(ind_smap)
        MAP_L   = NaN*ice_mask2;    
    else
        MAP_L   = squeeze(tbv_wref_diff(   :,:,ind_smap)).*ice_mask2;    
    end
    MAP_6   = squeeze(tbv_06_wref_diff(:,:,ind_time2(ii))).*ice_mask2;
    MAP_10  = squeeze(tbv_10_wref_diff(:,:,ind_time2(ii))).*ice_mask2;
    MAP_18  = squeeze(tbv_18_wref_diff(:,:,ind_time2(ii))).*ice_mask2;
    MAP_36  = squeeze(tbv_36_wref_diff(:,:,ind_time2(ii))).*ice_mask2;
    
    TRAN_L  = tbv_tran_diff;
    TRAN_6  = tbv_06_tran_diff;
    TRAN_10 = tbv_10_tran_diff;
    TRAN_18 = tbv_18_tran_diff;
    TRAN_36 = tbv_36_tran_diff;
    
    BAR_L   = bar_Tair_L;
    BAR_M   = bar_Tair_M;
    BAR_U   = bar_Tair_U;
    BAR_B   = bar_Tair_B;
    
    % >>> plotting script
    close all


%  v v v v v v v v


    Transect_GrIs_AMSR2_AC_plotX_timeseries_MultiFreqMapPaper;
    

%  ^ ^ ^ ^ ^ ^ ^ ^       
   

    if savefig
        
        if ~exist(SAVEFIG_PATH, 'dir')
            mkdir(SAVEFIG_PATH);
        end
        
        saveas(gcf, sprintf('%sGrIs_%s_%s.png',...
            SAVEFIG_PATH, savefilename, datestr((dnum_amsr2_m(ind_time2(ii))), 'yyyymmdd')));
        
    end
    
end

%% TBH

VarNameMap   = '\DeltaTBH';
VarName      = '\DeltaTBH';% 'Scaled \DeltaTBV';

savefilename = 'TBH_SMAP_AMSR2';

Y_lim_map    = [-100 100];
Y_lim        = [-100 100];%[-1 1.2];

for ii=1:length(ind_time2)
    
    ind_smap                        = find(floor(dnum_amsr2_m(ind_time2(ii))) == floor(dnum_pm_m))

    dt_tmp                          = (cellstr(datetime(datevec(dnum_amsr2_m(ind_time2(ii))),'Format','MMMM d, y')));

    % Use SMAP measurement time stamps to find the closest to in other stations
    
    [close_check_U,closestIndex_U]  = min(abs(dnumu_KAN_U-floor(dnum_amsr2_m(ind_time2(ii)))'));
    [close_check_L,closestIndex_L]  = min(abs(dnumu_KAN_L-floor(dnum_amsr2_m(ind_time2(ii)))'));
    [close_check_M,closestIndex_M]  = min(abs(dnumu_KAN_M-floor(dnum_amsr2_m(ind_time2(ii)))'));
    [close_check_B,closestIndex_B]  = min(abs(dnumu_KAN_B-floor(dnum_amsr2_m(ind_time2(ii)))'));
   
    
    if close_check_U==0
        bar_Tair_U = (Tair_max_KAN_U(closestIndex_U));
    else
        bar_Tair_U = NaN;
    end
    
    if close_check_L==0
        bar_Tair_L = (Tair_max_KAN_L(closestIndex_L));
    else
        bar_Tair_L = NaN;
    end
    
    if close_check_M==0
        bar_Tair_M = (Tair_max_KAN_M(closestIndex_M));
    else
        bar_Tair_M = NaN;
    end
    
    if close_check_B==0
        bar_Tair_B = (Tair_max_KAN_B(closestIndex_B));
    else
        bar_Tair_B = NaN;
    end


    for jj=1:size(coords,1)
        
        if isempty(ind_smap)
            tbh_tran_diff(jj)    = NaN;
        else
            tbh_tran_diff(jj)    = tbh_wref_diff(   coords(jj,1), coords(jj,2), ind_smap); 
        end
        tbh_06_tran_diff(jj) = tbh_06_wref_diff(coords(jj,1), coords(jj,2), ind_time2(ii));
        tbh_10_tran_diff(jj) = tbh_10_wref_diff(coords(jj,1), coords(jj,2), ind_time2(ii));
        tbh_18_tran_diff(jj) = tbh_18_wref_diff(coords(jj,1), coords(jj,2), ind_time2(ii));
        tbh_36_tran_diff(jj) = tbh_36_wref_diff(coords(jj,1), coords(jj,2), ind_time2(ii));

    end
    
    if isempty(ind_smap)
        MAP_L   = NaN*ice_mask2;    
    else
        MAP_L   = squeeze(tbh_wref_diff(   :,:,ind_smap)).*ice_mask2;    
    end
    MAP_6   = squeeze(tbh_06_wref_diff(:,:,ind_time2(ii))).*ice_mask2;
    MAP_10  = squeeze(tbh_10_wref_diff(:,:,ind_time2(ii))).*ice_mask2;
    MAP_18  = squeeze(tbh_18_wref_diff(:,:,ind_time2(ii))).*ice_mask2;
    MAP_36  = squeeze(tbh_36_wref_diff(:,:,ind_time2(ii))).*ice_mask2;
    
    TRAN_L  = tbh_tran_diff;
    TRAN_6  = tbh_06_tran_diff;
    TRAN_10 = tbh_10_tran_diff;
    TRAN_18 = tbh_18_tran_diff;
    TRAN_36 = tbh_36_tran_diff;
    
    BAR_L   = bar_Tair_L;
    BAR_M   = bar_Tair_M;
    BAR_U   = bar_Tair_U;
    BAR_B   = bar_Tair_B;
    

    % >>> plotting script
    close all

%  v v v v v v v v


    Transect_GrIs_AMSR2_AC_plotX_timeseries_MultiFreqMapPaper;
    

%  ^ ^ ^ ^ ^ ^ ^ ^   


    if savefig
        
        if ~exist(SAVEFIG_PATH, 'dir')
            mkdir(SAVEFIG_PATH);
        end
        
        saveas(gcf, sprintf('%sGrIs_%s_%s.png',...
            SAVEFIG_PATH, savefilename, datestr((dnum_amsr2_m(ind_time2(ii))), 'yyyymmdd')));
        
    end
    
end

