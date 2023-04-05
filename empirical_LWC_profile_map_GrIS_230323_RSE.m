close all; clear;
%% this is the cleanup branch version for now

%%
global HOMEDIR

%% setup some parameters
%
savedailyfigs         = 0;

savesnowstatusprofile = 0; % load multiyear set

savemovie             = 0;
savefigyes            = 0;

driveletter           = 'Y';

curdate               = datestr(now, 'yymmdd');

%%

datadir1 = [driveletter, ':\IceSheets\data\'];

%% Load SMAP and AMSR2 data

% FROM: HOMEDIR\11 - _PROJECTS_\07 - MEASURES FT\04 - Analysis\06 - Multi-freq\
% FROM: combine_and_save_SMAP_and_AMSR2_over_GrIS_211109.m
% fname = 'SMAP_and_AMSR2_GrIS_2016_20160101_20161231_v220725.mat'
% fname = 'SMAP_and_AMSR2_GrIS_2017_20170101_20171231_v220725.mat'
% fname = 'SMAP_and_AMSR2_GrIS_2018_20180101_20181231_v220725.mat'
fname = '../SMAP_and_AMSR2_GrIS_2019_20190101_20191231_v220725.mat';
% FROM: combine_and_save_SMAP_and_AMSR2_over_GrIS_230323.m
% fname = '../SMAP_and_AMSR2_GrIS_20150401_20191231_v230323.mat';

x = load(fname);
%         AMorPM: [732×1 single]
%            TBH: {1×5 cell}
%            TBV: {1×5 cell}
%     dnum_AMSR2: [241×281×732 single]
%      dnum_SMAP: [241×281×732 single]
%          dnums: [732×1 single]
%       lat_area: [241×281 double]
%       lon_area: [241×281 double]
%           maxc: 980
%           maxr: 1280
%           minc: 700
%           minr: 1040

% remove funny data
dnum_funny = datenum([2016, 1, 16; 2017, 6, 19; 2017, 9, 6; 2017, 9, 8; 2018, 6, 22; 2019, 8, 12]);
ap_funny   = [0, 0, 0, 1, 0, 0];

indf = [];

k = 0;
for uu = 1:length(dnum_funny)
    
    ixi = find(x.dnums == dnum_funny(uu) & x.AMorPM == ap_funny(uu));
    
    if ~isempty(ixi)
        k = k + 1;
        indf(k) = ixi;
    end
    
end

if ~isempty(indf)
    x.AMorPM(indf) = [];
    for vv = 1:5
        x.TBH{vv}(:,:,indf) = [];
        x.TBV{vv}(:,:,indf) = [];
    end
    x.dnum_AMSR2(:,:,indf) = [];
    x.dnum_SMAP(:,:,indf)  = [];
    x.dnums(indf)          = [];
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
    
    %     save([promicedir, 'iceedge.mat'], 'iceedge', 'ice_mask', 'ice_lat', 'ice_lon');
    
end
disp('Ice edge loaded');

ice_mask2               = double(ice_mask);
ice_mask2(ice_mask2==0) = NaN;

%% Manipulate the TB data
[nr,nc,n3] = size(x.TBV{1});

for w = 1:5
    
    TBv(:,:,w)        = reshape(x.TBV{w}.*ice_mask2, nr*nc, n3);
    TBh(:,:,w)        = reshape(x.TBH{w}.*ice_mask2, nr*nc, n3);
    
    dnum_amsr2(:,:,w) = reshape(x.dnum_AMSR2, nr*nc, n3);
    dnum_smap(:,:,w)  = reshape(x.dnum_SMAP,  nr*nc, n3);
    
    TBv_max(:,1,w)    = zeros(nr*nc,1) + nanmax(nanmax(TBv(:,:,w)));
    TBh_max(:,1,w)    = zeros(nr*nc,1) + nanmax(nanmax(TBh(:,:,w)));
    
end

dnums                = double(x.dnums);

[yr, mm, dd, hr, mn] = datevec(dnums);

yrs                  = unique(yr);

%% Find transect grid indices
lon_min   = -50;
lon_max   = -43;
lon       = lon_min:0.01:lon_max;

rowx      = NaN(length(lon),1);
colx      = NaN(length(lon),1);

for ii=1:length(lon)
    [col1, row1]  = get_grid_col_row_for_lat_lon(67, lon(ii), 9, 'N');
    rowx(ii)      = row1-x.minr+1;
    colx(ii)      = col1-x.minc+1;
end

tran_colrow = [rowx,colx];
tran_colrow = unique(tran_colrow, 'rows', 'stable');

tran_inds   = sub2ind(size(x.lat_area), tran_colrow(:,1), tran_colrow(:,2));

lons_area   = x.lon_area(:); 
lon_trans   = lons_area(tran_inds);

%%
for n = 1:length(yrs)
    
    dn1             = datenum(yrs(n),  1,  1);
    dn2             = datenum(yrs(n), 12, 31);
    
    fr1             = datenum(yrs(n),  1,  1);
    fr2             = datenum(yrs(n),  1, 31);
    fr3             = datenum(yrs(n), 12,  1);
    fr4             = datenum(yrs(n), 12, 31);
    
    mr1             = datenum(yrs(n), 5,  1);
    mr2             = datenum(yrs(n), 9, 30);
        
    mri             = dnums >= mr1 & dnums <= mr2;
    
    dni             = dnums >= dn1 & dnums <= dn2;

    nz              = sum(dni);
    
    TBV_n           = TBv(:,dni,:);
    TBH_n           = TBh(:,dni,:);
    
    AMorPM_n{n}     = x.AMorPM(dni);

    dnums_n{n}      = dnums(dni);

    seasoni_n       = dnums_n{n} >= fr2 & dnums_n{n} <= fr3;

    fri             = (dnums_n{n} >= fr1 & dnums_n{n} <= fr2) | (dnums_n{n} >= fr3 & dnums_n{n} <= fr4);

    friApril        = (dnums_n{n} >= datenum(yrs(n),  4,  1) & dnums_n{n} <= datenum(yrs(n),  4,  6));%

    TBVr_F{n}       = nanmean(TBV_n(:, fri, :), 2);
    TBHr_F{n}       = nanmean(TBH_n(:, fri, :), 2);
    
    TBVr_F_april{n} = nanmean(TBV_n(:, friApril, :), 2);
    TBHr_F_april{n} = nanmean(TBH_n(:, friApril, :), 2);    

    SD_TBV{n}       = nanstd(TBV_n(:,fri,:), [], 2);
    SD_TBH{n}       = nanstd(TBH_n(:,fri,:), [], 2);    

    SD_TBVall{n}    = nanstd(TBV_n(:,:,:), [], 2);
    SD_TBHall{n}    = nanstd(TBH_n(:,:,:), [], 2);  
    
    min_TBVall{n}   = nanmin(TBV_n(:,:,:), [], 2);
    min_TBHall{n}   = nanmin(TBH_n(:,:,:), [], 2);  

    max_TBVall{n}   = nanmax(TBV_n(:,:,:), [], 2);
    max_TBHall{n}   = nanmax(TBH_n(:,:,:), [], 2);  
    
    mm_TBVall{n}    = max_TBVall{n} - min_TBVall{n};
    mm_TBHall{n}    = max_TBHall{n} - min_TBHall{n};
    
    E_SD_TBV{n}     = squeeze(nanmean(SD_TBV{n}));
    E_SD_TBH{n}     = squeeze(nanmean(SD_TBH{n}));  

    dTBV{n}         = TBV_n - TBVr_F{n};
    dTBH{n}         = TBH_n - TBHr_F{n};
    
    TBVr_maxd{n}    = NaN(nr*nc, 1, 5);
    TBHr_maxd{n}    = NaN(nr*nc, 1, 5);
    
    % Find max deviation and TB at max deviation
    for w = 1:5
        [dTBvmax(:,:,w), dTBvmax_inds(:,1,w)] = nanmax(abs(dTBV{n}(:,seasoni_n,w)), [], 2);
        [dTBhmax(:,:,w), dTBhmax_inds(:,1,w)] = nanmax(abs(dTBH{n}(:,seasoni_n,w)), [], 2);
        
        indsv               = (dTBvmax_inds(:,1,w)-1)*nr*nc + [1:nr*nc]';
        indsh               = (dTBhmax_inds(:,1,w)-1)*nr*nc + [1:nr*nc]';
        
        TBVx                = TBV_n(:,seasoni_n,w);
        TBHx                = TBH_n(:,seasoni_n,w);
        
        TBVr_maxd{n}(:,1,w) = TBVx(indsv);
        TBHr_maxd{n}(:,1,w) = TBHx(indsv);
    end
    % Find L-band negative devations area    
    negdevL_V   = TBVr_maxd{n}(:,1,1)-TBVr_F{n}(:,1,1);
    negdevL_H   = TBHr_maxd{n}(:,1,1)-TBHr_F{n}(:,1,1);
    
    negdevL_V_b = negdevL_V < -5*E_SD_TBV{n}(1);
    negdevL_H_b = negdevL_H < -5*E_SD_TBH{n}(1);
    
    % Melt reference
    % == 1.4 GHz - for positive and negative change
    TBVr_M{n}(:, 1,   1)         = 260 + zeros(nr*nc, 1, 1); % TBv_max;%
    TBHr_M{n}(:, 1,   1)         = 250 + zeros(nr*nc, 1, 1); % TBh_max;%
   
    TBVr_M{n}(negdevL_V_b, 1, 1) = TBVr_maxd{n}(negdevL_V_b,1,1);
    TBHr_M{n}(negdevL_H_b, 1, 1) = TBHr_maxd{n}(negdevL_H_b,1,1);
    
    % == 6.9-36.5 GHz
    TBVr_M{n}(:, 1, 2:5)         = 273 + zeros(nr*nc, 1, 4); % TBv_max;%
    TBHr_M{n}(:, 1, 2:5)         = 273 + zeros(nr*nc, 1, 4); % TBh_max;%
    
    % Melt factor
    SFV{n}             = (TBV_n-TBVr_F{n})./(TBVr_M{n}-TBVr_F{n});
    SFH{n}             = (TBH_n-TBHr_F{n})./(TBHr_M{n}-TBHr_F{n});
    
    SD_SFV{n}          = nanstd(SFV{n}(:,fri,:), [], 2);
    SD_SFH{n}          = nanstd(SFH{n}(:,fri,:), [], 2);

    E_SD_SFV{n}        = squeeze(nanmean(SD_SFV{n}));
    E_SD_SFH{n}        = squeeze(nanmean(SD_SFH{n}));
    
    SFV{n}(SFV{n} < 0) = 0;
    SFH{n}(SFH{n} < 0) = 0;
    
    Zv   = [6,6,6,6,6];
    for w = 1:5
        MbV{n}(:,:,w) = SFV{n}(:,:,w) > Zv(w)*E_SD_SFV{n}(w);
        MbH{n}(:,:,w) = SFH{n}(:,:,w) > Zv(w)*E_SD_SFH{n}(w);
    end
    
    [Nr, Nc, Zz] = size(MbV{n});

    % classification
    lyrsV{1}  = MbV{n}(:,:,5)       == 1     & all(MbV{n}(:,:,1:4) == 0, 3);                                                   % 1d1: 1st layer M, others F
    lyrsV{2}  = all(MbV{n}(:,:,4:5) == 1, 3) & all(MbV{n}(:,:,1:3) == 0, 3);                                                   % 1d2: the first two layers M, others F
    lyrsV{3}  = MbV{n}(:,:,5)       == 1     & MbV{n}(:,:,3)       == 1     & all(MbV{n}(:,:,1:2) == 0, 3);                    % 1d3: 1st and 3rd layer M; 4th and 5th F; 2nd don't care
    lyrsV{4}  = MbV{n}(:,:,5)       == 1     & MbV{n}(:,:,2)       == 1     & MbV{n}(:,:,1)       == 0;                        % 1d4: 1st and 4th layer M; 5th F; 2nd and 3rd don't care
    lyrsV{5}  = MbV{n}(:,:,5)       == 1     & MbV{n}(:,:,1)       == 1;                                                       % 1d5: 1st and 4th layer M; 2nd, 3rd, 4th don't care
    lyrsV{6}  = MbV{n}(:,:,5)       == 0     & MbV{n}(:,:,4)       == 1     & MbV{n}(:,:,2)       == 1 & MbV{n}(:,:,1)   == 0; % 2d4
    lyrsV{7}  = all(MbV{n}(:,:,4:5) == 0, 3) & all(MbV{n}(:,:,2:3) == 1, 3) & MbV{n}(:,:,1)       == 0;                        % 3d4
    lyrsV{8}  = MbV{n}(:,:,5)       == 0     & MbV{n}(:,:,4)       == 1     & MbV{n}(:,:,1)       == 1;                        % 2d5:
    lyrsV{9}  = all(MbV{n}(:,:,4:5) == 0, 3) & MbV{n}(:,:,3)       == 1     & MbV{n}(:,:,1)       == 1;                        % 3d5: 
    lyrsV{10} = all(MbV{n}(:,:,3:5) == 0, 3) & MbV{n}(:,:,2)       == 1     & MbV{n}(:,:,1)       == 1;                        % 4d5:
    lyrsV{11} = all(MbV{n}(:,:,2:5) == 0, 3) & MbV{n}(:,:,1)       == 1;                                                       % 5d5: 
    
    lyrsH{1}  = MbH{n}(:,:,5)       == 1     & all(MbH{n}(:,:,1:4) == 0, 3);                                                   % 1d1: 1st layer M, others F
    lyrsH{2}  = all(MbH{n}(:,:,4:5) == 1, 3) & all(MbH{n}(:,:,1:3) == 0, 3);                                                   % 1d2: the first two layers M, others F
    lyrsH{3}  = MbH{n}(:,:,5)       == 1     & MbH{n}(:,:,3)       == 1     & all(MbH{n}(:,:,1:2) == 0, 3);                    % 1d3: 1st and 3rd layer M; 4th and 5th F; 2nd don't care
    lyrsH{4}  = MbH{n}(:,:,5)       == 1     & MbH{n}(:,:,2)       == 1     & MbH{n}(:,:,1)       == 0;                        % 1d4: 1st and 4th layer M; 5th F; 2nd and 3rd don't care
    lyrsH{5}  = MbH{n}(:,:,5)       == 1     & MbH{n}(:,:,1)       == 1;                                                       % 1d5: 1st and 4th layer M; 2nd, 3rd, 4th don't care
    lyrsH{6}  = MbH{n}(:,:,5)       == 0     & MbH{n}(:,:,4)       == 1     & MbH{n}(:,:,2)       == 1 & MbH{n}(:,:,1)   == 0; % 2d4
    lyrsH{7}  = all(MbH{n}(:,:,4:5) == 0, 3) & all(MbH{n}(:,:,2:3) == 1, 3) & MbH{n}(:,:,1)       == 0;                        % 3d4
    lyrsH{8}  = MbH{n}(:,:,5)       == 0     & MbH{n}(:,:,4)       == 1     & MbH{n}(:,:,1)       == 1;                        % 2d5:
    lyrsH{9}  = all(MbH{n}(:,:,4:5) == 0, 3) & MbH{n}(:,:,3)       == 1     & MbH{n}(:,:,1)       == 1;                        % 3d5: 
    lyrsH{10} = all(MbH{n}(:,:,3:5) == 0, 3) & MbH{n}(:,:,2)       == 1     & MbH{n}(:,:,1)       == 1;                        % 4d5:
    lyrsH{11} = all(MbH{n}(:,:,2:5) == 0, 3) & MbH{n}(:,:,1)       == 1;                                                       % 5d5: 
    
    MSP_V{n}  = NaN(Nr, Nc);
    MSP_H{n}  = NaN(Nr, Nc);
    
    for q = 1:length(lyrsV)
        MSP_V{n}(lyrsV{q}) = q;
        MSP_H{n}(lyrsH{q}) = q;
    end    
    
    % difference between frequencies
    subM_KaC_V{n}      = SFV{n}(:,:,5) - SFV{n}(:,:,2);
    subM_KaC_H{n}      = SFH{n}(:,:,5) - SFH{n}(:,:,2);
    
    % binary state for the difference
    subMd_KaC_V{n}     = subM_KaC_V{n} < -0.2;
    subMd_KaC_H{n}     = subM_KaC_H{n} < -0.2;
    
    % "number of days" for different cases
    NsubMd_KuC_V_AM{n} = nansum(subMd_KaC_V{n}(:,AMorPM_n{n} == 0), 2);
    NsubMd_KuC_V_PM{n} = nansum(subMd_KaC_V{n}(:,AMorPM_n{n} == 1), 2);
    NsubMd_KuC_H_AM{n} = nansum(subMd_KaC_H{n}(:,AMorPM_n{n} == 0), 2);
    NsubMd_KuC_H_PM{n} = nansum(subMd_KaC_H{n}(:,AMorPM_n{n} == 1), 2);
    
    NMbV_AM{n}         = nansum(MbV{n}(:,AMorPM_n{n} == 0, :), 2);
    NMbV_PM{n}         = nansum(MbV{n}(:,AMorPM_n{n} == 1, :), 2);
    NMbH_AM{n}         = nansum(MbH{n}(:,AMorPM_n{n} == 0, :), 2);
    NMbH_PM{n}         = nansum(MbH{n}(:,AMorPM_n{n} == 1, :), 2);
    
    % freezup date
    FdateV{n} = NaN(Nr, 1, 5);
    for ii = 1:Nr
        for jj = 1:5
            b_all = MbV{1}(ii, :, jj);
            b_all(fri) = 0;
            indall = find(b_all);
            if ~isempty(indall)
                indlast = indall(end);
                FdateV{n}(ii,1,jj) = dnums(indlast)-datenum(yrs(n), 1, 1);
            end
        end
    end
    
    %
    
    % reshape matrizes for plotting
    for w = 1:5
        
        TBVp{n}{w}          = reshape(TBV_n(:,:,w), nr, nc, nz);
        TBHp{n}{w}          = reshape(TBH_n(:,:,w), nr, nc, nz);
        
        SFVp{n}{w}          = reshape(SFV{n}(:,:,w), nr, nc, nz);
        SFHp{n}{w}          = reshape(SFH{n}(:,:,w), nr, nc, nz);
        
        TBVr_Mp{n}{w}       = reshape(TBVr_M{n}(:,:,w), nr, nc, 1);
        TBHr_Mp{n}{w}       = reshape(TBHr_M{n}(:,:,w), nr, nc, 1);

        TBVr_Fp{n}{w}       = reshape(TBVr_F{n}(:,:,w), nr, nc, 1);
        TBHr_Fp{n}{w}       = reshape(TBHr_F{n}(:,:,w), nr, nc, 1);
        
        TBVr_F_aprilp{n}{w} = reshape(TBVr_F_april{n}(:,:,w), nr, nc, 1);
        TBHr_F_aprilp{n}{w} = reshape(TBHr_F_april{n}(:,:,w), nr, nc, 1);        

        Fdatep{n}{w}        = reshape(FdateV{n}(:,:,w), nr, nc, 1);
        Fdatep{n}{w}        = reshape(FdateV{n}(:,:,w), nr, nc, 1);
        
        SD_TBVp{n}{w}       = reshape(SD_TBV{n}(:,:,w), nr, nc, 1);
        SD_TBHp{n}{w}       = reshape(SD_TBH{n}(:,:,w), nr, nc, 1);        
        
        SD_TBVallp{n}{w}    = reshape(SD_TBVall{n}(:,:,w), nr, nc, 1);
        SD_TBHallp{n}{w}    = reshape(SD_TBHall{n}(:,:,w), nr, nc, 1);    
        
        mm_TBVallp{n}{w}    = reshape(mm_TBVall{n}(:,:,w), nr, nc, 1);
        mm_TBHallp{n}{w}    = reshape(mm_TBHall{n}(:,:,w), nr, nc, 1);            
        
        MbHp{n}{w}          = reshape(MbH{n}(:,:,w), nr, nc, nz);
        MbVp{n}{w}          = reshape(MbV{n}(:,:,w), nr, nc, nz);
        
        NMbH_AMp{n}{w}      = reshape(NMbH_AM{n}(:,1,w), nr, nc);
        NMbH_PMp{n}{w}      = reshape(NMbH_PM{n}(:,1,w), nr, nc);
        
        NMbV_AMp{n}{w}      = reshape(NMbV_AM{n}(:,1,w), nr, nc);
        NMbV_PMp{n}{w}      = reshape(NMbV_PM{n}(:,1,w), nr, nc);        
    end
    
    MSP_Vp{n}           = reshape(MSP_V{n}(:,:), nr, nc, nz);
    MSP_Hp{n}           = reshape(MSP_H{n}(:,:), nr, nc, nz);
    
    subM_KuC_Vp{n}      = reshape(subM_KaC_V{n}(:,:), nr, nc, nz);
    subM_KuC_Hp{n}      = reshape(subM_KaC_H{n}(:,:), nr, nc, nz);
    
    NsubMd_KuC_H_AMp{n} = reshape(NsubMd_KuC_H_AM{n}, nr, nc);
    NsubMd_KuC_H_PMp{n} = reshape(NsubMd_KuC_H_PM{n}, nr, nc);
    
    negdevL_Vp          = reshape(double(negdevL_V), nr, nc);
    negdevL_Hp          = reshape(double(negdevL_H), nr, nc);
    
    negdevL_V_bp        = reshape(double(negdevL_V_b), nr, nc);
    negdevL_H_bp        = reshape(double(negdevL_H_b), nr, nc);    

end

%% SAVE SNOW STATUS PROFILE (var MSP) INTO A FILE
if savesnowstatusprofile

    k = 0;

    for n = 1:length(yrs)

        xi1 = find(datenum(yrs(n),  4, 15) == dnums_n{n})
        xi2 = find(datenum(yrs(n), 11, 30) == dnums_n{n})

        for seld = xi1(1):xi2(end) %1:nz

            k = k + 1;

            if AMorPM_n{n}(seld) == 0
                AMPM(k,:) = 'AM';
            elseif AMorPM_n{n}(seld) == 1
                AMPM(k,:) = 'PM';
            end

            datnum(k,1) = dnums_n{n}(seld);

            for w = 1:5
                snowstatus_V(:,:,k,w) = MbVp{n}{w}(:,:,seld).*ice_mask2;
                snowstatus_H(:,:,k,w) = MbHp{n}{w}(:,:,seld).*ice_mask2;
            end

        end

        LAT = x.lat_area;
        LON = x.lon_area;

    end

    save('greenland_snow_status_profile_index_v1.mat', ...
        'snowstatus_V', 'snowstatus_H', ...
        'datnum', 'AMPM', 'LAT', 'LON', '-v7.3');

end

%% PLOTTING

addpath([HOMEDIR, '00 - General\matlab\funcs_ext\Arctic Mapping Tools v1.04']);

%

freqs = {'1.4 GHz', '6.9 GHz', '10.7 GHz', '18.9 GHz', '36.5 GHz'};

fntsz = 14;

%%
cmap_blue = rgb2rgb_colormap([255,255,255], [0, 0, 255], 20)/255;

figure; set(gcf, 'position', [19         686        1447         638]);
for w = 1:5
    axs(w) = subplot(1,5,w);
    set(axs(w), 'fontsize', 14);
    greenland('color', colo('grey1'));
    pcolorpsn(x.lat_area, x.lon_area, Fdatep{1}{w});
    plotpsn(ice_lat, ice_lon, 'color', colo('k0'), 'linewidth', 1);
    set(axs(w), 'clim', [180 300]);
    title(freqs{w}, 'fontsize', 14)
    colormap(cmap_blue);
    xlim([-0.7e6 0.9e6]);
    ylim([-3.5e6 -0.5e6 ]);
    if w == 5; colorbar; end
    set(axs(w), 'position', [0.03+(w-1)*0.18, 0.1, 0.20, 0.75]);
    if w > 1; set(axs(w), 'yticklabel', []); end
    xlabel('Easting [m]');
    if w == 1; ylabel('Northing [m]'); end
end
text(axs(3), 0, -200000, sprintf('Refreeze Timing (V-pol): %.0f', yrs), 'fontsize', 16, 'horizontalalignment', 'center');

% saveas(gcf, ['figs220901_refreeze_timing\', sprintf('refreeze_timing_V_%.0f_v%s.png', yrs, curdate) ]);
%%

cmap_red = rgb2rgb_colormap([255,0,0], [255, 255, 255], 20)/255;

cmap_r2b = [cmap_red; cmap_blue];

figure; set(gcf, 'position', [138   304   744   638]);
for w = 1:2
    axs(w) = subplot(1,2,w);
    set(axs(w), 'fontsize', fntsz);
    greenland('color', colo('grey1'));
    if w == 1
        pcolorpsn(x.lat_area, x.lon_area, Fdatep{1}{1}-Fdatep{1}{5});
        title('1.4-36.5 GHz', 'fontsize', fntsz);
    elseif w == 2
        pcolorpsn(x.lat_area, x.lon_area, Fdatep{1}{2}-Fdatep{1}{5});
        title('6.9-36.5 GHz', 'fontsize', fntsz)
    end
    plotpsn(ice_lat, ice_lon, 'color', colo('k0'), 'linewidth', 1);
    set(gca, 'clim', [-50 50]);
    colormap(cmap_r2b);
    xlim([-0.7e6  0.9e6]);
    ylim([-3.5e6 -0.5e6]);
    if w == 2
        cb = colorbar;
        set(cb, 'fontsize', fntsz);
    end
    set(axs(w), 'position', [0.03+(w-1)*0.4, 0.1, 0.50, 0.75]);
    if w > 1
        set(axs(w), 'yticklabel', []);
    else
        ylabel('Northing [km]');
    end
    xlabel('Easting [km]');
end

text(axs(1), 12e5, -0.2e6, sprintf('Difference in Refreeze Timing (V-pol): %.0f', yrs(1)), 'fontsize', 17, 'horizontalalignment', 'center');

% saveas(gcf, ['figs220901_refreeze_timing\', sprintf('refreeze_timing_difference_V_%.0f_v%s.png', yrs, curdate) ]);

%% Winter reference (April)

figure; set(gcf, 'position', [19         134        1447        1190]);%[19         686        1447         638]);
for w = 1:5
    axs(w) = subplot(2,5,w);
    greenland('color', colo('grey1'));
    pcolorpsn(x.lat_area, x.lon_area, TBVr_F_aprilp{n}{w}.*ice_mask2);
    plotpsn(ice_lat, ice_lon, 'color', colo('k0'), 'linewidth', 1);
    plotpsn([67,67], [-50,-43], 'color', colo('k0'), 'linewidth', 1, 'linewidth', 3);
    set(axs(w), 'clim', [100 270]);
    title(freqs{w}, 'fontsize', fntsz+1)
    xlim([-0.7e6 0.9e6]);
    ylim([-3.5e6 -0.5e6 ]);
    colormap(turbo);
    if w == 5
        cb1 = colorbar;
        set(cb1, 'fontsize', 14);
    end
    
    set(axs(w), 'position', [0.03+(w-1)*0.18, 0.55, 0.18, 0.4]);
    
    if w > 1
        set(axs(w), 'yticklabel', []);
    end
    
end

for w = 6:10
    axs(w) = subplot(2,5,w);
    greenland('color', colo('grey1'));
    pcolorpsn(x.lat_area, x.lon_area, TBHr_F_aprilp{n}{w-5}.*ice_mask2);
    plotpsn(ice_lat, ice_lon, 'color', colo('k0'), 'linewidth', 1);
    plotpsn([67,67], [-50,-43], 'color', colo('k0'), 'linewidth', 1, 'linewidth', 3);
    set(axs(w), 'clim', [100 270]);
    title(freqs{w-5}, 'fontsize', fntsz+1)
    xlim([-0.7e6 0.9e6]);
    ylim([-3.5e6 -0.5e6 ]);
    colormap(turbo);
    if w == 10
        cb2 = colorbar;
        set(cb2, 'fontsize', 14);
    end
    
    set(axs(w), 'position', [0.03+(w-5-1)*0.18, 0.05, 0.18, 0.40]);
    
    if w > 6
        set(axs(w), 'yticklabel', []);
    end
    
end

text(axs(3), 0, -200000, sprintf('Brighntess Temperature V-pol (April 1-6)'), 'fontsize', 17, 'horizontalalignment', 'center');
text(axs(8), 0, -200000, sprintf('Brighntess Temperature H-pol (April 1-6)'), 'fontsize', 17, 'horizontalalignment', 'center');

% saveas(gcf, ['figs220825_winter_TB\', sprintf('TB_April_%.0f_v%s.png', yrs, curdate) ]);

%%

XTBV  = double(squeeze(TBVr_F_april{n}(tran_inds,1,:))');
XTBH  = double(squeeze(TBHr_F_april{n}(tran_inds,1,:))');


figure; hold all
set(gcf, 'position', [489         275        1031         584]);

ax1 = subplot(2,1,1); hold all;
set(ax1, 'fontsize', 18);
L1 = plot(lon_trans, flipud(XTBV), 'linewidth', 3);
lg = legend(fliplr(freqs), 'location', 'northeastoutside');
set(ax1, 'xticklabel', []);
title('\rmBrightness Temperature Along Focus Transect (April 1-6)');
ylabel('TBV [K]');
ylim([100, 270]);
grid on; grid minor; zoom on;

ax2 = subplot(2,1,2); hold all;
set(ax2, 'fontsize', 18);
L1 = plot(lon_trans, flipud(XTBH), 'linewidth', 3);
ylabel('TBH [K]');
xlabel('Longitude [^o]');
ylim([100, 270]);
grid on; grid minor; zoom on;

set(ax1, 'position', [0.1, 0.55, 0.75, 0.35]);
set(ax2, 'position', [0.1, 0.18, 0.75, 0.35]);

% saveas(gcf, ['figs220825_winter_TB\', sprintf('TB_April_Transect_%.0f_v%s.png', yrs, curdate) ]);

%% Ku vs C: Number of Days with Subsurface Meltwater and Frozen Surface

cmapsn = colormap(parula);
cmapsn(1,:) = [1,1,1];

figure; set(gcf, 'position', [138   304   744   638]);
for w = 1:2
    axs(w) = subplot(1,2,w);
    set(axs(w), 'fontsize', fntsz);
    greenland('color', colo('grey1'), 'km');
    if w == 1
        pcolorpsn(x.lat_area, x.lon_area, NsubMd_KuC_H_AMp{n}.*ice_mask2, 'km');
        title('\rmAM', 'fontsize', fntsz+1);
    elseif w == 2
        pcolorpsn(x.lat_area, x.lon_area, NsubMd_KuC_H_PMp{n}.*ice_mask2, 'km');
        title('\rmPM', 'fontsize', fntsz+1)
    end
    plotpsn(ice_lat, ice_lon, 'color', colo('k0'), 'linewidth', 1, 'km');
    colormap(cmapsn);
    set(axs(w), 'clim', [0 100]);
    xlim([-0.7e3 0.9e3]);
    ylim([-3.5e3 -0.5e3 ]);
    
    if w == 2
        cb = colorbar;
        set(cb, 'fontsize', fntsz);
    end
    
    set(axs(w), 'position', [0.03+(w-1)*0.4, 0.1, 0.50, 0.75]);
    
    if w > 1
        set(axs(w), 'yticklabel', []);
    else
        ylabel('Northing [km]');
    end
    
    xlabel('Easting [km]');
    
end

text(axs(1), 12e2, -0.1e3, sprintf('Number of Days with Subsurface Melt (C) but Frozen Surface (Ku) (%.0f)', yrs(1)), 'fontsize', 16, 'horizontalalignment', 'center');

% saveas(gcf, ['figs220901_number_of_melt_days\', sprintf('Number_of_submelt_with_frozen_surf_KuC_%.0f_v%s.png', yrs(1), curdate ) ]);


%% Number of Ka V-pol melt days
cmapn = colormap(spring);
cmapn = flipud(cmapn);
cmapn(1,:) = [1,1,1];

figure; set(gcf, 'position', [138   304   744   638]);
for w = 1:2
    axs(w) = subplot(1,2,w);
    set(axs(w), 'fontsize', fntsz);
    greenland('color', colo('grey1'), 'km');
    if w == 1
        pcolorpsn(x.lat_area, x.lon_area, NMbV_AMp{n}{5}.*ice_mask2, 'km');
        title('\rmAM', 'fontsize', fntsz+1);
    elseif w == 2
        pcolorpsn(x.lat_area, x.lon_area, NMbV_PMp{n}{5}.*ice_mask2, 'km');
        title('\rmPM', 'fontsize', fntsz+1)
    end
    plotpsn(ice_lat, ice_lon, 'color', colo('k0'), 'linewidth', 1, 'km');
    colormap(cmapn);
    set(axs(w), 'clim', [0 150]);
    xlim([-0.7e3 0.9e3]);
    ylim([-3.5e3 -0.5e3 ]);
    
    if w == 2
        cb = colorbar;
        set(cb, 'fontsize', fntsz);
    end
    
    set(axs(w), 'position', [0.03+(w-1)*0.4, 0.1, 0.50, 0.75]);
    
    if w > 1
        set(axs(w), 'yticklabel', []);
    else
        ylabel('Northing [km]');
    end
    
    xlabel('Easting [km]');
    
end

text(axs(1), 12e2, -0.1e3, sprintf('Number of Melt Days: 36.5 GHz V-pol (%.0f)', yrs(1)), 'fontsize', 16, 'horizontalalignment', 'center');

% if savefigyes
%     saveas(gcf, ['figs220901_number_of_melt_days\', sprintf('Number_of_melt_days_KaV_%.0f_v%s.png', yrs(1), curdate ) ]); 
% end

%% Number of C V-pol melt days
cmapn = colormap(spring);
cmapn = flipud(cmapn);
cmapn(1,:) = [1,1,1];

figure; set(gcf, 'position', [138   304   744   638]);
for w = 1:2
    axs(w) = subplot(1,2,w);
    set(axs(w), 'fontsize', fntsz);
    greenland('color', colo('grey1'), 'km');
    if w == 1
        pcolorpsn(x.lat_area, x.lon_area, NMbV_AMp{n}{2}.*ice_mask2, 'km');
        title('\rmAM', 'fontsize', fntsz+1);
    elseif w == 2
        pcolorpsn(x.lat_area, x.lon_area, NMbV_PMp{n}{2}.*ice_mask2, 'km');
        title('\rmPM', 'fontsize', fntsz+1)
    end
    plotpsn(ice_lat, ice_lon, 'color', colo('k0'), 'linewidth', 1, 'km');
    colormap(cmapn);
    set(axs(w), 'clim', [0 150]);
    xlim([-0.7e3 0.9e3]);
    ylim([-3.5e3 -0.5e3 ]);
    
    if w == 2
        cb = colorbar;
        set(cb, 'fontsize', fntsz);
    end
    
    set(axs(w), 'position', [0.03+(w-1)*0.4, 0.1, 0.50, 0.75]);
    
    if w > 1
        set(axs(w), 'yticklabel', []);
    else
        ylabel('Northing [km]');
    end
    
    xlabel('Easting [km]');
    
end

text(axs(1), 12e2, -0.1e3, sprintf('Number of Melt Days: 6.9 GHz V-pol (%.0f)', yrs(1)), 'fontsize', 16, 'horizontalalignment', 'center');

% if savefigyes
%     saveas(gcf, ['figs220901_number_of_melt_days\', sprintf('Number_of_melt_days_CV_%.0f_v%s.png', yrs(1), curdate ) ]); 
% end


%% Number of L V-pol melt days
cmapn = colormap(spring);
cmapn = flipud(cmapn);
cmapn(1,:) = [1,1,1];

figure; set(gcf, 'position', [138   304   744   638]);
for w = 1:2
    axs(w) = subplot(1,2,w);
    set(axs(w), 'fontsize', fntsz);
    greenland('color', colo('grey1'), 'km');
    if w == 1
        pcolorpsn(x.lat_area, x.lon_area, NMbV_AMp{n}{1}.*ice_mask2, 'km');
        title('\rmAM', 'fontsize', fntsz+1);
    elseif w == 2
        pcolorpsn(x.lat_area, x.lon_area, NMbV_PMp{n}{1}.*ice_mask2, 'km');
        title('\rmPM', 'fontsize', fntsz+1)
    end
    plotpsn(ice_lat, ice_lon, 'color', colo('k0'), 'linewidth', 1, 'km');
    colormap(cmapn);
    set(axs(w), 'clim', [0 150]);
    xlim([-0.7e3 0.9e3]);
    ylim([-3.5e3 -0.5e3 ]);
    
    if w == 2
        cb = colorbar;
        set(cb, 'fontsize', fntsz);
    end
    
    set(axs(w), 'position', [0.03+(w-1)*0.4, 0.1, 0.50, 0.75]);
    
    if w > 1
        set(axs(w), 'yticklabel', []);
    else
        ylabel('Northing [km]');
    end
    
    xlabel('Easting [km]');
    
end

text(axs(1), 12e2, -0.1e3, sprintf('Number of Melt Days: 1.4 GHz V-pol (%.0f)', yrs(1)), 'fontsize', 16, 'horizontalalignment', 'center');

% if savefigyes
%     saveas(gcf, ['figs220901_number_of_melt_days\', sprintf('Number_of_melt_days_LV_%.0f_v%s.png', yrs(1), curdate ) ]); 
% end

%%
for seld = 1%:length(x.AMorPM)
    
    if x.AMorPM(seld) == 0
        ampmstr = 'AM';
    elseif x.AMorPM(seld) == 1
        ampmstr = 'PM';
    end
    
    figure; set(gcf, 'position', [19         134        1447        1190]);
    for w = 1:5
        axs(w) = subplot(2,5,w);
        greenland('color', colo('grey1'));
        pcolorpsn(x.lat_area, x.lon_area, TBVp{1}{w}(:,:,seld).*ice_mask2);
        plotpsn(ice_lat, ice_lon, 'color', colo('k0'), 'linewidth', 1);
        set(axs(w), 'clim', [100 300]);
        title(freqs{w}, 'fontsize', fntsz+1)
        xlim([-0.7e6 0.9e6]);
        ylim([-3.5e6 -0.5e6 ]);
        colormap('jet');
        if w == 5
            colorbar;
        end
        
        set(axs(w), 'position', [0.03+(w-1)*0.18, 0.55, 0.18, 0.4]);
        
        if w > 1
            set(axs(w), 'yticklabel', []);
        end
        
    end
        
    for w = 6:10
        axs(w) = subplot(2,5,w);
        greenland('color', colo('grey1'));
        pcolorpsn(x.lat_area, x.lon_area, TBHp{1}{w-5}(:,:,seld).*ice_mask2);
        plotpsn(ice_lat, ice_lon, 'color', colo('k0'), 'linewidth', 1);
        set(axs(w), 'clim', [100 300]);
        title(freqs{w-5}, 'fontsize', fntsz+1)
        xlim([-0.7e6 0.9e6]);
        ylim([-3.5e6 -0.5e6 ]);
        colormap('jet');
        if w == 10
            colorbar;
        end
        
        set(axs(w), 'position', [0.03+(w-5-1)*0.18, 0.05, 0.18, 0.40]);
        
        if w > 6
            set(axs(w), 'yticklabel', []);
        end
        
    end    
    
    text(axs(3), 0, -200000, sprintf('TB V-pol (%s %s)', datestr(dnums(seld)), ampmstr), 'fontsize', 16, 'horizontalalignment', 'center');
    text(axs(8), 0, -200000, sprintf('TB H-pol (%s %s)', datestr(dnums(seld)), ampmstr), 'fontsize', 16, 'horizontalalignment', 'center');
    
%             saveas(gcf, ['figs220819_daily_TB\', sprintf('TBV_%s_%s.png', datestr(dnums(seld), 'yyyymmdd'), ampmstr) ]);
%             close(gcf);
    
end

%% Melt profile for PM and AM on certain days
cmap_bp = flipud(bluepurple_colormap);

datesel = [...
    2019, 4, 25;
    2019, 5, 28;
    2019, 6,  9;
    2019, 8,  8;
    2019, 8, 10;
    2019, 8, 17;
    2019, 9, 25;
    2019, 10, 30];

dnumsel = datenum(datesel);

for ii = 1:length(datesel) %450%1:nz
    
    xix = find(dnumsel(ii) == dnums);
    
    if ~isempty(xix)
        
        seld = xix(end);
        
        
        if x.AMorPM(seld) == 0
            ampmstr1 = 'AM';
        elseif x.AMorPM(seld) == 1
            ampmstr1 = 'PM';
        end
        
        if x.AMorPM(seld+1) == 0
            ampmstr2 = 'AM';
        elseif x.AMorPM(seld+1) == 1
            ampmstr2 = 'PM';
        end
        
        figure; set(gcf, 'position', [138   304   744   638]);
        
        for w = 1:2
            axs(w) = subplot(1,2,w);
            set(axs(w), 'fontsize', fntsz);
            greenland('color', colo('grey1'));%, 'km');
            if w == 1
                pcolorpsn(x.lat_area, x.lon_area, MSP_Vp{1}(:,:,seld).*ice_mask2);%, 'km');
                title({'Snow Status Profile V-pol', sprintf('%s %s', datestr(dnums(seld)), ampmstr1)}, 'fontsize', fntsz)
            elseif w == 2
                pcolorpsn(x.lat_area, x.lon_area, MSP_Vp{1}(:,:,seld+1).*ice_mask2);%, 'km');
                title({'Snow Status Profile V-pol', sprintf('%s %s', datestr(dnums(seld+1)), ampmstr2)}, 'fontsize', fntsz)
            end
            plotpsn(ice_lat, ice_lon, 'color', colo('k0'), 'linewidth', 1);%, 'km');
            plotpsn([67,67], [-50,-43], 'color', colo('r0'), 'linewidth', 1, 'linewidth', 3);
            colormap(cmap_bp);
            set(axs(w), 'clim', [1 11]);
            xlim([-0.7e6 0.9e6]);
            ylim([-3.5e6 -0.5e6 ]);
            if w == 2
                cb = colorbar;
                set(cb, 'fontsize', fntsz);
                set(cb, 'ytick', linspace(1.5,10.5, length(lyrsV)))
                set(cb, 'yticklabel', 1:11);
            end
            set(axs(w), 'position', [0.03+(w-1)*0.4, 0.1, 0.50, 0.82]);
            if w > 1
                set(axs(w), 'yticklabel', []);
            else
                ylabel('Northing [m]');
            end
            xlabel('Easting [m]');
        end
        % text(axs(3), 0, -0.2, sprintf('Saturation Factor V-pol (%s %s)', datestr(dnums(seld)), ampmstr), 'fontsize', 16, 'horizontalalignment', 'center');
%             if savefigyes
%            saveas(gcf, ['figs230123_melt_status_map\', sprintf('Melt_status_profile_map_V_%s_%s_v%s.png', datestr(dnums(seld), 'yyyymmdd'), ampmstr, curdate ) ]); 
%            close(gcf);
%             end
    end
end


%% Melt status transect

cmap_d = [colo('grey1'); colo('blue1')];

for seld = 450%1:nz
    
    if x.AMorPM(seld) == 0
        ampmstr = 'AM';
    elseif x.AMorPM(seld) == 1
        ampmstr = 'PM';
    end
    
    XDX  = double(squeeze(MbV{n}(tran_inds,seld,:))');
    XDX2 = [XDX;XDX(5,:)];
    
    XSFV  = double(squeeze(SFV{n}(tran_inds,seld,:))');
      
    figure; hold all
    set(gcf, 'position', [489         502        1031         357]);
    
    ax1 = subplot(2,1,1); hold all;
    set(ax1, 'fontsize', 14);
    L1 = plot(lon_trans, flipud(XSFV), 'linewidth', 3);
    for jj = 1:length(L1)
        L2 = plot(ax1, [lon_trans(1), lon_trans(end)], [Zv(jj)*E_SD_SFV{n}(jj), Zv(jj)*E_SD_SFV{n}(jj)], '--', 'linewidth', 2, 'color', L1(6-jj).Color);
    end  
    title(['Melt Status Profile V-pol:', sprintf(' %s %s', datestr(dnums(seld)), ampmstr)], 'fontsize', fntsz)
    ylim([0 1]);
    xlim([lon_trans(1), lon_trans(end)+0.011]);
    lg = legend(fliplr(freqs), 'location', 'northeastoutside');
    set(ax1, 'xticklabel', []);
    ylabel(ax1, 'SF [-]');
    grid on;
   
    ax2 = subplot(2,1,2); hold all;
    set(ax2, 'fontsize', 14);
    s = pcolor(ones(6,1)*lon_trans', [0:1:5]'*ones(1,length(lon_trans)), XDX2);
    set(s, 'edgecolor', 'none');
    set(s, 'facealpha', 0.5);
    set(gca, 'clim', [0 1]);
    colormap(cmap_d);
    for jj = 1:4
        plot([lon_trans(1), lon_trans(end)], [jj,jj], 'color', colo('grey1'));
    end
    ylim([0 5]);
    xlim([lon_trans(1), lon_trans(end)+0.01]);
    set(ax2, 'ytick', [0.5:4.5]);
    set(ax2, 'yticklabel', [5:-1:1]); 
    xlabel(ax2, 'Longitude [^o]');
    ylabel(ax2, 'Layer');
    set(ax2, 'xgrid', 'on'); 
    
    set(ax1, 'position', [0.1, 0.55, 0.75, 0.36]);
    set(ax2, 'position', [0.1, 0.18, 0.75, 0.36]);
 
%     saveas(gcf, ['figs220901_melt_status_transect\', sprintf('Melt_status_transect_V_%s_%s_v%s.png', datestr(dnums(seld), 'yyyymmdd'), ampmstr, curdate ) ]); 
%     close(gcf);
    
end


