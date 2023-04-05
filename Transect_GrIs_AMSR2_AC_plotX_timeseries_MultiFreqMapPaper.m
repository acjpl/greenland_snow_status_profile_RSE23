
linclrs         =  linspecer(5);

[tansect_lat,tansect_long] ...
                = track2('gc',lat_avg,lon_min,lat_avg,lon_max);
            
barWidth        = 0.15;%bar width
alfa            = 0.2;%bar transparency

tair_lim_min    = -25;
tair_lim_max    =  25;

fntsz           = 24;
fntst           = 24;
%%
fig             = figure; hold all

set(fig, 'defaultAxesColorOrder', [left_color; right_color]);
set(fig, 'units','normalized','position', [0.0462    0.5840    0.6143    0.3174]);

hold all
set(gca,'fontsize',fntsz,'FontWeight','bold')
yyaxis left
h(1) = plot(lons, TRAN_L,  '--', 'color', linclrs(1,:), 'linewidth', 3); %-k
h(2) = plot(lons, TRAN_6,  '--', 'color', linclrs(2,:), 'linewidth', 3); %--r
h(3) = plot(lons, TRAN_10, '--', 'color', linclrs(3,:), 'linewidth', 3); %--g
h(4) = plot(lons, TRAN_18, '--', 'color', linclrs(4,:), 'linewidth', 3); %--k
h(5) = plot(lons, TRAN_36, '--', 'color', linclrs(5,:), 'linewidth', 3); %--c
h(6) = plot([ice_mask_lon_edge, ice_mask_lon_edge],Y_lim,'--','Color',[0.71, 0.59, 0.47],'linewidth', 1);

text(ice_mask_lon_edge, Y_lim(1)+0.01, 'Ice Edge',...
    'Color',[0.71, 0.59, 0.47],'fontsize',14,'FontWeight','bold','VerticalAlignment','bottom','horizontalalign','right')

text(floor(lon_min)+(ceil(lon_max)-floor(lon_min))*0.5, Y_lim(1)+(Y_lim(2)-Y_lim(1))*4, [VarName, ' ', dt_tmp{1}], ...
    'fontsize', 16, 'FontWeight', 'bold');

ylabel(VarName)
ylim(Y_lim)
xlim([floor(lon_min) ,ceil(lon_max)])
yyaxis right

if BAR_L < 0
    col_str_L   = 'b0';
    M1          = 0;
elseif BAR_L >= 0
    col_str_L   = 'r0';
    M1          = 1;
end
if BAR_M < 0
    col_str_M   = 'b0';
    M2          = 0;
elseif BAR_M >= 0
    col_str_M   = 'r0';
    M2          = 1;
end
if BAR_U < 0
    col_str_U   = 'b0';
    M3          = 0;
elseif BAR_U >= 0
    col_str_U   = 'r0';
    M3          = 1;
end
if BAR_B < 0
    col_str_B   = 'b0';
    M4          = 0;
elseif BAR_B >= 0
    col_str_B   = 'r0';
    M4          = 1;
end

left_color = [0 0 0];

if M1==1 || M2==1 || M3==1 || M4==1
    right_color = [1 0 0];
else
    right_color = [0 0 1];
end

% Temperature bars
h(7)  = bar(lon_KAN_L, BAR_L, 'FaceColor',colo(col_str_L),'EdgeColor',colo(col_str_L),'BarWidth', barWidth,'FaceAlpha',alfa);
text(lon_KAN_L, tair_lim_max, '        KAN_L',...
    'fontsize',fntst,'FontWeight','bold','VerticalAlignment','top','horizontalalign','center', 'interpreter', 'none')
hold on
h(8)  = bar(lon_KAN_M, BAR_M, 'FaceColor',colo(col_str_M),'EdgeColor',colo(col_str_M),'BarWidth', barWidth,'FaceAlpha',alfa);
text(lon_KAN_M, tair_lim_max, 'KAN_M',...
    'fontsize',fntst,'FontWeight','bold','VerticalAlignment','top','horizontalalign','center', 'interpreter', 'none')

h(9)  = bar(lon_KAN_U, BAR_U, 'FaceColor',colo(col_str_U),'EdgeColor',colo(col_str_U),'BarWidth', barWidth,'FaceAlpha',alfa);
text(lon_KAN_U, tair_lim_max, 'KAN_U',...
    'fontsize',fntst,'FontWeight','bold','VerticalAlignment','top','horizontalalign','center', 'interpreter', 'none')

h(10) = bar(lon_KAN_B, BAR_B, 'FaceColor',colo(col_str_B),'EdgeColor',colo(col_str_B),'BarWidth', barWidth,'FaceAlpha',alfa);
text(lon_KAN_B, tair_lim_max, 'KAN_B',...
    'fontsize',fntst,'FontWeight','bold','VerticalAlignment','top','horizontalalign','right', 'interpreter', 'none')

legend(h(1:5),'1.4 GHz (SMAP)',' 6.9 GHz (AMSR2)','10.7 GHz (AMSR2)','18.9 GHz (AMSR2)','36.5 GHz (AMSR2)',...
    'Location','northeastoutside')

grid on; box on;
grid minor
ylim([tair_lim_min tair_lim_max])
ylabel('Air Temp [^oC]');

xlabel('Longtitude (degree)')

title([dt_tmp{1}, ' PM']);

