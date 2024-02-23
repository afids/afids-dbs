%% AC/PC comparison - b/w lead DBS and Raters
clear all
clc

%% Load Files
main_dir = 'afids-dbs'; % path to cloned repo
dir_dat = fullfile(main_dir,'data');
dir_out = fullfile(main_dir,'Output');

load(fullfile(dir_dat,'input_acpc_native','acpc_native.mat'));
load(fullfile(dir_dat,'input_acpc_native','dbs_cont_native.mat'));
load(fullfile(dir_dat,'input_acpc_MNI_leaddbs','acpc_mni.mat'));
load(fullfile(dir_dat,'input_acpc_MNI_leaddbs','dbs_cont_mni.mat'));

%% Localization Error calculations

% Inter-rater displacement between rater placements (i.e. LE) for AC and PC
acpc_rat_diff = acpc_native(:,:,1) - acpc_native(:,:,2);

% Inter-rater displacement between rater placements (i.e. LE) for distal contacts
cont_rat_diff = dbs_cont_native(:,1:4,:,1) - dbs_cont_native(:,1:4,:,2);

% Generate matrix, dim 1 = pts, dim2 = R cont, L cont, AC, PC, dim 3 = X, Y, Z, Euc
Err_sum = zeros(size(cont_rat_diff,1),4,4);
for d = 1:3
    Err_sum(:,1,d) = abs(cont_rat_diff(:,d+1,1));
    Err_sum(:,2,d) = abs(cont_rat_diff(:,d+1,2));

    Err_sum(:,3,d) = abs(acpc_rat_diff(:,d+1));
    Err_sum(:,4,d) = abs(acpc_rat_diff(:,d+4));
end
Err_sum(:,:,4) = sqrt(sum(Err_sum(:,:,1:3).^2,3));

%% Summary Localization Error Metrics (Table S2)

tbl_cname = ["Metric","Min","LQ","Median","UQ","Max","Mean","Std"];
tbl_rname = ["RE_x","RE_y","RE_z","RE_Euc","LE_x","LE_y","LE_z","LE_Euc","AC_x","AC_y","AC_z","AC_Euc","PC_x","PC_y","PC_z","PC_Euc"];

FLE_calc = [Err_sum(:,:,1),Err_sum(:,:,2),Err_sum(:,:,3),Err_sum(:,:,4)];
FLE_tb = compute_stats_sum(FLE_calc);
FLE_tb = round(FLE_tb,2);
FLE_tb_fin = array2table([tbl_rname',FLE_tb],'VariableNames',tbl_cname);

writetable(FLE_tb_fin,fullfile(dir_out,'Tables','Table_S2.csv'));
%% Figure 1b - Localization Error
mkdir(fullfile(dir_out,'Figures','Fig1'));

% Box plot of XYZ and euc error of contact localization vs acpc localization
ax_lab = {'X-Axis','Y-Axis','Z-Axis','Euclidean Error'};
xlab = {'L Electrode','R Electrode','AC','PC'};

% p-val = signrank, row = L vs R, L vs AC, L vs PC, R vs AC, R vs PC, AC vs PC
p_sirank = zeros(6,4);

figure;
for ax = 1:4
    subplot(1,4,ax)
    boxchart(Err_sum(:,:,ax))

    ylim([0 5])
    title(ax_lab(ax))
    set(gca,'xticklabel',xlab)
    
    p_sirank(1,ax) = signrank(Err_sum(:,1,ax),Err_sum(:,2,ax));
    p_sirank(2,ax) = signrank(Err_sum(:,1,ax),Err_sum(:,3,ax));
    p_sirank(3,ax) = signrank(Err_sum(:,1,ax),Err_sum(:,4,ax));
    p_sirank(4,ax) = signrank(Err_sum(:,2,ax),Err_sum(:,3,ax));
    p_sirank(5,ax) = signrank(Err_sum(:,2,ax),Err_sum(:,4,ax));
    p_sirank(6,ax) = signrank(Err_sum(:,3,ax),Err_sum(:,4,ax));
end

subplot(1,4,1)
ylabel('Localization Error (mm)')

set(gcf,'color','w');
set(gcf, 'PaperSize', [15 7]);
set(gcf, 'Position',  [100, 100, 1400, 600])

saveas(gcf,fullfile(dir_out,'Figures','Fig1','Fig1b.pdf'));

%% Registration Error Calculations
% Displacement between transformed AC and PC coordinates and ground truth AC/PC coordinates (note this is "real world AFRE)

% Ground truth AC/PC coords, as defined by consensus among authors
ACPC_gt = [- 0.239,   1.880,   -4.750, -0.058, -24.683,  -2.360];

% Displacement between transformed AC/PC coords and ground truth coords
acpc_disp = acpc_mni(:,2:7,:) - repmat(ACPC_gt,size(acpc_mni,1),1,2);
acpc_disp(:,end+1,:) = vecnorm(acpc_disp(:,1:3,:),2,2);
acpc_disp(:,end+1,:) = vecnorm(acpc_disp(:,4:6,:),2,2);

% Mean AFRE across both raters for each subject (absolute error in each
% axis and Euclidean Error)
AFRE_m = mean(abs(acpc_disp),3);
AFRE_cat = [acpc_disp(:,:,1);acpc_disp(:,:,2)];

% RMSE error across axes (as reported by Schönecker et al. 2009)
RMSE_ax = [sqrt(mean(AFRE_m(:,1:3).^2,2)),sqrt(mean(AFRE_m(:,4:6).^2,2))];

%% Summary Registration Error Metrics (Table S3)

tbl_cname = ["Metric","Min","LQ","Median","UQ","Max","Mean","Std"];
tbl_rname = ["AC_x","AC_y","AC_z","AC_Euc","AC_RMSE","PC_x","PC_y","PC_z","PC_Euc","PC_RMSE","ACPC_x","ACPC_y","ACPC_z","ACPC_Euc","ACPC_RMSE"];

FRE_calc = [AFRE_m(:,1:3),AFRE_m(:,7),RMSE_ax(:,1),AFRE_m(:,4:6),AFRE_m(:,8),RMSE_ax(:,2)];
FRE_tb = round(compute_stats_sum(FRE_calc),2);

FRE2_calc = [FRE_calc(:,1:5);FRE_calc(:,6:10)];
FRE2_tb = round(compute_stats_sum(FRE2_calc),2);

FRE_tb_fin = array2table([tbl_rname',[FRE_tb;FRE2_tb]],'VariableNames',tbl_cname);

RMSE_acpc_dat = squeeze(mean(acpc_disp(:,1:6),3));
RMSE_comb_dat = [RMSE_acpc_dat(:,1:3);RMSE_acpc_dat(:,4:6)];
RMSE_acpc = sqrt(mean(RMSE_acpc_dat.^2))';
RMSE_comb = sqrt(mean(RMSE_comb_dat.^2))';

FRE_tb_fin.RMSE = round([RMSE_acpc(1:3);NaN;NaN;RMSE_acpc(4:6);NaN;NaN;RMSE_comb;NaN;NaN],2);

writetable(FRE_tb_fin,fullfile(dir_out,'Tables','Table_S3.csv'));
%% Figure 2a/b

mkdir(fullfile(dir_out,'Figures','Fig2'));

figure;
subplot(2,1,1)
scatter3(AFRE_cat(:,1),AFRE_cat(:,2),AFRE_cat(:,3),'filled','MarkerfaceColor',[0.8500, 0.3250, 0.0980])
hold on
scatter3(0,0,0,400,'filled','MarkerFaceColor',[0.9290    0.6940    0.1250],'MarkerEdgeColor',[0, 0, 0])
xlabel('X-Axis (mm)')
ylabel('Y-Axis (mm)')
zlabel('Z-Axis (mm)')
axis equal
view(310,20)

subplot(2,1,2)
scatter3(AFRE_cat(:,4),AFRE_cat(:,5),AFRE_cat(:,6),'filled','MarkerfaceColor',[0.8500, 0.3250, 0.0980])
hold on
scatter3(0,0,0,400,'filled','MarkerFaceColor',[0.9290    0.6940    0.1250],'MarkerEdgeColor',[0, 0, 0])
xlabel('X-Axis (mm)')
ylabel('Y-Axis (mm)')
zlabel('Z-Axis (mm)')
axis equal
view(310,20)

set(gcf, 'color','w','Renderer','Painter','PaperSize', [7 13],'Position',  [100, 100, 600, 1200])
saveas(gcf,fullfile(dir_out,'Figures','Fig2','Fig2ab.pdf'));

%% Statistics - Sign rank between AFRE and LE at AC and PC

abs_AFLE = abs(acpc_rat_diff(:,2:7));
abs_AFLE(:,end+1) = vecnorm(abs_AFLE(:,1:3),2,2);
abs_AFLE(:,end+1) = vecnorm(abs_AFLE(:,4:6),2,2);

% Signrank, paired comparison between AFRE and LE
p_sirank_bw = zeros(8,1);
for ax = 1:8
    p_sirank_bw(ax) = signrank(AFRE_m(:,ax),abs_AFLE(:,ax));
end


%% Figure 2c/d

ax_lab = {'X-Axis','Y-Axis','Z-Axis','All Axes','X-Axis','Y-Axis','Z-Axis','All Axes'};

acpc_lab = zeros(1,length(AFRE_m));
acpc_lab(end+1:end+length(AFRE_m)) = ones(1,length(AFRE_m));

figure;
count = 1;
for ax = [1:3,7,4:6,8]
    subplot(2,4,count)
    
    boxchart([abs_AFLE(:,ax);AFRE_m(:,ax)],'GroupbyColor',acpc_lab)
    title(ax_lab(count))
    set(gca,'xtick',[])
    ylim([0 5])
    xlabel(p_sirank_bw(ax))

    count = count+1;
end

set(gcf,'color','w','PaperSize', [11 11],'Position',  [100, 100, 1000, 1000]);
saveas(gcf,fullfile(dir_out,'Figures','Fig2','Fig2cd.pdf'));


%% Separate AFRE by raters
% 
% rat = strcmp(table2array(Rat_Mod_QC(:,2)),'MA');
% 
% acpc_dif_rat{1} = acpc_disp(rat,:,1);
% acpc_dif_rat{2} = acpc_disp(rat,:,2);
% acpc_dif_rat{3} = acpc_disp(~rat,:,1);
% acpc_dif_rat{4} = acpc_disp(~rat,:,2);
% 
% %% Figure S1 a/b - AFRE by rater
% 
% figure;
% 
% subplot(2,1,1)
% 
% for r = 1:4
%     rat2pl = acpc_dif_rat{r};
%     scatter3(rat2pl(:,1),rat2pl(:,2),rat2pl(:,3),'filled')
%     hold on
% end
% 
% scatter3(0,0,0,400,'filled','MarkerFaceColor',[0    0    0],'MarkerEdgeColor',[0, 0, 0])
% xlabel('X-Axis (mm)')
% ylabel('Y-Axis (mm)')
% zlabel('Z-Axis (mm)')
% axis equal
% view(310,20)
% 
% 
% subplot(2,1,2)
% 
% for r = 1:4
%     rat2pl = acpc_dif_rat{r};
%     scatter3(rat2pl(:,4),rat2pl(:,5),rat2pl(:,6),'filled')
%     hold on
% end
% 
% scatter3(0,0,0,400,'filled','MarkerFaceColor',[0    0    0],'MarkerEdgeColor',[0, 0, 0])
% xlabel('X-Axis (mm)')
% ylabel('Y-Axis (mm)')
% zlabel('Z-Axis (mm)')
% axis equal
% view(310,20)
% 
% 
% set(gcf,'color','w');
% set(gcf, 'PaperSize', [7 13]);
% set(gcf, 'Position',  [100, 100, 600, 1200])
% 
% 
% 
% %% Figure S1 c/d - Box AFRE by Rater
% 
% ax_lab = {'X-Axis','Y-Axis','Z-Axis','All Axes','X-Axis','Y-Axis','Z-Axis','All Axes'};
% 
% rat_cat = [rat;rat+2];
% 
% figure;
% 
% clear p_rater
% count = 1;
% for ax = [1:3,7,4:6,8]
%     subplot(2,4,count)
%     
% %     boxchart(abs(AFRE_cat(:,ax)),'GroupbyColor',rat_cat)
%     boxchart(AFRE_cat(:,ax),'GroupbyColor',rat_cat)
%     title(ax_lab(count))
%     set(gca,'xtick',[])
% %     ylim([0 5])
% %     xlabel(p_sirank_bw(ax))
% 
%     
%     p_rater(1,ax) = ranksum(AFRE_cat(rat_cat==0,1),AFRE_cat(rat_cat==1,1));
%     p_rater(2,ax) = ranksum(AFRE_cat(rat_cat==0,1),AFRE_cat(rat_cat==2,1));
%     p_rater(3,ax) = ranksum(AFRE_cat(rat_cat==0,1),AFRE_cat(rat_cat==3,1));
%     p_rater(4,ax) = ranksum(AFRE_cat(rat_cat==1,1),AFRE_cat(rat_cat==2,1));
%     p_rater(5,ax) = ranksum(AFRE_cat(rat_cat==1,1),AFRE_cat(rat_cat==3,1));
%     p_rater(6,ax) = ranksum(AFRE_cat(rat_cat==2,1),AFRE_cat(rat_cat==3,1));
% 
% 
%     count = count+1;
% end   
% set(gcf,'color','w');
% set(gcf, 'PaperSize', [11 11]);
% set(gcf, 'Position',  [100, 100, 1000, 1000])
