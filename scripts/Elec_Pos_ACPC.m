%% AC/PC comarison - b/w lead DBS and Raters
clear all
clc

%% Load files
main_dir = 'C:\Research\DBS_Retro\Final';
dir_dat = fullfile(main_dir,'data');
dir_out = fullfile(main_dir,'Output');

load(fullfile(dir_dat,'demogr_sum.mat'));
load(fullfile(dir_dat,'input_acpc_MNI_leaddbs','acpc_mni.mat'));
load(fullfile(dir_dat,'input_acpc_MNI_leaddbs','dbs_cont_mni.mat'));

% Load STL meshes for plots
brain_stl = stlread(fullfile(main_dir,'etc','Segmentation_brain.stl'));
L_stn = stlread(fullfile(main_dir,'etc','STN_Segment_L.stl'));
R_stn = stlread(fullfile(main_dir,'etc','STN_Segment_R.stl'));

%% Registration Error Calculations
% Displacement between transformed AC and PC coordinates and ground truth AC/PC coordinates (note this is "real world AFRE)

% Ground truth AC/PC coords, as defined by consensus among authors
ACPC_gt = [- 0.239,   1.880,   -4.750, -0.058, -24.683,  -2.360];

% Displacement between transformed AC/PC coords and ground truth coords
acpc_disp = acpc_mni(:,2:7,:) - repmat(ACPC_gt,size(acpc_mni,1),1,2);
acpc_disp(:,end+1,:) = vecnorm(acpc_disp(:,1:3,:),2,2);
acpc_disp(:,end+1,:) = vecnorm(acpc_disp(:,4:6,:),2,2);

% Mean AFRE across both raters for each subject (error in each
% axis and Euclidean Error)
AFRE_m = mean(acpc_disp,3);
AFRE_cat = [acpc_disp(:,:,1);acpc_disp(:,:,2)];

%% Calculate mean cont (electrode tip) position and deviation

% mean coordinate of distal contact across raters
cont_mni_m = mean(dbs_cont_mni(:,2:4,:,:),4);
% mean coordinate of left and right distal contact across all subjects
cont_mni_cen = mean(cont_mni_m,1);

% Contact deviation from mean position
cont_mc = cont_mni_m - repmat(cont_mni_cen,size(cont_mni_m,1),1,1);
cont_mc(:,end+1,:) = vecnorm(cont_mc,2,2);

% Determine outlier electrodes
cont_mc_cat = [cont_mc(:,:,1);cont_mc(:,:,2)];
% outliers >1.5 IQR of Euclidean error across all contacts
cont_out_cat = isoutlier(abs(cont_mc_cat(:,4)),"quartiles");
% Subjects to remove (if either left or right contact is an outlier)
cont_out = cont_out_cat(1:length(cont_mc)) + cont_out_cat(1+length(cont_mc):end)>0;

%% Recompute tip displacement and AFRE after removing outliers
cont_mni_mqc = mean(dbs_cont_mni(~cont_out,2:4,:,:),4);
cont_mni_cenqc = mean(cont_mni_mqc,1);
cont_mc_qc = cont_mni_mqc - repmat(cont_mni_cenqc,size(cont_mni_mqc,1),1,1);
cont_mc_qc(:,end+1,:) = vecnorm(cont_mc_qc,2,2);

AFRE_m_qc = AFRE_m(~cont_out,:);

cont_stats = compute_stats_sum(abs(cont_mc_qc(:,:,1)));
cont_stats = [cont_stats;compute_stats_sum(abs(cont_mc_qc(:,:,2)))];
cont_stats(:,end+1) = cont_stats(:,end).^2;
%% Figure 4a (top) - 3d plot of all electrode tips in MNI space
mkdir(fullfile(dir_out,'Figures','Fig4'));

cord = colororder;

figure
trimesh(R_stn,'FaceColor','k','EdgeColor','k','FaceAlpha',0.05,'EdgeAlpha',0.1)
hold on
scatter3(cont_mni_m(:,1,1),cont_mni_m(:,2,1),cont_mni_m(:,3,1),'filled','markeredgecolor','k','markerfacecolor',cord(1,:))
scatter3(cont_mni_cen(1,1,1),cont_mni_cen(1,2,1),cont_mni_cen(1,3,1),400,'filled','markeredgecolor','k','markerfacecolor',cord(2,:))

trimesh(L_stn,'FaceColor','k','EdgeColor','k','FaceAlpha',0.05,'EdgeAlpha',0.1)
hold on
scatter3(cont_mni_m(:,1,2),cont_mni_m(:,2,2),cont_mni_m(:,3,2),'filled','markeredgecolor','k','markerfacecolor',cord(1,:))
scatter3(cont_mni_cen(1,1,2),cont_mni_cen(1,2,2),cont_mni_cen(1,3,2),400,'filled','markeredgecolor','k','markerfacecolor',cord(2,:))
axis equal

set(gca,'fontsize',16)
axis equal

view(5,25)
set(gcf,'color','w','Position',  [100, 100, 1500, 800],'PaperSize', [16 9]);
% set(gcf,'color','w','Renderer','Painter','Position',  [100, 100, 1500, 800],'PaperSize', [16 9]);

saveas(gcf,fullfile(dir_out,'Figures','Fig4','Fig4a_mesh.pdf'));

%% Figure 3b (bottom) - boxchart of contact deviation

xlab = {'X-Axis','Y-Axis','Z-Axis','All Axes'};
si_or = [2,1];

figure
for s = 1:2
    si = si_or(s);
    subplot(1,2,s)
    boxchart(abs(cont_mc(:,:,si)))
    set(gca,'xticklabel',xlab)
    ylabel('Electrode Tip Displacement (mm)')

    ylim([0 10])
end
set(gcf,'color','w','Position',  [100, 100, 1000, 400],'PaperSize', [11 5]);
saveas(gcf,fullfile(dir_out,'Figures','Fig4','Fig4a_hist.pdf'));

%% Figure 4b - Correlating contact tip deviation with ACPC FRE

%fid2pl 0 = ac, 1 = pc);
fid2pl = 1;
fid2pl = fid2pl*3;

ax_lab = {'X-Axis','Y-Axis','Z-Axis'};
clear p_frecor
clear r2_frecor

pl_or = [2,4,6,1,3,5];

figure
count = 1;
for si = 1:2
    for ax = 1:3
        subplot(3,2,pl_or(count))
        count = count+1;

        mdl = fitlm(AFRE_m_qc(:,ax+fid2pl),cont_mc_qc(:,ax,si));
        p_frecor(si,ax) = mdl.Coefficients.pValue(2);
        r2_frecor(si,ax) = mdl.Rsquared.Ordinary;

        h = plot(mdl);
        dataHandle = findobj(h,'DisplayName','Data');
        fitHandle = findobj(h,'DisplayName','Fit');
        cbHandles = findobj(h,'DisplayName','Confidence bounds');
        cbHandles = findobj(h,'LineStyle',cbHandles.LineStyle, 'Color', cbHandles.Color);

        dataHandle.Color = [0 0.2 0.6];
        fitHandle.Color = [0 0 0];
        set(cbHandles, 'Color', [0,0,0])

        axis equal
        xlim([-5 5])
        ylim([-5 5])
        xticks(-4:2:4)
        yticks(-4:2:4)

        xline(0,'--')
        yline(0,'--')

        legend('hide')

        ylabel('Contact Displacement (mm)')
        xlabel('AFRE (mm)')
        title(['r2 = ', num2str(r2_frecor(si,ax)),' P-val = ',num2str(p_frecor(si,ax))])
    end
end

set(gcf,'color','w','Position',  [100, 100, 600, 800],'PaperSize', [7 9]);
saveas(gcf,fullfile(dir_out,'Figures','Fig4','Fig4b.pdf'));

%% Statistics

% pred_var - summary of indep variables, cols: Age, Disease
% duration, Sex (0=M), Implantation order (0 = L), rater (1 = MA/BS),
% modality (1 = mri), Side (0 = L), first logical (0 = yes)

cont_mcRL = [cont_mc_qc(:,:,1);cont_mc_qc(:,:,2)];
AFRE_dub = [AFRE_m_qc;AFRE_m_qc];

pred_var = demogr_sum(~cont_out,2:7);
pred2pl = [pred_var;pred_var];
pred2pl(1:sum(~cont_out),end+1) = 1;
pred2pl(:,end+1) = ~(pred2pl(:,4) == pred2pl(:,7));
pred_lab = {demogr_lab{2:end},'Side','Implanted First'};

pred2pl = pred2pl(:,[1:3,5:end]);
pred_lab = pred_lab([1:3,5:end]);


%% Histogram of electrode tip displacement (Fig S3)
mkdir(fullfile(dir_out,'Figures','FigS3'));

figure
for ax2pl = 1:3
    subplot(3,1,ax2pl)
    histogram(cont_mcRL(:,ax2pl),'binwidth',1,'orientation','horizontal')
end
set(gcf,'color','w','Position',  [100, 100, 200, 1200],'PaperSize', [3 13]);
saveas(gcf,fullfile(dir_out,'Figures','FigS3','FigS3_hist.pdf'));

%%

p_mvar = nan(9,3);
tstat_mvar = nan(9,3);
nonpar_stat = nan(9,3);
p_preds = nan(9,3);

for ax2pl = 1:3

    fre2pl = [ax2pl,ax2pl+3];
    
    figure;
    
    % Spearman rank correlation for continuous var (age and dz duration)
    for p = 1:2
        subplot(2,5,p)
        scatter(pred2pl(:,p),cont_mcRL(:,ax2pl))
        title(pred_lab{p})
        [nonpar_stat(p,ax2pl),p_preds(p,ax2pl)] = corr(pred2pl(:,p),cont_mcRL(:,ax2pl));
        yline(0)
        xlabel(num2str(p_preds(p,ax2pl)))
    end
    
        % Wilcoxon rank sum for binary var (rater pair, modality,side,implanted first)
    for p = 3:7
        subplot(2,5,p)
        boxchart(cont_mcRL(:,ax2pl), 'GroupByColor', pred2pl(:,p))
        title(pred_lab{p})
        [p_preds(p,ax2pl),~,stat] = ranksum(cont_mcRL(pred2pl(:,p)==1,ax2pl),cont_mcRL(pred2pl(:,p)==0,ax2pl));
        nonpar_stat(p,ax2pl) = stat.ranksum;
        yline(0)
        xlabel(num2str(p_preds(p,ax2pl)))
    end
    
    % Spearman rank correlation for AFRE at AC and PC
    subplot(2,5,8)
    scatter(AFRE_dub(:,fre2pl(1)),cont_mcRL(:,ax2pl))
    title('AFRE at AC')
    [nonpar_stat(8,ax2pl),p_preds(8,ax2pl)] = corr(AFRE_dub(:,fre2pl(1)),cont_mcRL(:,ax2pl));
    yline(0)
    xlabel(num2str(p_preds(8,ax2pl)))
    
    subplot(2,5,9)
    scatter(AFRE_dub(:,fre2pl(2)),cont_mcRL(:,ax2pl))
    title('AFRE at PC')
    [nonpar_stat(9,ax2pl),p_preds(9,ax2pl)] = corr(AFRE_dub(:,fre2pl(2)),cont_mcRL(:,ax2pl));
    yline(0)
    xlabel(num2str(p_preds(9,ax2pl)))
    
    % Multivariate linear regression including all variables
    subplot(2,5,10)
    mdl_fin = fitlm([pred2pl,AFRE_dub(:,fre2pl)],cont_mcRL(:,ax2pl));
    plot(mdl_fin)
    legend('hide')

    p_mvar(:,ax2pl) = mdl_fin.Coefficients.pValue(2:end);
    tstat_mvar(:,ax2pl) = mdl_fin.Coefficients.tStat(2:end);

    set(gcf,'color','w','Position',  [100, 100, 1400, 600],'PaperSize', [15 7]);
    saveas(gcf,fullfile(dir_out,'Figures','FigS3',['FigS3_',ax_lab{ax2pl},'.pdf']));

end
