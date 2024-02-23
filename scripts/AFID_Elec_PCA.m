%% AC/PC comarison - b/w lead DBS and Raters

clear all
clc


%% Load files


dir_o = 'C:\Research\DBS_Retro\Output\';


load(strcat(dir_o,'Patient_info.mat'));
load(strcat(dir_o,'Cont_coor.mat'));
load(strcat(dir_o,'acpc_MNI_QC.mat'));


%% ACPC AFRE
% Inter-rater errors:
% Diff within ldbs and raters (i.e. AFLE)
acpc_mni_diff = acpc_MNI_QC(:,:,1) - acpc_MNI_QC(:,:,2);

% Use contact coor in native space, last contact coor
cont_rat_diff = cont_rat_QC(:,1:4,:,1) - cont_rat_QC(:,1:4,:,2);

% Diff b/w mean ldbs and mean rater placement (AFRE, note this is "real
% world AFRE)

% ACPC_gold = [- 0.25,   1.298,   -5.003, -0.188, -24.756,  -2.376];
% ACPC_gold = mean(acpc_MNI_QC(:,2:7,:),[1,3]);
ACPC_gold = [- 0.07,   1.833,   -4.65, -0.058, -24.683,  -2.360];

acpc_bw_diff = acpc_MNI_QC(:,2:7,:) - repmat(ACPC_gold,size(acpc_MNI_QC,1),1,2);
acpc_bw_diff(:,end+1,:) = vecnorm(acpc_bw_diff(:,1:3,:),2,2);
acpc_bw_diff(:,end+1,:) = vecnorm(acpc_bw_diff(:,4:6,:),2,2);

AFRE_m = squeeze(mean(acpc_bw_diff,3));
AFRE_cat = [acpc_bw_diff(:,:,1);acpc_bw_diff(:,:,2)];

%% Calculate mean cont (electrode tip) position and deviation

cont_mni_m = mean(cont_mni_QC(:,2:4,:,:),4);

% Use contact coor in mni space, last contact coor
cont_cat = [cont_mni_QC(:,2:4,:,1);cont_mni_QC(:,2:4,:,2)];
cont_mni_cen = mean(cont_cat,1);

cont_cen = cont_cat - repmat(cont_mni_cen,size(cont_cat,1),1,1);
cont_cen(:,end+1,:) = vecnorm(cont_cen,2,2);

cont_mc = cont_mni_m - repmat(cont_mni_cen,size(cont_mni_m,1),1,1);
cont_mc(:,end+1,:) = vecnorm(cont_mc,2,2);

% cont_mc_med = squeeze(median(abs(cont_mc),1))';
% cont_mc_iqr = quantile(abs(cont_mc),[0.25 0.75],1);

cont_mc_med = squeeze(median(cont_mc,1))';
cont_mc_iqr = quantile(cont_mc,[0.25 0.75],1);

% Determine outlier electrodes
cont_mc_cat = [cont_mc(:,:,1);cont_mc(:,:,2)];

cont_out_cat = isoutlier(abs(cont_mc_cat(:,4)),"quartiles");
% cont_out_cat = sum(isoutlier(abs(cont_mc_cat(:,1:3)),"quartiles"),2)>0;

cont_out = cont_out_cat(1:length(cont_mc)) + cont_out_cat(1+length(cont_mc):end)>0;


% cont_out = or(cont_out,isoutlier(AFRE_m(:,8),"quartiles"));
% cont_out = or(cont_out,sum(isoutlier(AFRE_m(:,7:8),'quartiles'),2)>0);
% cont_out = or(cont_out,sum(isoutlier(AFRE_m(:,1:6),'quartiles'),2)>0);

% cont_out = sum(isoutlier(abs(cont_mc(:,4,:)),"quartiles"),3)>0;


cont_mc_qc = cont_mc(~cont_out,:,:);

AFRE_m_qc = AFRE_m(~cont_out,:);


%% load STLs

brain_stl = stlread('C:\Research\DBS_Retro\Input\maps\Brain\Segmentation_brain.stl');

L_stn = stlread('C:\Research\DBS_Retro\Input\maps\Left\STN_Segment_1.stl');
R_stn = stlread('C:\Research\DBS_Retro\Input\maps\Right\STN_Segment_1.stl');


R_stn_cen = triangulation(R_stn.ConnectivityList,R_stn.Points-squeeze(cont_mni_cen(:,:,1)));
L_stn_cen = triangulation(L_stn.ConnectivityList,L_stn.Points-squeeze(cont_mni_cen(:,:,2)));

%% Figure 3a (top) - 3d plot of all electrode tips in MNI space
cord = colororder;

figure
trimesh(R_stn,'FaceColor','k','EdgeColor','k','FaceAlpha',0.05,'EdgeAlpha',0.1)
hold on
scatter3(cont_mni_m(:,1,1),cont_mni_m(:,2,1),cont_mni_m(:,3,1),'filled','markeredgecolor','k','markerfacecolor','k','markerfacealpha',0.4)
scatter3(cont_mni_cen(1,1,1),cont_mni_cen(1,2,1),cont_mni_cen(1,3,1),400,'filled','markeredgecolor','k','markerfacecolor','k')


set(gca,'fontsize',16)
axis equal

% subplot(1,2,1)
trimesh(L_stn,'FaceColor','k','EdgeColor','k','FaceAlpha',0.05,'EdgeAlpha',0.1)
hold on
scatter3(cont_mni_m(:,1,2),cont_mni_m(:,2,2),cont_mni_m(:,3,2),'filled','markeredgecolor','k','markerfacecolor','k','markerfacealpha',0.4)
scatter3(cont_mni_cen(1,1,2),cont_mni_cen(1,2,2),cont_mni_cen(1,3,2),400,'filled','markeredgecolor','k','markerfacecolor','k')
axis equal

acpc_MNIm = squeeze(mean(acpc_MNI_QC(:,2:7),3));

scatter3(acpc_MNIm(:,1),acpc_MNIm(:,2),acpc_MNIm(:,3),'filled','markeredgecolor','k','markerfacecolor',cord(2,:),'markerfacealpha',0.6)
scatter3(ACPC_gold(1),ACPC_gold(2),ACPC_gold(3),400,'filled','markeredgecolor','k','markerfacecolor',cord(2,:))

scatter3(acpc_MNIm(:,4),acpc_MNIm(:,5),acpc_MNIm(:,6),'filled','markeredgecolor','k','markerfacecolor',cord(1,:),'markerfacealpha',0.6)
scatter3(ACPC_gold(4),ACPC_gold(5),ACPC_gold(6),400,'filled','markeredgecolor','k','markerfacecolor',cord(1,:))



view(5,25)
set(gcf,'color','w','Position',  [100, 100, 1500, 800]);
set(gcf, 'PaperSize', [16 9]);





%%

cont_RL_dif = [cont_mc_qc(:,1,1)-cont_mc_qc(:,1,2),cont_mc_qc(:,2,1)-cont_mc_qc(:,2,2),cont_mc_qc(:,3,1)-cont_mc_qc(:,3,2)];
[coeff_contd,~,~,~,exp_contd] = pca(cont_RL_dif);
pc1_contdline = [coeff_contd(:,1)*-4,coeff_contd(:,1),coeff_contd(:,1)*4];
pc2_contdline = [coeff_contd(:,2)*-2,coeff_contd(:,2),coeff_contd(:,2)*2];
pc3_contdline = [coeff_contd(:,3)*-2,coeff_contd(:,3),coeff_contd(:,3)*2];

fid_dev = AFRE_m_qc(:,1:6) - repmat(mean(AFRE_m_qc(:,1:6)),length(AFRE_m_qc),1);

[coeff_ac,~,~,~,exp_ac] = pca(fid_dev(:,1:3));
[coeff_pc,~,~,~,exp_pc] = pca(fid_dev(:,4:6));


cont_mcat = [cont_mc_qc(:,1:3,1),cont_mc_qc(:,1:3,2)];
[coeff_cont,~,~,~,exp_cont] = pca(cont_mcat(:,1:6));
[coeff_contR,~,~,~,exp_contR] = pca(cont_mc_qc(:,1:3,1));
[coeff_contL,~,~,~,exp_contL] = pca(cont_mc_qc(:,1:3,2));

pc1_contlr = [coeff_contR(:,1)*-4,coeff_contR(:,1),coeff_contR(:,1)*4];
pc2_contlr = [coeff_contR(:,2)*-3,coeff_contR(:,2),coeff_contR(:,2)*3];
pc3_contlr = [coeff_contR(:,3)*-3,coeff_contR(:,3),coeff_contR(:,3)*3];

pc1_contll = [coeff_contL(:,1)*-4,coeff_contL(:,1),coeff_contL(:,1)*4];
pc2_contll = [coeff_contL(:,2)*-3,coeff_contL(:,2),coeff_contL(:,2)*3];
pc3_contll = [coeff_contL(:,3)*-3,coeff_contL(:,3),coeff_contL(:,3)*3];

pc1_acline = [coeff_ac(:,1)*-4,coeff_ac(:,1),coeff_ac(:,1)*4];
pc2_acline = [coeff_ac(:,2)*-3,coeff_ac(:,2),coeff_ac(:,2)*3];

pc1_pcline = [coeff_pc(:,1)*-4,coeff_pc(:,1),coeff_pc(:,1)*4];
pc2_pcline = [coeff_pc(:,2)*-3,coeff_pc(:,2),coeff_pc(:,2)*3];


%%

ac_rot = fid_dev(:,1:3)*coeff_ac;
pc_rot = fid_dev(:,4:6)*coeff_pc;

ac_rot_cor = corr(ac_rot).^2;
pc_rot_cor = corr(pc_rot).^2;
fid_rot_cor = corr([ac_rot,pc_rot]).^2;

ac_dev_cor = corr(fid_dev(:,1:3)).^2;
pc_dev_cor = corr(fid_dev(:,1:3)).^2;
fid_dev_cor = corr(fid_dev).^2;

[coeff_fids,~,~,~,exp_fids] = pca(fid_dev);

figure
% imagesc(fids_multicor)
% subplot(1,2,1)
% heatmap(ac_dev_cor,'ColorLimits',[0 1])
% subplot(1,2,2)
% heatmap(pc_dev_cor,'ColorLimits',[0 1])
subplot(1,3,1)
heatmap(fid_dev_cor,'ColorLimits',[0 1])
subplot(1,3,2)
heatmap(fid_rot_cor,'ColorLimits',[0 1])
subplot(1,3,3)
heatmap(coeff_fids.^2,'ColorLimits',[0 1])

set(gcf,'color','w','Position',  [100, 100, 1800, 400]);





%%

fids_rot = fid_dev*coeff_fids;

cont_R_rot = cont_mc_qc(:,1:3,1)*coeff_contR;
cont_L_rot = cont_mc_qc(:,1:3,2)*coeff_contL;

% cont_cor = corr(fid_dev(:,1:6),[cont_mc_qc(:,1:3,2),cont_mc_qc(:,1:3,1)]).^2
% cont_cor = corr(fids_rot,[cont_mc_qc(:,1:3,2),cont_mc_qc(:,1:3,1)]).^2
% cont_cor = corr(fids_rot,[cont_L_rot,cont_R_rot]).^2
cont_cor = corr([ac_rot,pc_rot],[cont_L_rot,cont_R_rot]).^2

figure
heatmap(round(cont_cor,2),'ColorLimits',[0 0.25])
set(gcf,'color','w','Position',  [100, 100, 700, 600]);


%%

% heatmap(corr(cont_mc_qc(:,1:3,1),fid_dev(:,4:6)).^2)
% heatmap(corr(cont_L_rot,pc_rot).^2)
% heatmap(corr(cont_mc_qc(:,1:3,1),cont_mc_qc(:,1:3,2)).^2)
% heatmap(corr(cont_R_rot,cont_L_rot).^2)

%%


figure

scatter3(fid_dev(:,1),fid_dev(:,2),fid_dev(:,3),'MarkerEdgeColor',cord(2,:))
hold on
scatter3(fid_dev(:,4),fid_dev(:,5),fid_dev(:,6),'MarkerEdgeColor',cord(1,:))

plot3(pc1_acline(1,:),pc1_acline(2,:),pc1_acline(3,:),'color',cord(2,:),'linewidth',2)
plot3(pc1_pcline(1,:),pc1_pcline(2,:),pc1_pcline(3,:),'color',cord(1,:),'linewidth',2)

plot3(pc2_acline(1,:),pc2_acline(2,:),pc2_acline(3,:),'color',cord(2,:),'linewidth',2)
plot3(pc2_pcline(1,:),pc2_pcline(2,:),pc2_pcline(3,:),'color',cord(1,:),'linewidth',2)

axis equal





%%

ac_PC1cor = sum(fid_dev(:,1:3).*repmat(coeff_ac(:,1)',length(fid_dev),1),2);
ac_PC2cor = sum(fid_dev(:,1:3).*repmat(coeff_ac(:,2)',length(fid_dev),1),2);
ac_PC3cor = sum(fid_dev(:,1:3).*repmat(coeff_ac(:,3)',length(fid_dev),1),2);

pc_PC1cor = sum(fid_dev(:,4:6).*repmat(coeff_pc(:,1)',length(fid_dev),1),2);
pc_PC2cor = sum(fid_dev(:,4:6).*repmat(coeff_pc(:,2)',length(fid_dev),1),2);
pc_PC3cor = sum(fid_dev(:,4:6).*repmat(coeff_pc(:,3)',length(fid_dev),1),2);


% figure
% scatter(ac_PC3cor,pc_PC3cor)

resid_ac = [];
mdl = fitlm([pc_PC1cor,pc_PC2cor,pc_PC3cor],ac_PC1cor);
resid_ac = [resid_ac,mdl.Residuals.Raw];
mdl = fitlm([pc_PC1cor,pc_PC2cor,pc_PC3cor],ac_PC2cor);
resid_ac = [resid_ac,mdl.Residuals.Raw];
mdl = fitlm([pc_PC1cor,pc_PC2cor,pc_PC3cor],ac_PC3cor);
resid_ac = [resid_ac,mdl.Residuals.Raw];

figure
scatter3(resid_ac(:,1),resid_ac(:,2),resid_ac(:,3),'MarkerEdgeColor','k')
hold on
scatter3(ac_PC1cor,ac_PC2cor,ac_PC3cor,'MarkerEdgeColor','b')
axis equal

%%


% figure
% scatter(ac_PC3cor,pc_PC3cor)

Rc_PC1cor = sum(cont_mc_qc(:,1:3,1).*repmat(coeff_contR(:,1)',length(fid_dev),1),2);
Rc_PC2cor = sum(cont_mc_qc(:,1:3,1).*repmat(coeff_contR(:,2)',length(fid_dev),1),2);
Rc_PC3cor = sum(cont_mc_qc(:,1:3,1).*repmat(coeff_contR(:,3)',length(fid_dev),1),2);

resid_contR = [];
mdl = fitlm([pc_PC1cor,pc_PC2cor,pc_PC3cor],Rc_PC1cor);
resid_contR = [resid_ac,mdl.Residuals.Raw];
mdl = fitlm([pc_PC1cor,pc_PC2cor,pc_PC3cor],Rc_PC2cor);
resid_contR = [resid_ac,mdl.Residuals.Raw];
mdl = fitlm([pc_PC1cor,pc_PC2cor,pc_PC3cor],Rc_PC3cor);
resid_contR = [resid_ac,mdl.Residuals.Raw];

figure
scatter3(resid_contR(:,1),resid_contR(:,2),resid_contR(:,3),'MarkerEdgeColor','k')
hold on
scatter3(Rc_PC1cor,Rc_PC2cor,Rc_PC3cor,'MarkerEdgeColor','b')
axis equal

%%
% figure
% scatter(ac_PC3cor,pc_PC3cor)

Lc_PC1cor = sum(cont_mc_qc(:,1:3,2).*repmat(coeff_contL(:,1)',length(fid_dev),1),2);
Lc_PC2cor = sum(cont_mc_qc(:,1:3,2).*repmat(coeff_contL(:,2)',length(fid_dev),1),2);
Lc_PC3cor = sum(cont_mc_qc(:,1:3,2).*repmat(coeff_contL(:,3)',length(fid_dev),1),2);

resid_contL = [];
mdl = fitlm([pc_PC1cor,pc_PC2cor,pc_PC3cor],Lc_PC1cor);
resid_contL = [resid_ac,mdl.Residuals.Raw];
mdl = fitlm([pc_PC1cor,pc_PC2cor,pc_PC3cor],Lc_PC2cor);
resid_contL = [resid_ac,mdl.Residuals.Raw];
mdl = fitlm([pc_PC1cor,pc_PC2cor,pc_PC3cor],Lc_PC3cor);
resid_contL = [resid_ac,mdl.Residuals.Raw];

figure
scatter3(resid_contL(:,1),resid_contL(:,2),resid_contL(:,3),'MarkerEdgeColor','k')
hold on
scatter3(Lc_PC1cor,Lc_PC2cor,Lc_PC3cor,'MarkerEdgeColor','b')
axis equal



%% Generate unit vecs
z_vals = -1:0.05:1;
sc_res = sqrt(1 - abs(z_vals).^2);

theta = 0:2*pi/(length(z_vals)-1):2*pi;

unit_vecs = [];
for z = 1:length(z_vals)

    sc2use = sc_res(z);
    m = sc2use'.*[cos(theta') sin(theta')];
    temp_vecs = [m,repmat(z_vals(z),length(theta),1)];

    unit_vecs = [unit_vecs;temp_vecs];

end


%%

num_units = length(unit_vecs);

clear r2_vecs
clear p_vecs
for u = 1:num_units
    ac_u = sum(fid_dev(:,1:3).*repmat(unit_vecs(u,:),length(fid_dev),1),2);
    pc_u = sum(fid_dev(:,4:6).*repmat(unit_vecs(u,:),length(fid_dev),1),2);

    mdl = fitlm(ac_u,pc_u);

    r2_vecs(u) = mdl.Rsquared.Ordinary;
    p_vecs(u) = mdl.Coefficients.pValue(2);

end
%%


cmap = parula(256);
v = rescale(r2_vecs, 1, 256,'InputMin',0,'InputMax',0.5); 
numValues = length(r2_vecs);
markerColors = zeros(numValues, 3);

% Now assign marker colors according to the value of the data.
for k = 1 : numValues
    row = round(v(k));
    markerColors(k, :) = cmap(row, :);
end

figure
% subplot(1,2,2)
scatter3(unit_vecs(:,1)*1,unit_vecs(:,2)*1,unit_vecs(:,3)*1, 200, markerColors,'filled','markerfacealpha',0.5);
hold on
scatter3(0,0,0,400,'filled','markerfacecolor','k')
axis equal

xlabel('X Axis')

[~,idx_maxr2] = max(r2_vecs);
[~,idx_minr2] = min(r2_vecs);

max_vcline = [unit_vecs(idx_maxr2,:)'*-5,unit_vecs(idx_maxr2,:)',unit_vecs(idx_maxr2,:)'*5];
min_vcline = [unit_vecs(idx_minr2,:)'*-5,unit_vecs(idx_minr2,:)',unit_vecs(idx_minr2,:)'*5];
plot3(max_vcline(1,:),max_vcline(2,:),max_vcline(3,:),'g','linewidth',2)
plot3(min_vcline(1,:),min_vcline(2,:),min_vcline(3,:),'r','linewidth',2)


%%
max_vec = unit_vecs(idx_maxr2,:);
ac_max = sum(fid_dev(:,1:3).*repmat(max_vec,length(fid_dev),1),2);
pc_max = sum(fid_dev(:,4:6).*repmat(max_vec,length(fid_dev),1),2);

figure
mdl = fitlm(ac_max,pc_max);
plot(mdl)
pc_res = mdl.Residuals.Raw;


%%


% fid_dev_dif = [zscore(fid_dev(:,1))-zscore(fid_dev(:,4)),zscore(fid_dev(:,2))-zscore(fid_dev(:,5)),zscore(fid_dev(:,3))-zscore(fid_dev(:,6))];
fid_dev_dif = [zscore(fid_dev(:,1)-fid_dev(:,4)),zscore(fid_dev(:,2)-fid_dev(:,5)),zscore(fid_dev(:,3)-fid_dev(:,6))];
[coeff_fidd,~,~,~,exp_fidd] = pca(fid_dev_dif);
pc1_fidline = [coeff_fidd(:,1)*-4,coeff_fidd(:,1),coeff_fidd(:,1)*4];
pc2_fidline = [coeff_fidd(:,2)*-2,coeff_fidd(:,2),coeff_fidd(:,2)*2];
pc3_fidline = [coeff_fidd(:,3)*-2,coeff_fidd(:,3),coeff_fidd(:,3)*2];

figure
scatter3(fid_dev_dif(:,1),fid_dev_dif(:,2),fid_dev_dif(:,3),'MarkerEdgeColor','k')
hold on
plot3(pc1_fidline(1,:),pc1_fidline(2,:),pc1_fidline(3,:),'k','linewidth',2)
plot3(pc2_fidline(1,:),pc2_fidline(2,:),pc2_fidline(3,:),'k','linewidth',2)

axis equal

%%
















%%


figure
subplot(1,2,1)
scatter3(cont_mc_qc(:,1,1),cont_mc_qc(:,2,1),cont_mc_qc(:,3,1),'MarkerEdgeColor','k')
hold on
scatter3(cont_mc_qc(:,1,2),cont_mc_qc(:,2,2),cont_mc_qc(:,3,2),'MarkerEdgeColor','b')
plot3(pc1_contlr(1,:),pc1_contlr(2,:),pc1_contlr(3,:),'k','linewidth',2)
plot3(pc1_contll(1,:),pc1_contll(2,:),pc1_contll(3,:),'b','linewidth',2)
plot3(pc2_contlr(1,:),pc2_contlr(2,:),pc2_contlr(3,:),'k','linewidth',2)
plot3(pc2_contll(1,:),pc2_contll(2,:),pc2_contll(3,:),'b','linewidth',2)
% 
% plot3(pc1_pcline(1,:),pc1_pcline(2,:),pc1_pcline(3,:),'g','linewidth',2)
% plot3(pc2_pcline(1,:),pc2_pcline(2,:),pc2_pcline(3,:),'g','linewidth',2)

axis equal


subplot(1,2,2)
scatter3(cont_RL_dif(:,1),cont_RL_dif(:,2),cont_RL_dif(:,3),'MarkerEdgeColor','k')
hold on
plot3(pc1_contdline(1,:),pc1_contdline(2,:),pc1_contdline(3,:),'k','linewidth',2)
plot3(pc2_contdline(1,:),pc2_contdline(2,:),pc2_contdline(3,:),'k','linewidth',2)
axis equal
xlabel('X-Axis')







%%












%%


figure

% subplot(2,2,1)

% scatter3(cont_mcRL(:,1),cont_mcRL(:,2),cont_mcRL(:,3),'MarkerEdgeColor','k')
scatter3(cont_mc_qc(:,1,si),cont_mc_qc(:,2,si),cont_mc_qc(:,3,si),'MarkerEdgeColor','k')
hold on
plot3(pc1_contline(1,:),pc1_contline(2,:),pc1_contline(3,:),'k','linewidth',2)
plot3(pc2_contline(1,:),pc2_contline(2,:),pc2_contline(3,:),'k','linewidth',2)
plot3(pc3_contline(1,:),pc3_contline(2,:),pc3_contline(3,:),'k','linewidth',2)
scatter3(0,0,0,400,'filled','MarkerFaceColor','k')
axis equal

trimesh(R_stn_cen,'FaceColor','k','EdgeColor','k','FaceAlpha',0.05,'EdgeAlpha',0.1)

% 
% 
% 
% subplot(2,2,2)
% 
% scatter3(fid_dev(:,1),fid_dev(:,2),fid_dev(:,3),'MarkerEdgeColor','r')
% hold on
% plot3(pc1_acline(1,:),pc1_acline(2,:),pc1_acline(3,:),'r','linewidth',2)
% plot3(pc2_acline(1,:),pc2_acline(2,:),pc2_acline(3,:),'r','linewidth',2)
% scatter3(0,0,0,400,'filled','MarkerFaceColor','k')
% axis equal
% 
% 
% subplot(2,2,3)
% scatter3(fid_dev(:,4),fid_dev(:,5),fid_dev(:,6),'MarkerEdgeColor','b')
% hold on
% plot3(pc1_pcline(1,:),pc1_pcline(2,:),pc1_pcline(3,:),'b','linewidth',2)
% plot3(pc2_pcline(1,:),pc2_pcline(2,:),pc2_pcline(3,:),'b','linewidth',2)
% scatter3(0,0,0,400,'filled','MarkerFaceColor','k')
% axis equal
% 
% 
% subplot(2,2,4)
% scatter3(0,0,0,400,'filled','MarkerFaceColor','k')
% hold on
% plot3(pc1_contline(1,:),pc1_contline(2,:),pc1_contline(3,:),'k','linewidth',2)
% plot3(pc2_contline(1,:),pc2_contline(2,:),pc2_contline(3,:),'k','linewidth',2)
% plot3(pc1_acline(1,:),pc1_acline(2,:),pc1_acline(3,:),'r','linewidth',2)
% plot3(pc2_acline(1,:),pc2_acline(2,:),pc2_acline(3,:),'r','linewidth',2)
% plot3(pc1_pcline(1,:),pc1_pcline(2,:),pc1_pcline(3,:),'b','linewidth',2)
% plot3(pc2_pcline(1,:),pc2_pcline(2,:),pc2_pcline(3,:),'b','linewidth',2)
% scatter3(0,0,0,400,'filled','MarkerFaceColor','k')

axis equal

set(gcf,'color','w','Position',  [100, 100, 700, 800],'PaperSize', [3 13]);


%%


pc2pl = 3;

% cont_devpc1 = sum(cont_mcRL(:,1:3).*repmat(coeff_cont(:,pc2pl)',length(cont_mcRL),1),2);
cont_devpc1 = sum(cont_mc_qc(:,1:3,si).*repmat(coeff_cont(:,pc2pl)',length(cont_mc_qc),1),2);

% ac_devpc1 = sum(fid_dev(:,1:3).*repmat(coeff_ac(:,3)',length(cont_mc_qc),1),2);
% pc_devpc1 = sum(fid_dev(:,4:6).*repmat(coeff_pc(:,3)',length(cont_mc_qc),1),2);

ac_devpc1 = sum(fid_dev(:,1:3).*repmat(coeff_cont(:,pc2pl)',length(fid_dev),1),2);
pc_devpc1 = sum(fid_dev(:,4:6).*repmat(coeff_cont(:,pc2pl)',length(fid_dev),1),2);

figure
% mdl = fitlm([pc_devpc1;pc_devpc1],cont_devpc1)
mdl = fitlm(pc_devpc1,cont_devpc1)
% scatter([pc_devpc1;pc_devpc1],cont_devpc1)
plot(mdl)
axis equal


%%

z_vals = -1:0.05:1;
sc_res = sqrt(1 - abs(z_vals).^2);

theta = 0:2*pi/(length(z_vals)-1):2*pi;

unit_vecs = [];
for z = 1:length(z_vals)

    sc2use = sc_res(z);
    m = sc2use'.*[cos(theta') sin(theta')];
    temp_vecs = [m,repmat(z_vals(z),length(theta),1)];

    unit_vecs = [unit_vecs;temp_vecs];

end
% scatter3(unit_vecs(:,1),unit_vecs(:,2),unit_vecs(:,3))
% axis equal

%%

num_units = length(unit_vecs);

clear r2_vecs
clear p_vecs
for u = 1:num_units
%     cont_u = sum(cont_mcRL(:,1:3).*repmat(unit_vecs(u,:),length(cont_mcRL),1),2);
    cont_u = sum(cont_mc_qc(:,1:3,si).*repmat(unit_vecs(u,:),length(cont_mc_qc),1),2);

%     cont_u = sum(cont_mc_qc(:,1:3,si).*repmat(coeff_cont(:,1)',length(cont_mc_qc),1),2);

%     pc_u = sum(fid_dev(:,4:6).*repmat(unit_vecs(u,:),length(fid_dev),1),2);
%     pc_u = sum(fid_dev(:,4:6).*repmat(coeff_pc(:,1)',length(fid_dev),1),2);
    pc_u = fid_dev(:,4:6)*coeff_pc(:,1);

%     mdl = fitlm([pc_u;pc_u],cont_u);
    mdl = fitlm(pc_u,cont_u);

    r2_vecs(u) = mdl.Rsquared.Ordinary;
    p_vecs(u) = mdl.Coefficients.pValue(2);

end

%%

cmap = parula(256);
v = rescale(r2_vecs, 1, 256,'InputMin',0,'InputMax',0.25); 
numValues = length(r2_vecs);
markerColors = zeros(numValues, 3);

% Now assign marker colors according to the value of the data.
for k = 1 : numValues
    row = round(v(k));
    markerColors(k, :) = cmap(row, :);
end

figure
% subplot(1,2,2)
scatter3(unit_vecs(:,1)*1,unit_vecs(:,2)*1,unit_vecs(:,3)*1, 200, markerColors,'filled','markerfacealpha',0.5);
hold on
scatter3(0,0,0,400,'filled','markerfacecolor','k')
axis equal

xlabel('X Axis')

[~,idx_maxr2] = max(r2_vecs);
[~,idx_minr2] = min(r2_vecs);

max_vcline = [unit_vecs(idx_maxr2,:)'*-5,unit_vecs(idx_maxr2,:)',unit_vecs(idx_maxr2,:)'*5];
min_vcline = [unit_vecs(idx_minr2,:)'*-5,unit_vecs(idx_minr2,:)',unit_vecs(idx_minr2,:)'*5];
plot3(max_vcline(1,:),max_vcline(2,:),max_vcline(3,:),'g','linewidth',2)
plot3(min_vcline(1,:),min_vcline(2,:),min_vcline(3,:),'r','linewidth',2)

plot3(pc1_contline(1,:),pc1_contline(2,:),pc1_contline(3,:),'k','linewidth',2)
plot3(pc2_contline(1,:),pc2_contline(2,:),pc2_contline(3,:),'--','color','k','linewidth',2)

plot3(pc1_pcline(1,:),pc1_pcline(2,:),pc1_pcline(3,:),'b','linewidth',2)
plot3(pc2_pcline(1,:),pc2_pcline(2,:),pc2_pcline(3,:),'--','color','b','linewidth',2)

% plot3(pc1_acline(1,:),pc1_acline(2,:),pc1_acline(3,:),'r','linewidth',2)
% plot3(pc2_acline(1,:),pc2_acline(2,:),pc2_acline(3,:),'--','color','r','linewidth',2)

% plot3([-2,0,2],[0,0,0],[0,0,0],'k')
% plot3([0,0,0],[-2,0,2],[0,0,0],'k')

% scatter3(cont_mc_qc(:,1,si),cont_mc_qc(:,2,si),cont_mc_qc(:,3,si),'filled','MarkerEdgeColor','k','markerfacecolor','k')
% scatter3(cont_mcRL(:,1),cont_mcRL(:,2),cont_mcRL(:,3),'filled','MarkerEdgeColor','k','markerfacecolor','k')



trimesh(L_stn_cen,'FaceColor','k','EdgeColor','k','FaceAlpha',0.05,'EdgeAlpha',0.1)



%%

% ax2pl = 1;
% fre2pl = ax2pl+3;
% fre2pl = [ax2pl,ax2pl+3];

% mdl_fin = fitlm([pred2pl,AFRE_dub(:,fre2pl)],cont_mcRL(:,ax2pl))
% mdl_fin = fitlm([pred2pl],cont_mcRL(:,ax2pl))


% plot(mdl_fin)



%% Figure __ heatmap showing contact placements in MNI space, color scaled by AFRE in a given axis

cmap = parula(256);
v = rescale(AFRE_m(:,5), 1, 256); 
numValues = length(AFRE_m);
markerColors = zeros(numValues, 3);

% Now assign marker colors according to the value of the data.
for k = 1 : numValues
    row = round(v(k));
    markerColors(k, :) = cmap(row, :);
end

figure

% subplot(1,2,2)
trimesh(R_stn,'FaceColor','k','EdgeColor','k','FaceAlpha',0.05,'EdgeAlpha',0.1)
hold on
scatter3(cont_mni_m(:,1,1),cont_mni_m(:,2,1),cont_mni_m(:,3,1), 100, markerColors,'filled');
% scatter3(cont_mni_cen(1,1,1),cont_mni_cen(1,2,1),cont_mni_cen(1,3,1),200,'filled','MarkerFaceColor',[0, 0, 0],'MarkerEdgeColor',[0, 0, 0],'MarkerFaceAlpha',0.5)
% xlim([-5 5])
% ylim([-5 5])
% zlim([-5 5])
axis equal


% subplot(1,2,1)
trimesh(L_stn,'FaceColor','k','EdgeColor','k','FaceAlpha',0.05,'EdgeAlpha',0.1)
hold on
scatter3(cont_mni_m(:,1,2),cont_mni_m(:,2,2),cont_mni_m(:,3,2), 100, markerColors,'filled');
% scatter3(cont_mni_cen(1,1,1),cont_mni_cen(1,2,1),cont_mni_cen(1,3,1),200,'filled','MarkerFaceColor',[0, 0, 0],'MarkerEdgeColor',[0, 0, 0],'MarkerFaceAlpha',0.5)
% xlim([-5 5])
% ylim([-5 5])
% zlim([-5 5])
axis equal

set(gcf,'color','w');
set(gcf, 'PaperSize', [7 6]);
set(gcf, 'Position',  [100, 100, 1500, 500])



%%

sub_pl = Rat_Mod_QC.subject;

sub_pl_qc = sub_pl(~cont_out);

pred_var = zeros(length(sub_pl_qc),5);
for s = 1:length(sub_pl_qc)
        s_idx = find(PtID == sub_pl_qc(s));

        pred_var(s,1:5) = Demogr(s_idx,[1:4,8]);

end

cont_mcRL = [cont_mc_qc(:,:,1);cont_mc_qc(:,:,2)];

pred_RL = [pred_var;pred_var];
% pred_RL(length(pred_var)+1:end,5) = pred_RL(length(pred_var)+1:end,5)==0;

%%





acpc_dot = [];

acpc_dot = [acpc_dot,sum(cont_mc(:,1:3,1).*AFRE_m(:,1:3),2)];
acpc_dot = [acpc_dot,sum(cont_mc(:,1:3,2).*AFRE_m(:,1:3),2)];
acpc_dot = [acpc_dot,sum(cont_mc(:,1:3,1).*AFRE_m(:,4:6),2)];
acpc_dot = [acpc_dot,sum(cont_mc(:,1:3,2).*AFRE_m(:,4:6),2)];

% figure
histogram(acpc_dot(:,3),'binwidth',2)
hold on
histogram(acpc_dot(:,4),'binwidth',2)
