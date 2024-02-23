%% AC/PC comarison - b/w lead DBS and Raters
clear all
clc

%% Load files
main_dir = 'C:\Research\DBS_Retro\Final';
dir_dat = fullfile(main_dir,'data');
dir_out = fullfile(main_dir,'Output');

% Load contact coors in MNI space
load(fullfile(dir_dat,'input_acpc_MNI_leaddbs','dbs_cont_mni.mat'));

% Load transformed AFID coors using Lead-DBS transform
load(fullfile(dir_dat,'input_fid_MNI_leaddbs','afids_mni.mat'));

% for labels (note, add afid_lab)
afid_lab = readtable(fullfile(main_dir,'etc','afid_lab.xlsx'));
afid_lab = afid_lab{:,1};

% ground truth AFID placement in MNI09basym
load(fullfile(dir_dat,'fid_standards','MNI152NLin2009bAsym_rater_standard','MNI152NLin2009bAsym_desc-raterstandard_afids.mat'))
gold2use = mni_rater_standard;

% Load STL meshes for plots
brain_stl = stlread(fullfile(main_dir,'etc','Segmentation_brain.stl'));
L_stn = stlread(fullfile(main_dir,'etc','STN_Segment_L.stl'));
R_stn = stlread(fullfile(main_dir,'etc','STN_Segment_R.stl'));

%%
% Determine which transformed FID to use (6 = mean AFID coor transformed,
% i.e. consensus AFRE)
fid2use = squeeze(fcormni(:,:,:,6,:));

% Compute AFRE for each AFID
fid_diff = fid2use - repmat(gold2use(:,2:4),1,1,24,2);
fid_AFRE = squeeze(sqrt(sum(fid_diff.^2,2)));

% Each rater has a unique transform, AFRE for each subject was computed as
% the mean AFRE across the two raters
AFRE_ldbs = squeeze(mean(fid_AFRE,3));

% Global AFRE computed as final row (33)
AFRE_ldbs(end+1,:) = mean(AFRE_ldbs,1);

% Contact coords in mni for subset with AFIDs, and compute displacement
% from the mean contact position
cont_fids = dbs_cont_mni(fidsinpd,2:4,:,:);
cont_err = cont_fids - mean(cont_fids,1);

%% Reorg data - Compute mean AFRE and contact displacement across raters, and combine R/L contacts

% Generate fidcat_m: x,y,z AFRE for each AFID
fid_reorg = permute(fid_diff,[3,2,1,4]);
fid_reorgm = mean(fid_reorg,4);
fidcat_m = [fid_reorgm;fid_reorgm];

% Generate fidcat_m: x,y,z contact displacement for R and L contacts
cont_err_m = mean(cont_err,4);
cont_catm = [cont_err_m(:,:,1);cont_err_m(:,:,2)];

%% Simple linear regression for each AFID and axis

r2 = nan(32,3);
p = nan(32,3);

for afid = 1:32
    for ax = 1:3
        mdl = fitlm(fidcat_m(:,ax,afid),cont_catm(:,ax));
    
        r2(afid,ax) = mdl.Rsquared.Ordinary;
        p(afid,ax) = mdl.Coefficients.pValue(2);
    end
end

%% Figure 5 - Heatmap
mkdir(fullfile(dir_out,'Figures','Fig5'));

fidlabpl = strcat(afid_lab, ' -  ',cellfun(@(x) num2str(x), num2cell(1:32),'UniformOutput',false)');

figure
heatmap(round(r2,2),'CellLabelColor','none','YData',fidlabpl,'XData',{'X','Y','Z'})
set(gcf,'color','w','Position',  [100, 100, 350, 1200],'PaperSize', [6 13]);
saveas(gcf,fullfile(dir_out,'Figures','Fig5','Fig5_heatmap.pdf'));

%% Figure 5 - brain mesh
cmap(1:256,1) = 0.05:0.95/256:1-0.05/256;
cmap(1:256,2) = 0.45:0.55/256:1-0.45/256;
cmap(1:256,3) = 1;
cmap = flipud(cmap);

v = rescale(r2(:,2), 1, 256); 
numValues = length(r2);
markerColors = zeros(numValues, 3);

% Now assign marker colors according to the value of the data.
for k = 1 : numValues
    row = round(v(k));
    markerColors(k, :) = cmap(row, :);
end

figure
trimesh(brain_stl,'FaceColor','k','EdgeColor','none','FaceAlpha',0.03)
axis equal

hold on

trimesh(L_stn,'FaceColor','k','EdgeColor','none','FaceAlpha',0.3)
trimesh(R_stn,'FaceColor','k','EdgeColor','none','FaceAlpha',0.3)

scatter3(gold2use(:,2),gold2use(:,3),gold2use(:,4), 200, markerColors,'filled','MarkerEdgeColor','k');

zlim([-60 50])
ylim([-100 40])
view(220,20)

xlabel('X-Axis (mm)')
ylabel('Y-Axis (mm)')
zlabel('Z-Axis (mm)')

set(gcf,'color','w','Position',  [100, 100, 1500, 1200],'PaperSize', [16 13]);
% set(gcf,'color','w','Renderer','Painter','Position',  [100, 100, 1500, 1200],'PaperSize', [16 13]);

saveas(gcf,fullfile(dir_out,'Figures','Fig5','Fig5_brainmesh.pdf'));

%% Figure 6a - 3d scatter plot of select AFIDs (2,3,4,14)
mkdir(fullfile(dir_out,'Figures','Fig6'));

c_ord = colororder;

cont_fidsm = [squeeze(mean(cont_fids(:,:,1,:),4));squeeze(mean(cont_fids(:,:,2,:),4))];

figure
trimesh(L_stn,'FaceColor','k','EdgeColor','k','FaceAlpha',0.03,'edgealpha',0.05)
hold on
trimesh(R_stn,'FaceColor','k','EdgeColor','k','FaceAlpha',0.03,'edgealpha',0.05)
axis equal

count = 1;
for afid = [2,3,14,4]
    fidcat = [squeeze(fid2use(afid,:,:,1)),squeeze(fid2use(afid,:,:,2))];
    scatter3(fidcat(1,:),fidcat(2,:),fidcat(3,:),'MarkerEdgeColor',c_ord(count,:),'MarkerFaceColor',c_ord(count,:),'markeredgealpha',0.6,'markerfacealpha',0.4)
    scatter3(gold2use(afid,2),gold2use(afid,3),gold2use(afid,4),200,'filled','MarkerEdgeColor','k','MarkerFaceColor',c_ord(count,:))
    count = count+1;
end

scatter3(cont_fidsm(:,1),cont_fidsm(:,2),cont_fidsm(:,3),'MarkerEdgeColor','k','MarkerFaceColor','k','markeredgealpha',0.6,'markerfacealpha',0.4)

% xlabel('X-Axis')
% ylabel('Y-Axis')
% zlabel('Z-Axis')

view(75,15)
set(gcf,'color','w','Position',  [100, 100, 800, 600],'PaperSize', [9 7]);
saveas(gcf,fullfile(dir_out,'Figures','Fig6','Fig6a.pdf'));

%%
fid2pl = [2,3,14,4];

fid_mcoor = squeeze(mean(fid2use,[3,4]));
fid_dev = fid2use - repmat(fid_mcoor,1,1,24,2);
fid_devm = mean(fid_dev,4);

fid_devcat = [];
fids_multidim = [];
fid_devpcs = [];
clear coeff_fid_ind
clear exp_fid_indiv
for afid = 1:length(fid2pl)
    temp_fidev = squeeze(fid_devm(fid2pl(afid),:,:))';
    fid_devcat = [fid_devcat;temp_fidev];
    fids_multidim = [fids_multidim,temp_fidev];
    [coeff_fid_ind(:,:,afid),~,~,~,exp_fid_indiv(:,:,afid)] = pca(temp_fidev);
    fid_devpcs = [fid_devpcs,temp_fidev*coeff_fid_ind(:,:,afid)];
end

[coeff_cont,~,~,~,exp_cont] = pca(cont_catm);
cont_catpc = cont_catm*coeff_cont;

[coeff_fmult,~,~,~,exp_fmult] = pca(fids_multidim);
fid_mdimpc = fids_multidim*coeff_fmult;

%% Figure 6b, Correlation matricies for AFREs, along Cartesian and PC axes

% Correlation matrix for AFREs along Cartesian axes
[fids_multicor,p_multcor] = corr(fids_multidim);
fids_multicor = fids_multicor.^2;

% Correlation matrix for AFREs along PCs
[fipc_multicor,p_mpcor] = corr(fid_devpcs);
fipc_multicor = fipc_multicor.^2;

% Triangular matrix, as shown in figure 6c
fids_corcomb = fids_multicor;
fids_corcomb(tril(ones(12))>0) = fipc_multicor(tril(ones(12))>0);

% for 2 separate matricies
% figure
% subplot(1,2,1)
% heatmap(fids_multicor,'ColorLimits',[0 1])
% subplot(1,2,2)
% heatmap(fipc_multicor,'ColorLimits',[0 1])
% set(gcf,'color','w','Position',  [100, 100, 1600, 500],'PaperSize', [16 13]);

figure
heatmap(fids_corcomb,'ColorLimits',[0 1])
set(gcf,'color','w','Position',  [100, 100, 600, 500],'PaperSize', [7 6]);
saveas(gcf,fullfile(dir_out,'Figures','Fig6','Fig6b.pdf'));

%% Figure 6c, matrix showing total variance explained of each feature for top 4 PCs

cmap = parula(256);
v = rescale([fid_mdimpc(:,1);fid_mdimpc(:,1)], 1, 256); 
numValues = length([fid_mdimpc(:,1);fid_mdimpc(:,1)]);
markerColors = zeros(numValues, 3);

% Now assign marker colors according to the value of the data.
for k = 1 : numValues
    row = round(v(k));
    markerColors(k, :) = cmap(row, :);
end
cont_fidsm = [squeeze(mean(cont_fids(:,:,1,:),4));squeeze(mean(cont_fids(:,:,2,:),4))];

figure
trimesh(L_stn,'FaceColor','k','EdgeColor','k','FaceAlpha',0.03,'edgealpha',0.01)
hold on
trimesh(R_stn,'FaceColor','k','EdgeColor','k','FaceAlpha',0.03,'edgealpha',0.02)
axis equal

for afid = [2,3,14,4]
    count = 1;
    fidcat = [squeeze(fid2use(afid,:,:,1)),squeeze(fid2use(afid,:,:,2))];

    for i = 1:24
        scatter3(fidcat(1,i),fidcat(2,i),fidcat(3,i),50,'MarkerEdgeColor',markerColors(i,:),'MarkerFaceColor',markerColors(i,:),'markeredgealpha',0.6,'markerfacealpha',0.6)
        count = count+1;
    end

    scatter3(gold2use(afid,2),gold2use(afid,3),gold2use(afid,4),100,'filled','MarkerEdgeColor','k','MarkerFaceColor','k')
end

for i = 1:48
    scatter3(cont_fidsm(i,1),cont_fidsm(i,2),cont_fidsm(i,3),'MarkerEdgeColor',markerColors(i,:),'MarkerFaceColor',markerColors(i,:),'markeredgealpha',0.6,'markerfacealpha',0.6)
end

colorbar
caxis([min(fid_mdimpc(:,1)) max(fid_mdimpc(:,1))])

% xlabel('X-Axis')
% ylabel('Y-Axis')
% zlabel('Z-Axis')

view(75,15)

set(gcf,'color','w','Position',  [100, 100, 800, 600],'PaperSize', [9 7]);
saveas(gcf,fullfile(dir_out,'Figures','Fig6','Fig6c.pdf'));

%% Figure 6d, correlating AFREs along PCs to electrode displacement along PCs

fidpc_contcor = corr([fid_mdimpc(:,1:4);fid_mdimpc(:,1:4)],cont_catpc).^2;

figure
subplot(3,1,1)
heatmap(round(fidpc_contcor(1:4,:),3)','ColorLimits',[0 0.3])
xlabel('AFRE PrCs')
ylabel('Electrode Disp PrCs')

subplot(3,1,2:3)
mdl = fitlm([fid_mdimpc(:,1);fid_mdimpc(:,1)],cont_catpc(:,1));
h = plot(mdl);
dataHandle = findobj(h,'DisplayName','Data');
fitHandle = findobj(h,'DisplayName','Fit');
cbHandles = findobj(h,'DisplayName','Confidence bounds');
cbHandles = findobj(h,'LineStyle',cbHandles.LineStyle, 'Color', cbHandles.Color);
dataHandle.Color = [0 0.2 0.6];
fitHandle.Color = [0 0 0];
set(cbHandles, 'Color', [0,0,0])

hold on
scatter([fid_mdimpc(:,1);fid_mdimpc(:,1)],cont_catpc(:,1),50,markerColors,'filled')
axis equal


xlabel('AFRE (PrC1)')
ylabel('Electrode Displacement (PrC1)')
title([])
legend('off')

set(gcf,'color','w','Position',  [100, 100, 600, 900],'PaperSize', [7 10]);
saveas(gcf,fullfile(dir_out,'Figures','Fig6','Fig6d.pdf'));

%%

mdl_fit = mdl.Fitted;
mdl_exp_mm = prctile(abs(mdl_fit),[25,50,75]);
