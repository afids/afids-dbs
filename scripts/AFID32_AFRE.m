%% AFID comarison - b/w lead DBS and Raters
clear all 
clc
    
%% Load files
main_dir = 'C:\Research\DBS_Retro\Final';
dir_dat = fullfile(main_dir,'data');
dir_out = fullfile(main_dir,'Output');

% Load AFLE computed for 32 AFIDs from Abbass et al., 2022
load(fullfile(dir_dat,'input_fid_native','AFLE_native.mat'));
% Load transformed AFID Coors from Abbass et al., 2022
load(fullfile(dir_dat,'input_fid_MNI_fmriprep','AFID_Sum.mat'));
% Load computed consensus AFRE from Abbass et al., 2022
load(fullfile(dir_dat,'input_fid_MNI_fmriprep','AFRE_con.mat'));
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
%% Compute AFREs

% Determine which transformed FID to use (6 = mean AFID coor transformed,
% i.e. consensus AFRE)
fid2use = squeeze(fcormni(:,:,:,6,:));

% Compute AFRE for each AFID
fid_diff = fid2use - repmat(gold2use(:,2:4),1,1,24,2);
fid_AFRE = squeeze(sqrt(sum(fid_diff.^2,2)));

% Find common patients from current study used in Abbass et al. 2022
pd_pts = squeeze(Tot_Data(1,5,:,1));
pd_common = ismember(pd_pts,subs2use);

% AFLE for common patients
mAFLE_nat = mean(AFLE_native(:,pd_common,:),3);

% AFRE using fmriprep
AFRE_fmri = AFRE_con(:,pd_common);

% Each rater has a unique transform, AFRE for each subject was computed as
% the mean AFRE across the two raters
AFRE_ldbs = squeeze(mean(fid_AFRE,3));

% Global AFRE computed as final row (33)
AFRE_fmri(end+1,:) = mean(AFRE_fmri,1);
AFRE_ldbs(end+1,:) = mean(AFRE_ldbs,1);

%% Statistics

% Statistics for AFLE and AFRE using fmriprep and lead-dbs
AFLE_stats = compute_stats_sum(mAFLE_nat');
AFRE_fmri_stats = compute_stats_sum(AFRE_fmri');
AFRE_ldbs_stats = compute_stats_sum(AFRE_ldbs');

% Sign rank between AFREs obtained using fmriprep and lead-dbs
clear AFRE_sr
for afid = 1:33
  AFRE_sr(afid,1) = signrank(AFRE_fmri(afid,:),AFRE_ldbs(afid,:));
end

%% Figure 3a - boxchart of AFREs
mkdir(fullfile(dir_out,'Figures','Fig3'));

cord_fid = colormap(jet(32));
close Figure 1

fidlabpl = strcat(afid_lab, ' -  ',cellfun(@(x) num2str(x), num2cell(1:32),'UniformOutput',false)');
x = zeros(24,1);
fidplord = fliplr(1:32);

for afid = 1:32
    boxchart(x+fidplord(afid),AFRE_ldbs(afid,:),'BoxFaceColor',cord_fid(afid,:),'MarkerColor',cord_fid(afid,:),'BoxWidth',0.6,'orientation','horizontal')
    hold on
end

yticks(1:32)
yticklabels(flipud(fidlabpl))
ylim([0.5 32.5])

xlabel('Fiducial Registration Error (mm)')
xticks([0 2.5 5 10 15])
xline([0 2.5 5 10 15],'--')

set(gcf,'color','w','Position',  [100, 100, 500, 800],'PaperSize', [6 9]);
saveas(gcf,fullfile(dir_out,'Figures','Fig3','Fig3a.pdf'));

%% Figure 3b - 3D mesh with scatter plot of registered fiducials
figure

trimesh(brain_stl,'FaceColor','k','EdgeColor','none','FaceAlpha',0.03)
axis equal
hold on
for afid = 1:32

    fidcat = [squeeze(fid2use(afid,:,:,1)),squeeze(fid2use(afid,:,:,1))];
    scatter3(fidcat(1,:),fidcat(2,:),fidcat(3,:),'MarkerfaceColor',cord_fid(afid,:),'markerfacealpha',0.2,'markeredgecolor',cord_fid(afid,:))
    scatter3(gold2use(afid,2),gold2use(afid,3),gold2use(afid,4),100,'filled','MarkerEdgeColor','k','MarkerFaceColor',cord_fid(afid,:))
end

% xlim([-40 40])
% zlim([-40 50])
% ylim([-90 40])
view(220,10)

xlabel('X-Axis (mm)')
ylabel('Y-Axis (mm)')
zlabel('Z-Axis (mm)')
% grid off
set(gcf,'color','w','Position',  [0, 0, 1500, 1200],'PaperSize', [16 13]);

% set(gcf,'color','w','Renderer','Painter','Position',  [0, 0, 1500, 1200],'PaperSize', [16 13]);
saveas(gcf,fullfile(dir_out,'Figures','Fig3','Fig3b.pdf'));

%% Figure S2a - Boxplot comparing AFREs between fmriprep and lead-dbs
mkdir(fullfile(dir_out,'Figures','FigS2'));

x = zeros(24,1);
fidplord = fliplr(1:32);

c_ord = colororder;

figure
hold on

for afid = 1:32
    boxchart(x+fidplord(afid)-0.2,AFRE_ldbs(afid,:),'BoxFaceColor',c_ord(1,:),'MarkerColor',c_ord(1,:),'BoxWidth',0.35,'orientation','horizontal')
    boxchart(x+fidplord(afid)+0.2,AFRE_fmri(afid,:),'BoxFaceColor',c_ord(2,:),'MarkerColor',c_ord(2,:),'BoxWidth',0.35,'orientation','horizontal')
end

xline([5,10],'--')
yline(1:32,'--')

xlabel('Fiducial Registration Error (mm)')
yticks(1:32)
yticklabels(flipud(fidlabpl))
set(gcf,'color','w','Position',  [100, 100, 600, 1200],'PaperSize', [7 13]);

saveas(gcf,fullfile(dir_out,'Figures','FigS2','FigS2a.pdf'));
%% Figure S2b

figure
trimesh(brain_stl,'FaceColor','k','EdgeColor','none','FaceAlpha',0.03)
axis equal
hold on
trimesh(L_stn,'FaceColor','k','EdgeColor','k','FaceAlpha',0.3,'edgealpha',0.05)
trimesh(R_stn,'FaceColor','k','EdgeColor','k','FaceAlpha',0.3,'edgealpha',0.05)

alpha = 0.05/32;

AFRE_sig = AFRE_sr<alpha;
idx_fmri_hi = and(AFRE_sig(1:32),AFRE_ldbs_stats(1:32,3)<AFRE_fmri_stats(1:32,3));
idx_ldbs_hi = and(AFRE_sig(1:32),AFRE_ldbs_stats(1:32,3)>AFRE_fmri_stats(1:32,3));
idx_ns = ~or(idx_fmri_hi,idx_ldbs_hi);

scatter3(gold2use(idx_ns,2),gold2use(idx_ns,3),gold2use(idx_ns,4),400,'filled','MarkerFaceColor','k')
scatter3(gold2use(idx_ldbs_hi,2),gold2use(idx_ldbs_hi,3),gold2use(idx_ldbs_hi,4),400,'filled','MarkerFaceColor',c_ord(2,:),'markeredgecolor','k','markerfacealpha',0.5)
scatter3(gold2use(idx_fmri_hi,2),gold2use(idx_fmri_hi,3),gold2use(idx_fmri_hi,4),400,'filled','MarkerFaceColor',c_ord(1,:),'markeredgecolor','k','markerfacealpha',0.5)

view(230,20)
set(gcf,'color','w','Position',  [100, 100, 1500, 1200],'PaperSize', [16 13]);

saveas(gcf,fullfile(dir_out,'Figures','FigS2','FigS2b.pdf'));
