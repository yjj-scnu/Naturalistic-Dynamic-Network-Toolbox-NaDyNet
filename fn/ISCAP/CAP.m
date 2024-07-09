clc
clear
rootdir = 'G:\小王子数据\derivation\FR\Runs';
cd(rootdir)
storyList = dir('*run*');
for stotyIndex = 1:9

    runname = storyList(stotyIndex).name
    workingDir = [rootdir filesep runname '\LOO_ResReg'];
    desdir = [workingDir filesep 'ISCAP'];
    if ~exist(desdir, 'dir')
        mkdir(desdir)
    end
    cd(workingDir)
    prefix = 'sub';
    TR = 2;
    Tmot = 0.5;
    % Editing of all the CAP generation parameters
    K = 4;
    n_rep = 20;
    Pp = 100;
    Pn = 100;
    % Resets the CAP parameters (CAPs, standard deviation within CAPs and
    % indices of the CAPs to which all retained frames were assigned)
    CAP = [];
    STDCAP = [];

    
    if prefix == '*'
        subList = dir('*');
    else
        subList = dir([prefix '*']);
    end
    % loading hrd
    cd([workingDir filesep subList(1).name])
    firstNIIFile = dir('*.nii');

    allHead = spm_vol([workingDir filesep subList(1).name filesep firstNIIFile(1).name]);
    brain_info = {};
    brain_info{1}=allHead(1);


    % loading mask
    [pathstr, name, ext]=fileparts(which('CAP_TB.m'));
    a = spm_vol(fullfile(pathstr,'DefaultData','Default_mask.nii'));
    b = spm_read_vols(a);
    b(b < 0.9) = 0;
    b(b >= 0.9) = 1;
    maskf = CAP_V2V(b,a.dim,a.mat,brain_info{1}.dim,brain_info{1}.mat);
    mask = {};
    mask{1} = logical(maskf(:));

    %% loading subject and FD
    TC = {};
    FD = {};
    for i = 1:size(subList,1)

        disp(['Currently loading subject ',num2str(i),'...']);
        cd([workingDir filesep subList(i).name])
        NIIFile = dir('*.nii');

        tmp_data = [];
        [d,h]=y_Read(NIIFile(1).name);
        [s1,s2,s3,s4]=size(d);
        temp=reshape(d,[s1*s2*s3,s4]);
        tmp_data=temp(mask{1},:)';

        % Z-scoring is performed within the toolbox
        % tmp_data = detrend(tmp_data);
        % tmp_data = zscore(tmp_data);
        tmp_data = (tmp_data-repmat(mean(tmp_data),size(tmp_data,1),1)) ./ repmat(std(tmp_data),size(tmp_data,1),1);
        tmp_data(isnan(tmp_data)) = 0;

        % The ready-to-analyse data is put in TC

        TC{1}{i} = tmp_data;
        TXTFile = dir('*.txt');

        if size(TXTFile, 1) == 0
            FD{1}(:,i) = zeros(size(tmp_data, 1),1);
            disp(['Could not process motion text file of ' subList(i).name '; assuming zero movement...']);
        else
            FD{1}(:,i) = CAP_ComputeFD(TXTFile(1).name);
        end
    end
    
    SubjSize = {};
    SubjSize.VOX = size(TC{1}{1},2);
    SubjSize.TP = size(TC{1}{1},1);
    %         % Sets the text label about data dimensions
    %         set(Dimensionality_Text, 'String', [num2str(SubjSize.TP),...
    %             ' frames x ',num2str(SubjSize.VOX),' voxels (',...
    %             strjoin(arrayfun(@(x) num2str(x),cell2mat(n_subjects),...
    %             'UniformOutput',false),'+'),')']);
    %     Underlay_info = {};
    brain = importdata('brain.mat');
    Underlay = load_nii('Underlay.nii');

    Underlay_dim=[];Underlay_mat=[];Underlay_info=[];

    Underlay_mat = [Underlay.hdr.hist.srow_x; Underlay.hdr.hist.srow_y; Underlay.hdr.hist.srow_z; 0 0 0 1];
    Underlay_dim = Underlay.hdr.dime.dim;
    Underlay_dim = Underlay_dim(2:4);
    Underlay_info.dim = Underlay_dim;
    Underlay_info.mat = Underlay_mat;
    brain = CAP_V2V(brain,Underlay_info.dim,...
        Underlay_info.mat,brain_info{1}.dim,brain_info{1}.mat);

    %% Seed-free analysis
    % Performs the analysis to extract frames of activity
    % Xonp and Xonn contain the frames (deactivation frames have been
    % switched in sign, so that deactivation is positive)
    % Xonp 保存剩下的帧，大小为 1 * nsub的cell数组，每个cell里面是剩下的TR * nvoxel
    % p是3*nsub的double，第一第二行是相同的，提剔除的时间帧占比
    % 第三行 = 1-第一行
    % Indices 1*1的结构体，3个字段。 其中srubbed是nt * nsub  的logic，记录某个人被剔除的
    % 的时间点,Indices.scrubbed==Indices.scrubbedandactive
    % Indices.kept.active = ~Indices.scrubbed
    Xonp = {};
    [Xonp{1},p,Indices] = CAP_find_activity_SeedFree(TC{1},...
        FD{1},Tmot);

    % Percentage of retained frames across subjects
    RetainedPercentage = {};
    RetainedPercentage{1} = p(3,:);

    % Indices of the frames that have been retained (used later for metrics
    % computations)
    FrameIndices = {};
    FrameIndices{1} = Indices;
    %% running cluster

    % Indices of the CAP to which frames from the reference population and from
    % the other populations are assigned
    idx = {};
    [CAP,Disp,STDCAP,idx{1},...
        CorrDist,sfrac] = Run_Clustering(cell2mat(Xonp{1}),K,mask{1},brain_info{1},...
        Pp,Pn,n_rep,1,'SeedFree');

    %% computeMetric
    TPM = {};
    Counts = {}; Number = {}; Avg_Duration = {}; Duration = {}; TM={};
    From_Baseline = {}; To_Baseline = {}; Baseline_resilience = {};Resilience={};
    Betweenness = {}; kin = {}; kout = {}; SubjectEntries = {};
    [TPM{1},Counts{1},Number{1},Avg_Duration{1},...
        Duration{1},TM{1},From_Baseline{1},To_Baseline{1},...
        Baseline_resilience{1},Resilience{1},...
        Betweenness{1},kin{1},kout{1},SubjectEntries{1}] =...
        Compute_Metrics_simpler(idx{1},FrameIndices{1}.kept.active,...
        FrameIndices{1}.scrubbedandactive,...
        K,TR);


    % Concatenates information from the different datasets
    tmp_toplot = [];
    tmp_toplot = [tmp_toplot; TPM{1}; 0*ones(5,SubjSize.TP)];
    tmp_toplot = tmp_toplot(1:end-5,:);


    figure(1)
    imagesc(tmp_toplot);
    custom_cm = cbrewer('qual','Set1',4);
    custom_cm = [0.05,0.05,0.05;1,1,1;custom_cm];
    colormap((custom_cm));
    set(gca, 'FontName','Arial','FontSize',25,'LineWidth', 1.5);
    xlim([0 size(tmp_toplot, 2)])
    xlabel(gca,'Time [s]','FontSize',36);
    ylabel(gca,'Subjects','FontSize',36);
    clim([-1,4+1]);
    set(gcf,'Position',[100 100 1920*0.6 1080*0.4]);

    jianju = floor((size(tmp_toplot, 2) / 4));
    xticks([0:jianju:(size(tmp_toplot, 2)  - jianju), floor(size(tmp_toplot, 2))])
    xticklabels(floor([0:jianju:(size(tmp_toplot, 2)  - jianju), floor(size(tmp_toplot, 2))].*TR));


    if size(tmp_toplot, 1) < 4
        yticks(0:1:size(tmp_toplot, 1));
    end

    set(gca,'tickdir','in');

    filename=[desdir filesep runname '_stateTransition'];
    print(1,'-dtiff','-r300',filename);
    close(1)
    cd(desdir)
    save([desdir filesep runname '_StateTransition.mat'], "tmp_toplot");

    %%
    % Saves NIFTI files storing the CAPs in MNI space
    CAPToNIFTI(CAP,...
        mask{1},brain_info{1},...
        desdir,['CAP_NIFTI_',runname]);

    CAPToNIFTI(CAP_Zscore(CAP),...
        mask{1},brain_info{1},...
        desdir,['CAP_NIFTI_ZScored_',runname]);
    %% brainnet viewer
    cd(desdir)
    sublist=dir('*.nii');
    [pathstr, name, ext]=fileparts(which('BrainNet.m'));
    surfile=[pathstr filesep '\Data\SurfTemplate\BrainMesh_ICBM152_smoothed.nv'];
    cfgfile=[pathstr filesep '\brainnet\pro\cfg_brainnet_subisc.mat'];
    for f=1:length(sublist)
        niifile=[desdir filesep sublist(f).name];
        picname=[desdir filesep sublist(f).name(1:end-4) '.tif'];
        Yuan_BrainNet_MapCfg(surfile,cfgfile,niifile,picname);
    end
    %%
    clear TC;
    clear FD;
end
