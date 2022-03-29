% set different threhold
function threhold_cpm_network(cpm_path_result,op_dir,freq_thresh)

%% extract network under different threholds, save as a pos/neg_mask.txt(358 regions)
% extract network detailed information under different threholds and save
% as a xlsx document
% cpm_path_result -- the mat file of cpm_results including the path: 
% 'Z:\NN_FU\gPPI_Glasser_yihe\CPM_FU2\GO_RT\BL_GO_RT\output\Relate0.001\GO_RT_CPM.mat'
% op_dir -- output path: 'Z:\NN_FU\gPPI_Glasser_yihe\CPM_FU2\IRV\BL_IRV_717\output'
% freq_thresh -- threshold range:0-1

    %% 01 extract network under different threholds, save as a pos/neg_mask.txt(358 regions)
    if ~exist(op_dir,'dir')
        mkdir(op_dir)
    end
    
    load(cpm_path_result,'parameters');
    pos_mask_all = parameters.pos_mask_all;
    neg_mask_all = parameters.neg_mask_all;

    no_node = size(pos_mask_all,1);
    k_all = size(pos_mask_all,3);
    k = k_all/50;  % 50: iteration time for k fold to calculate mean r value

    [~, ~, pos_edges_thresh, neg_edges_thresh] = ...
        extract_edges_CPM(pos_mask_all, neg_mask_all, no_node, k_all, freq_thresh);

    % Assumes use of Shen atlas parcellation
    pos_edges_orig = pos_edges_thresh(:, [1 3 4]);
    neg_edges_orig = neg_edges_thresh(:, [1 3 4]);

    [pos_edge_mask, neg_edge_mask] = create_masks_CPM(pos_edges_orig, ...
        neg_edges_orig, no_node); 

    pos_mask_file = fullfile(op_dir, ['Threshold',num2str(freq_thresh),'_pos_mask.txt']);
    neg_mask_file = fullfile(op_dir, ['Threshold',num2str(freq_thresh),'_neg_mask.txt']);
    save(pos_mask_file, 'pos_edge_mask', '-ascii');
    save(neg_mask_file, 'neg_edge_mask', '-ascii');
    
    % node_degress
    pos_node_degree = sum(pos_edge_mask,2);
    neg_node_degree = sum(neg_edge_mask,2);
    
     %% extract and save network detailed information under different threholds
    regions = importdata('Z:\NN_FU\gPPI_Glasser_yihe\Glasser_atlas\glasser358.xlsx');
    regions_info = regions.textdata;
    regions_info(1,:) = [];

    for i = 1:size(pos_edges_thresh,1) 
        regions_name_pos{i,1} = cell2mat(regions_info(pos_edges_thresh(i,3),2));
        regions_name_pos{i,2} = cell2mat(regions_info(pos_edges_thresh(i,3),3));
        regions_name_pos{i,3} = cell2mat(regions_info(pos_edges_thresh(i,3),5));

        regions_name_pos{i,4} = cell2mat(regions_info(pos_edges_thresh(i,4),2));  
        regions_name_pos{i,5} = cell2mat(regions_info(pos_edges_thresh(i,4),3));
        regions_name_pos{i,6} = cell2mat(regions_info(pos_edges_thresh(i,4),5));
    end

    % for i = 1:length(edges.neg_edges_thresh) 
    for i = 1:size(neg_edges_thresh,1) 
        regions_name_neg{i,1} = cell2mat(regions_info(neg_edges_thresh(i,3),2));
        regions_name_neg{i,2} = cell2mat(regions_info(neg_edges_thresh(i,3),3));
        regions_name_neg{i,3} = cell2mat(regions_info(neg_edges_thresh(i,3),5));

        regions_name_neg{i,4} = cell2mat(regions_info(neg_edges_thresh(i,4),2));  
        regions_name_neg{i,5} = cell2mat(regions_info(neg_edges_thresh(i,4),3));
        regions_name_neg{i,6} = cell2mat(regions_info(neg_edges_thresh(i,4),5));
    end

    varibalename = {'region1_fullname','R/L_1','region1','region2_fullname','R/L_2','region2'};
    t_pos = cell2table(regions_name_pos,'VariableNames',varibalename);
    t_neg = cell2table(regions_name_neg,'VariableNames',varibalename);
    excelfile = fullfile(op_dir,['Threshold',num2str(freq_thresh),'CPM_edge_regions.xlsx']);
    writetable(t_pos,excelfile,'Sheet',['pos_mask regions ', num2str(k), 'k Success']);
    writetable(t_neg,excelfile,'Sheet',['neg_mask regions ', num2str(k), 'k Success']);
    
    
    %% 02 plot pos/neg network for cpm_results under different thresholds by using BrainnetViewer
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if size(pos_edge_mask,1) ~= 358
        error('A wrong matrix was imported')
    end

    pos_region = pos_edges_thresh(:,3:4);
    pos_region = pos_region(:);
    pos_id = unique(pos_region);

    neg_region = neg_edges_thresh(:,3:4);
    neg_region = neg_region(:);
    neg_id = unique(neg_region);
    
    remain_matrix_pos = pos_edge_mask(pos_id,pos_id);
    remain_matrix_neg = neg_edge_mask(neg_id,neg_id);
    
    % glasser 358 information
    glasser_info = readtable('Z:\NN_FU\gPPI_Glasser_yihe\Glasser_atlas\glasser358.txt');  % glasser358 atlas
    % update yihe
    glasser_info_pos = glasser_info;
    glasser_info_pos.Var5 = pos_node_degree;
    glasser_info_neg = glasser_info;
    glasser_info_neg.Var5 = neg_node_degree;    
    
    
    neg_T_region = glasser_info_neg(neg_id,:);
    pos_T_region = glasser_info_pos(pos_id,:);
    
%     neg_T_region = glasser_info(neg_id,:);
%     pos_T_region = glasser_info(pos_id,:);

    %% HERE ADD NODE!!

    % make edge file for brainnet viewer
    if ~exist(fullfile(op_dir,'pos_net.edge'),'file')       
        writematrix(remain_matrix_pos,fullfile(op_dir,'pos_net.txt'),'Delimiter','tab');
        cd(op_dir)
        eval(['!rename ',' pos_net.txt',' pos_net.edge']); 
    end 

    if ~exist(fullfile(op_dir,'neg_net.edge'),'file')   
        writematrix(remain_matrix_neg,fullfile(op_dir,'neg_net.txt'),'Delimiter','tab');
        cd(op_dir)
        eval(['!rename ',' neg_net.txt',' neg_net.edge']); 
    end

    % make node file for brainnet viewer
    if ~exist(fullfile(op_dir,'pos_node.node'),'file')   
        writetable(pos_T_region,fullfile(op_dir,'pos_node.txt'),'Delimiter','tab','WriteVariableNames',0);
        eval(['!rename ',' pos_node.txt',' pos_node.node']); 
    end

    if ~exist(fullfile(op_dir,'neg_node.node'),'file')   
        writetable(neg_T_region,fullfile(op_dir,'neg_node.txt'),'Delimiter','tab','WriteVariableNames',0);
        eval(['!rename ',' neg_node.txt',' neg_node.node']); 
    end
    %% 03 briannet viewer
    template = 'C:\software\BrainNetViewer_20191031\Data\SurfTemplate\BrainMesh_ICBM152.nv';

    % pos network setting
    pos_node_file = fullfile(op_dir,'pos_node.node');
    pos_edge_file = fullfile(op_dir,'pos_net.edge');
    pos_cfg = 'Z:\NN_FU\gPPI_Glasser_yihe\CPMresults\Connectviewer_cfg\pos_cfg.mat';
    pos_sav_name =  fullfile(op_dir,['Threshold',num2str(freq_thresh),'_pos_brainnet.jpg']);
    % neg network setting
    neg_node_file = fullfile(op_dir,'neg_node.node');
    neg_edge_file = fullfile(op_dir,'neg_net.edge');
    neg_cfg = 'Z:\NN_FU\gPPI_Glasser_yihe\CPMresults\Connectviewer_cfg\neg_cfg.mat';
    neg_sav_name =  fullfile(op_dir,['Threshold',num2str(freq_thresh),'neg_brainnet.jpg']);

    BrainNet_MapCfg(template,pos_node_file,pos_edge_file,pos_cfg,pos_sav_name);
    BrainNet_MapCfg(template,neg_node_file,neg_edge_file,neg_cfg,neg_sav_name);
    
    
end