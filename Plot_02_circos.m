clear;
clc;
rootdir = 'Z:\NN_FU\gPPI_Glasser_yihe\CPMresults\SSRT';
rootfolder = dir(rootdir); rootfolder(1:2,:) = [];
threshold = [1,0.98,0.95];
% plot setting
colmap = '';
clims = [.5 1]; %ignore
opt = 2;
% load('Z:\NN_FU\gPPI_Glasser_yihe\CPMresults\CPM_00_Plotscript\Connectome_cobar\Colorbar.mat')  % colorbar of network path
Glasser = readtable('Z:\NN_FU\gPPI_Glasser_yihe\Glasser_atlas\glasser358.xlsx');  % Glasser label path

for g = 1: length(rootfolder)
   
    %% Step 01: load network module index for BL/FU; make label 
    if strcmp(rootfolder(g).name(1:2),'BL')
        % BL Glasser labels
        netlabel = readtable('Z:\NN_FU\gPPI_Glasser_yihe\CPMresults\GlasserLabel_for_CPM\BL_1.1\Glasser_networklabels.csv');
        load('Z:\NN_FU\gPPI_Glasser_yihe\CPMresults\CPM_00_Plotscript\Connectome_cobar\Colorbar_BL.mat')  % colorbar of network path
    else  
        % FU_glasser network
        netlabel = readtable('Z:\NN_FU\gPPI_Glasser_yihe\CPMresults\GlasserLabel_for_CPM\FU2_1.08\Glasser_networklabels.csv');
        load('Z:\NN_FU\gPPI_Glasser_yihe\CPMresults\CPM_00_Plotscript\Connectome_cobar\Colorbar_FU.mat')  % colorbar of network path
    end

    label = [];set_rl=[];
    j = 1;  num_label = [];
    id_mod_all = [];
    id_all = [];
    set_rl = [];
    flip_id_mod_rl = [];
    id_mod_rl = [];
    sets = [];
    for s = 1:2
        for i = 1:max(netlabel.Network)  
            id_mod = find(netlabel.Network == i & Glasser.Hemisphere == s); 
            id_mod_rl{s,i} = id_mod';
            % generate sets for plot
            id_all = [id_all;id_mod];
            n = length(id_all);
            set_rl{s,i} = [j:n];
            j = n+1;  
%             %label as a number format
%             num_label = [num_label;repmat(i,length(id_mod),1)];
        end
    end

%     sets = [set_rl(1,:),set_rl(2,:)];
%     flip_id_mod_rl = [id_mod_rl(1,:),id_mod_rl(2,:)];
    sets = set_rl(1,:);
    %flip sets for the left hemisphere
    start_id = set_rl{1,size(set_rl,2)}(end)+1;
    left_n_mode = cellfun('length',flip(set_rl(2,:)));
    rest_node = start_id:1:n;  
    a = 1;
    b = 0;
    for ln = 1:size(set_rl,2)
        b = b+left_n_mode(ln);
        set_mod = {rest_node(a:b)};
        sets = [sets,set_mod];
        a = b+1;     
    end
    
    flip_id_mod_rl = [id_mod_rl(1,:),flip(id_mod_rl(2,:))];
    id_mod_all = cell2mat(flip_id_mod_rl(:)');
    
    label = Glasser.sRegionName(id_mod_all); 
    
    for e = 1:n
        num_label{e} = num2str(e);
    end
    
    %% Step 02
    group_dir = dir(fullfile(rootfolder(g).folder,rootfolder(g).name,'CPM*'));  
    for f = 1:length(group_dir)
        
        for t = 1:length(threshold)
            % 3 threshold folds
            matrix_dir =  fullfile(group_dir(f).folder,group_dir(f).name,'plot',['brainnet_',num2str(threshold(t))]);
            fpos_name = fullfile(matrix_dir,['Threshold',num2str(threshold(t)),'_pos_mask.txt']);
            fneg_name = fullfile(matrix_dir,['Threshold',num2str(threshold(t)),'_neg_mask.txt']);
            % Import edge
            pos = readmatrix(fpos_name);
            neg = readmatrix(fneg_name);
            
            % rearrange edges
            pos_mat = pos(id_mod_all,id_mod_all);
            neg_mat = neg(id_mod_all,id_mod_all);
            
            % op_name
            op_pos_name = fullfile(group_dir(f).folder,group_dir(f).name,'plot',['Tif_Pos_connectome_T',num2str(threshold(t)*100)]);
            op_neg_name = fullfile(group_dir(f).folder,group_dir(f).name,'plot',['Tif_Neg_connectome_T',num2str(threshold(t)*100)]);
            
            figure(1)
            h = SchemaBall(pos_mat, label, clims,colmap,opt,sets,colorbar,'#BC4749',2);
%             h = SchemaBall(pos_mat, label, clims,colmap,opt,sets,colorbar,'#BC4749');
            export_fig(gcf,'-tif', '-transparent','-nocrop','-r1200', 'filename',op_pos_name); 
            
            figure(2)
            fig = SchemaBall(neg_mat, label, clims,colmap,opt,sets,colorbar,'#457B9D',2);
            export_fig(gcf,'-tif', '-transparent','-nocrop','-r1200', 'filename',op_neg_name);
            
            % op_name
            op_pos_name_num = fullfile(group_dir(f).folder,group_dir(f).name,'plot',['Num_Pos_connectome_T',num2str(threshold(t)*100)]);
            op_neg_name_num = fullfile(group_dir(f).folder,group_dir(f).name,'plot',['Num_Neg_connectome_T',num2str(threshold(t)*100)]);
            
            
            figure(3)
            hh = SchemaBall(pos_mat, num_label, clims,colmap,opt,sets,colorbar,'#BC4749',2);
            export_fig(gcf,'-tif', '-transparent','-nocrop','-r1200', 'filename',op_pos_name_num); 
            
            figure(4)
            ff = SchemaBall(neg_mat, num_label, clims,colmap,opt,sets,colorbar,'#457B9D',2);
            export_fig(gcf,'-tif', '-transparent','-nocrop','-r1200', 'filename',op_neg_name_num);                      
         

        end
    end
end
