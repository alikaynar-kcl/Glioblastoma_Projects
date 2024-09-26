

cd  '{PATH}/Human-GEM'
load('Human-GEM.mat');  % loads model as structure named "ihuman"
%  model = ravenCobraWrapper(ihuman); %  Converting RAVEN structure to
%  COBRA structure 

ihuman = addBoundaryMets(ihuman); % Preparing model for tINIT, closing model
essentialTasks = parseTaskList('{PATH}/Human-GEM/data/metabolicTasks/metabolicTasks_Essential.txt');
%checkTasks(ihuman, [], true, false, false, essentialTasks);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Generating context-specific model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Transcriptome data
TCGA_tpm  = readtable('{PATH}/Model_input_T.txt'); %mean value of genes sppecific to patient group 
TCGA_tpm.Properties.VariableNames = ["ensemble_id","NT","TR","TP","Low TPS","Mid TPS","High TPS"];
%% Model generation 

[~, n] = size(TCGA_tpm);
numSamp = n-1;
Models_T = cell(numSamp, 1);
tic
for i = 1:numSamp
data_struct.genes = cellstr(TCGA_tpm{:, 1});
data_struct.tissues =  TCGA_tpm.Properties.VariableNames(i+1);
data_struct.levels = TCGA_tpm{:, i+1};
data_struct.threshold = 1;
Models_T{i} = getINITModel2(ihuman, data_struct.tissues{1}, [],[], data_struct, [], true, [], true, true, essentialTasks, [],[]);
Models_T{i}.id = data_struct.tissues{1};
end
toc

save('Models_T.mat', 'Models_T')% Models_T contain 6 model. But 2 of them (TR and mid_tps) are not considered in this study. 
% Checking each model for Essential Tasks (build-in file). Each model need
% to pass essential task
checkTasks(Models_T{6,1}, [], true, false, false, essentialTasks);

%% For structuraly model comparison 
%prepData = prepHumanModelForftINIT(ihuman, false, '/Users/alikaynar/Desktop/UpgrdPostDoc/GEM/GBM_GEM/Human-GEM/data/metabolicTasks/metabolicTasks_Essential.txt', '/Users/alikaynar/Desktop/UpgrdPostDoc/GEM/GBM_GEM/Human-GEM/model/reactions.tsv');
%save('prepData.mat', 'prepData') % prepData is useful for structural
%comparison
load('{PATH}/Human-GEM/prepData.mat')
%load('prepData.mat')
load('{PATH}/Models_T.mat')
for i=1:6
Models_T{i,1}.id = TCGA_tpm.Properties.VariableNames{i+1}; % Assign the model name
end

baseModel = prepData.refModel;
% build a matrix indicating which reactions are on the model
compMat = false(length(baseModel.rxns), length(Models_T));
for i = 1:size(compMat,2)
    compMat(:,i) = ismember(baseModel.rxns,Models_T{i}.rxns);
end
% run t-sne
rng(1);  %set random seed to make reproducible
% Without setting the random number generator's seed, the results of the algorithm could potentially vary from run to run, even with the same input data and parameters, due to the inherent randomness involved in some of the algorithm's steps.
proj_coords = tsne(double(compMat'), 'Distance', 'hamming', 'NumDimensions', 2, 'Exaggeration', 6, 'Perplexity', 6);% basically haming distances between models

% export to R
d = struct();
d.tsneX = proj_coords(:, 1);
d.tsneY = proj_coords(:, 2);

%save('TCGA_tsne.mat', 'd');
%save('KPS_TSNE.mat', 'd');
%save('Models_T_P.mat', 'd');

%%
taskFileName = '{PATH}/Human-GEM/data/metabolicTasks/metabolicTasks_Essential.txt';

model_ids = TCGA_tpm.Properties.VariableNames(2:7);
models = Models_T;

model_ids = TCGA_tpm.Properties.VariableNames(2:7); 
model_ids([5]) =[]
models = Models_T;
models([5]) =[] %(Removing Recurrent model and mid_TPS model since they are not included in this study)
%% Structural comparison of the models
% Taken from Tutorial (HUMAN-GEM)
model_ids = arrayfun(@(i) models{i}.id, (1:numel(models))', 'UniformOutput', false);
res = compareMultipleModels(models,false,false,[],false,taskFileName);
clustergram(res.structComp, 'Symmetric', false, 'Colormap', 'bone', 'RowLabels', res.modelIDs, 'ColumnLabels', res.modelIDs);
rxn2Dmap = tsne(res.reactions.matrix', 'Distance', 'hamming', 'NumDimensions', 2, 'Perplexity', 4);

% plot and label the GEMs in tSNE space
scatter(rxn2Dmap(:,1), rxn2Dmap(:,2));
hold on
text(rxn2Dmap(:,1), rxn2Dmap(:,2), res.modelIDs);

%% Subsystem-based Structural comparison reveals altered subsystems
useModels = res.modelIDs;
keep = ismember(res.modelIDs, useModels);
subMat = res.subsystems.matrix(:, keep);

subCoverage = (subMat - mean(subMat, 2)) ./ mean(subMat, 2) * 100;
% select subsystems to include in plot
inclSub = any(abs(subCoverage) > 25, 2); % 25 could be  negative or positive relative differences
subNames = res.subsystems.ID(inclSub);

% generate clustergram
cg = clustergram(subCoverage(inclSub,:), 'Colormap', redbluecmap, 'DisplayRange', 100, 'rowLabels', subNames, 'columnLabels', useModels, 'ShowDendrogram', 'OFF');
cg.DisplayRatio = [0.2,0.1];
cg.RowLabelsRotate = 60;

% cg.RowLabels{2,1} = 'Glycosphingolipid biosynthesis\nt-lacto and neolacto series'
% Increase space on the right side for Y-axis labels by adjusting figure position and size
% hFig = plot(cg);
% pos = get(hFig, 'Position');
% pos(3) = pos(3) + 200;  % Adjust the width as needed
% set(hFig, 'Position', pos);

%% Heatmap
figure('Position', [100, 100, 1200, 800]); % Adjust figure size 
h = heatmap(useModels, subNames, subCoverage(inclSub,:));

% Customizations
h.Colormap = redbluecmap;
h.ColorLimits = [-100, 100]; % Adjust range 
h.GridVisible = 'off'; % Hide grid lines
h.FontSize = 16; % Label rotation and font size adjustments
%%
% The larger metabolic task file "metabolicTasks_Full.txt" can be found in the /data/metabolicTasks/ 
% subdirectory of the Human-GEM repository.
taskFileName = '{PATH}/Human-GEM/data/metabolicTasks/metabolicTasks_Full.txt';
% Re-run the compareMultipleModels function, now including the location of the metabolic task file.
res_func = compareMultipleModels(models(keep), false, false, [], true, taskFileName);

% The ERROR messages here are not actual errors, but indicate that a task was 
% failed because the GEM was missing one or more metabolites involved in the test.

% Identify which tasks differed among the GEMs
isDiff = ~all(res_func.funcComp.matrix == 0, 2) & ~all(res_func.funcComp.matrix == 1, 2); 
% will give a logical vector indicaying zero and non-zero depending on res_func structure. 

% It includes all metabolic task which is checked one by one for each model
diffTasks = res_func.funcComp.tasks(isDiff);

% visualize the matrix
spy(res_func.funcComp.matrix(isDiff,:), 50);
% The spy() function is a built-in MATLAB function that is used to create 
% a plot of the sparsity pattern of a matrix. It creates a plot with black 
% pixels for the non-zero entries in the matrix and white pixels for the zero entries. 
% This function is particularly useful for visualizing sparse matrices.

% apply some formatting changes
set(gca, 'XTick', 1:numel(useModels), 'XTickLabel', useModels, 'XTickLabelRotation', 90, ...
    'YTick', 1:numel(diffTasks), 'YTickLabel', diffTasks, 'YAxisLocation', 'right');
xlabel(gca, '');

% gca is a function in MATLAB that stands for "get current axes." 
% It returns a handle to the current axes in the current figure.
%% Evaluating models for flux capacity to ascertain if they can generate biomass
load('Human-GEM.mat');  % loads model as structure named "ihuman"
ihuman = addBoundaryMets(ihuman);
ihuman.rxns(ihuman.c == 1)
sol = solveLP(ihuman)
Models_T{4,1}.rxns(Models_T{4,1}.c == 1) % ceheck biomass reaction. We used generic biomass fuction.

simplifyModel_T= simplifyModel(Models_T{6,1});
sol = solveLP(simplifyModel_T)

constructEquations(Models_T{4,1}, 'MAR13082')% MAR13082 is a generic biomass reaction.
% ihuman = setParam(ihuman, 'obj', 'MAR13082', 1); setting objective
% function
% Models_T{1,1} = setExchangeBounds(Models_T{1,1}, 'glucose', -1);  % negative flux indicates import
sol = solveLP(Models_T{4,1})
%%
ihuman = simplifyModel(ihuman);
sol = solveLP(ihuman); % FBA
printFluxes(ihuman,sol.x,true,10^-5,[],'%rxnID (%eqn):%flux\n'); %true is only for input and output
sol2 = solveLP(ihuman,1); % pFBA
printFluxes(ihuman,sol.x,true,10^-5,[],'%rxnID (%eqn):%flux\n'); %true is only for input and output


%% REporter Met Analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cd '{PSTH}/RAVEN_sysmedicine'

model_Low_TPS = simplifyModel(Models_T{4,1}); % remove bounday'


% Models_T{6,1}.rxns(Models_T{6,1}.c == 1)
% index = strfind(model_Low_TPS.rxns, 'MAR13082');
% constructEquations(model_Low_TPS, 'MAR13082')
% model_Low_TPS = setParam(model_Low_TPS, 'obj', 'MAR13082', 1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Model Simulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% constructEquations(Models_T{1,1}, 'MAR13082')
sol = solveLP(model_Low_TPS); % FBA
printFluxes(model_Low_TPS,sol.x,true,10^-5,[],'%rxnID (%eqn):%flux\n'); %true is only for input and output
sol2 = solveLP(model_Low_TPS,1); % pFBA
printFluxes(model_Low_TPS,sol.x,true,10^-5,[],'%rxnID (%eqn):%flux\n'); %true is only for input and output

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Run reporter metabolite analysis with DE genes from GBM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cd '{PATH}'; %  DEG_Res direction
% "ensemble_id","NT","TR","TP","low_tps","mid_tps","high_tps"
% 'NT_TR';modelT=simplifyModel(Models_T{2,1}); % remove bounday
% 'NT_TP';modelT=simplifyModel(Models_T{3,1}); % remove bounday
% 'NT_low_tps_count';modelT=simplifyModel(Models_T{4,1}); % remove bounday
% 'NT_high_tps_count';modelT=simplifyModel(Models_T{6,1}); % remove bounday
% 'low_tps_high_tps_count';modelT=simplifyModel(Models_T{6,1}); % remove bounday
% 'low_tps_high_tps_count_2';modelT=simplifyModel(Models_T{4,1}); % remove bounday

% Import the data via DEG_RES_import.m

dif_name='Low_vs_High_C';modelT=simplifyModel(Models_T{4,1}); % remove bounday
%DEG_Res = readtable(strcat('DESeq_Table_dds_',dif_name,'.xlsx'));
% Properties.VariableNames
DEG_Res = DESeq_Table_dds_high_low;
genes1=DEG_Res.rownames;
lgFC1=DEG_Res.log2FoldChange*(-1);
pvalues1=DEG_Res.padj;

% Run reporter metabolite analysis
cd '{PATH}/RAVEN_sysmedicine'
outFileName1 = strcat('reporterMetabolites_',dif_name,'.txt');
tic
repMets1 = reporterMets(modelT,genes1,pvalues1,true,outFileName1,lgFC1);
toc

all_1=repMets1(1,:);
only_up_1=repMets1(2,:);
only_down_1=repMets1(3,:);

colnames= {'mets' 'metNames' 'metZScores' 'metPValues' 'metNGenes' 'meanZ' 'stdZ'};
all_1_T = table(all_1.mets,all_1.metNames,all_1.metZScores,all_1.metPValues,all_1.metNGenes,all_1.meanZ, all_1.stdZ, 'VariableNames', colnames);
only_up_1_T = table(only_up_1.mets,only_up_1.metNames,only_up_1.metZScores,only_up_1.metPValues,only_up_1.metNGenes,only_up_1.meanZ, only_up_1.stdZ, 'VariableNames', colnames);
only_down_1_T = table(only_down_1.mets,only_down_1.metNames,only_down_1.metZScores,only_down_1.metPValues,only_down_1.metNGenes,only_down_1.meanZ, only_down_1.stdZ, 'VariableNames', colnames);

cd '{PATH}'; %  DEG_Res direction
writetable(all_1_T,strcat('all_1_T',dif_name,'.xlsx'));
writetable(only_up_1_T,strcat('only_up_1_T',dif_name,'.xlsx'));
writetable(only_down_1_T,strcat('only_down_1_T',dif_name,'.xlsx'));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Create a cell array of variable names and corresponding data
variables = {'NT', 'TR', 'TP', 'Low_TPS', 'Mid_TPS', 'High_TPS'};
%data = {NT, TR, TP, Low_TPS, Mid_TPS, High_TPS};

% Create a new Excel file
filename = 'Models_T_inside_2.xlsx';

for i = 1:numel(Models_T)
    sheetname = Models_T{i,1}.id;
    data1 = constructEquations(Models_T{i,1}, Models_T{i,1}.rxns);
    data2 = Models_T{i,1}.rxns;
    data3 = Models_T{i,1}.rxnNames;
    data4 = Models_T{i,1}.grRules;
    data5 = Models_T{i,1}.subSystems;
    data6 = Models_T{i,1}.eccodes;
    data7 = Models_T{i,1}.rxnReferences;
    
    % Convert cell column vectors to strings
    data5 = cellfun(@(x) join(string(x), ','), data5, 'UniformOutput', false);

    % Combine data into a single cell array
    data = [data2(:), data3(:), data1(:), data4(:), data5(:), data6(:), data7(:)];
    
    % Convert any numerical data to strings (if needed)
    data = cellfun(@string, data, 'UniformOutput', false);

    % Convert the cell array to a table
    T = cell2table(data, 'VariableNames', {'rxns', 'rxnNames', 'Reactions', 'grRules', 'subSystems', 'eccodes', 'rxnReferences'});

    % Write the table to the Excel file
    writetable(T, filename, 'Sheet', sheetname);
end