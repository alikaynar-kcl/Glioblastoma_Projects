%% Single Gene delation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load('/Users/alikaynar/Desktop/UpgrdPostDoc/GEM/GBM_GEM/Human-GEM/prepData.mat')
load('/Users/alikaynar/Desktop/UpgrdPostDoc/GEM/GEMs_R/Models_T.mat')
% cd '/Users/alikaynar/Desktop/UpgrdPostDoc/GEM/GBM_GEM/Lab'

%% if needed unlock chunk: removed selected reaction from ihuman, but grRule became problematic 
for i=1:6;
 model=Models_T{i,1};
 modelname=model.description;
model = simplifyModel(model); % remove bounday
model = generateRules(model);
model.rxns(model.c == 1)
constructEquations(model, 'MAR13082') % (MAR13082 biomass)
model= setParam(model,'obj','MAR13082',1); 

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[~, indexes] = ismember(model.rxns, ihuman.rxns);
 model.rxnReferences=ihuman.rxnReferences(indexes);
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Load constraints for RPMI-1640
% Specify the directory and file name
directory = '/Users/alikaynar/Desktop/UpgrdPostDoc/GEM/GBM_GEM/Lab';
filename = 'RPMI-1640 Media Formulation Model Input 2.xlsx';
filepath = fullfile(directory, filename);% Construct the full file path
[NUM,STR] = xlsread(filepath);% Read the data from the XLSX file
%[NUM,STR] = xlsread('RPMI-1640 Media Formulation Model Input 2.xlsx');
cRxns = STR(2:end,1);
MetS = STR(2:end,2);
cLB = NUM(:,1);
cUB = NUM(:,2);

index = find(model.metComps == 1);
modelMetName = model.metNames(index);
modelMetID = model.mets(index);

% Construct the empty cell array
RPMI_1640_Media = cell(length(NUM), 6);
RPMI_1640_Media(:, 1:4) = {''}; % Set the data type for the first 3 columns as strings
RPMI_1640_Media(:, 5:6) = {[]}; % Set the data type for the last 2 columns as numbers
disp(RPMI_1640_Media); % Display the cell array

sol = solveLP(model); % 
sol = solveLP(model,1); % 
FluxModel=sol.x;
modelF=model;
%model = setExchangeBounds(model, 'lipoic acid', -100,100, false);% lipoic acid prevents biomass.
for i = 1:length(MetS)
str = MetS(i);
splitStr = split(str, '[e]');
index2 = find(strcmp(modelMetName, splitStr(1)));
RPMI_1640_Media{i,1}= STR{i+1,1} ;
RPMI_1640_Media{i,2}= constructEquations(model, STR{i+1,1}) ;
RPMI_1640_Media{i,3}= modelMetID(index2) ;
RPMI_1640_Media{i,4}= splitStr(1) ;
RPMI_1640_Media{i,5}= NUM(i,1) ;
RPMI_1640_Media{i,6}= NUM(i,2) ;
modelF = setExchangeBounds(modelF, modelMetName{index2}, cLB(i),cUB(i), false);
%modelF = setExchangeBounds(modelF, 'lipoic acid', -1,1, false);
sol = solveLP(modelF);
sol = solveLP(modelF,1); % pFBA
FluxModel=horzcat(FluxModel,sol.x);
%modelF=model;
end

 modelF = addBoundaryMets(modelF)
% modelF = simplifyModel(modelF); % remove bounday
[outModel, modelChanges]=gapfill4EssentialTasks(modelF,ihuman,false); % it doesn't work
% [outModel, modelChanges]=gapfill4EssentialTasks(model,ihuman,false); % it doesn't work

%% Cobra Function
modelCobra=outModel;
modelCobra=simplifyModel(modelCobra)
tic
 [grRatio, grRateKO, grRateWT, hasEffect, delRxns, fluxSolution] = singleGeneDeletion(modelCobra);
toc
Res_SGD_C.grRatio=grRatio;
Res_SGD_C.grRateKO=grRateKO;
Res_SGD_C.grRateWT=grRateWT;
Res_SGD_C.hasEffect=hasEffect;
Res_SGD_C.delRxns=delRxns;
Res_SGD_C.fluxSolution=fluxSolution;
%% Save results to each model folder
% Specify the directory where the folder should be created
directory2 = '/Users/alikaynar/Desktop/UpgrdPostDoc/GEM/GEM_R_Second_Project/';

% Create the folder
splitStrings = split(modelname, ' ');
folderName = strcat(splitStrings{3},'_',splitStrings{5});
folderPath = fullfile(directory2, folderName);
if ~isfolder(folderPath)
    mkdir(folderPath);
end

save(strcat(directory2,folderName,'/Res_SGD_C.mat'), 'Res_SGD_C');

end


%%

% Initialize an empty cell array to store the genes
allGenes = {};

% Concatenate genes from Models_T{i,1}.genes
for i = [1, 3, 4, 6]
    allGenes = [allGenes; Models_T{i,1}.genes];
end

% Make the genes unique
MetabolicGenes = unique(allGenes);

dlmwrite('MetabolicGenes.txt', MetabolicGenes, 'delimiter', '\t');

% Save as a text file
fileID = fopen('MetabolicGenes.txt', 'w');
for i = 1:numel(MetabolicGenes)
    fprintf(fileID, '%s\n', MetabolicGenes{i});
end
fclose(fileID);




