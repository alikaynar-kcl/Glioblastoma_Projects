%% Single Gene delation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load('path/Human-GEM/prepData.mat')
load('path/Models_T.mat')
% cd 'path/Lab'

load('Human-GEM.mat');
ihuman = addBoundaryMets(ihuman); % Preparing model for tINIT, closing model

%% if needed unlock chunk: removed selected reaction from ihuman, but grRule became problematic 
for i=1:6;
 model=Models_T{i,1};
 modelname=model.description;
% model = simplifyModel(model); % remove bounday
% model = setExchangeBounds(model, 'glucose', cLB(1),cUB(1), false);  % negative flux indicates import
% sol = solveLP(model); % 
% FluxModel=sol.x;
% modelRXNS= model.rxns;
% Rebuild models (generated using tINIT ) using removeReaction implemented on ihuman
% modelRXNS= [model.rxns;{'MAR00021';'MAR03964'}];
% model.rxns(model.c == 1)
% diffContent = setdiff(ihuman.rxns, modelRXNS);
% reducedModel=removeReactions(ihuman,diffContent,true, true,true);
% reducedModel.rxns(reducedModel.c == 1)
% Example of addReraction
% model = addReaction(model,'rna','reactionFormula','0.18 atp[c] + 0.3 ctp[c] + 0.34 gtp[c] + 0.18 utp[c] -> ppi[c] + rna[c]');
%% Reaction addition is easy with COBRA structure
% ihuman = ravenCobraWrapper(ihuman); %  Converting RAVEN structure to COBRA..
% printRxnFormula(ihuman,'rxnAbbrList',{'MAR00021'});
% printRxnFormula(ihuman,'rxnAbbrList',{'MAR03964'});
% constructEquations(ihuman, 'MAR03964');%'MAR00021')% 
% ihuman = ravenCobraWrapper(ihuman);
% 
% model = ravenCobraWrapper(model); %  Converting RAVEN structure to COBRA..
% %model = addReaction(model,'MAR00021','reactionFormula','MAM01307c[c] + MAM01365c[c] + MAM01369c[c] + MAM01370c[c] + MAM01450c[c] + MAM01451r[r] + MAM01589c[c] + MAM01628c[c] + MAM01721n[n] + MAM01722n[n] + MAM01974c[c] + MAM01975c[c] + MAM01986c[c] + MAM02125c[c] + MAM02184c[c] + MAM02360c[c] + MAM02392c[c] + MAM02426c[c] + MAM02471c[c] + MAM02724c[c] + MAM02733c[c] + MAM02750c[c] + MAM02770c[c] + MAM02847c[c] + MAM02896c[c] + MAM02908c[c] + MAM02993c[c] + MAM03089c[c] + MAM03101c[c] + MAM03135c[c] + MAM03161c[c] + MAM01602c[c]  -> MAM03970c[c]');
% model = addReaction(model,'MAR03964','reactionFormula','MAM01371c[c] + MAM02040c[c]  -> MAM01285c[c] + MAM02039c[c] + MAM02751c[c]'); % reperesent ATP phosphohydrolase, some reaction that need ATP including transport activity
% % printRxnFormula(model,'rxnAbbrList',{'MAR00021'});
% printRxnFormula(model,'rxnAbbrList',{'MAR03964'});
% printRxnFormula(model,'rxnAbbrList',{'MAR13082'});
% % printRxnFormula(model,'rxnAbbrList',{'MAR06916'});
% 
% model = ravenCobraWrapper(model);
% % constructEquations(model, 'MAR00021')
% % constructEquations(model, 'MAR03964')
% % constructEquations(model, 'MAR13082')

essentialTasks = parseTaskList('path/Human-GEM/data/metabolicTasks/metabolicTasks_Essential.txt');
model = simplifyModel(model); % remove bounday
% [outModel, modelChanges]=gapfill4EssentialTasks(model,ihuman,false); % false means don't reset biomass reation
% model=outModel;
model = generateRules(model);
model.rxns(model.c == 1)
constructEquations(model, 'MAR13082') % (MAR13082 biomass)
model= setParam(model,'obj','MAR13082',1); 
% model= setParam(model,'obj','MAR00021',1); Biomass with expand amino
% acid: zoro flux 'MAR00021'
% model= setParam(model,'obj','MAR06916',1);% (MAR06916 Mitocondri) ATP production
% metabolites; -1000 fluxes
sol = solveLP(model)

detectDuplicateRxns(model)% Checking dublicate reaction


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[~, indexes] = ismember(model.rxns, ihuman.rxns);
% model.rxnReferences=ihuman.rxnReferences(indexes);
% emptyCell = cell(length(model.rxns)-length(model.rxnReferences), 1);
% emptyCell(:) = {''};
% model.rxnReferences=vertcat(model.rxnReferences,emptyCell);
% [essentialRxns, essentialRxnsIndexes]=getEssentialRxns(model); Doesn't
% give any results

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Load constraints for RPMI-1640
% Specify the directory and file name
directory = 'path/Lab';
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
% true: close exchange reactions for all other exchanged metabolites not present in the provided list

%modelF = setExchangeBounds(modelF, 'lipoic acid', -1000,1000, false);
sol = solveLP(modelF);
sol = solveLP(modelF,1); % pFBA
FluxModel=horzcat(FluxModel,sol.x);
%modelF=model;
end

[~, indexes] = ismember(vertcat(RPMI_1640_Media(:,1),{'MAR10024'}), model.rxns);
RPMI_1640_FluxEffect=FluxModel(indexes,:);
rowNames= vertcat(RPMI_1640_Media(:,1),{'MAR10024'});
colNames= horzcat({'wConstraint'},RPMI_1640_Media{:,4});
T = array2table(RPMI_1640_FluxEffect, 'RowNames', rowNames, 'VariableNames', colNames);

%T1_seperate_met_constraint=T;
%setHamsMedium(model)
[~, index] = ismember(rowNames, modelF.rxns);
modelF.ub(index)
modelF.lb(index)

%% % checkTasksGenes
%   Performs a set of simulations as defined in a task file. This function
%   is identical to "checkTasks", except it allows for the determination of
%   essential genes, rather than essential reactions.
% Define the directory path

inputFile= 'path/Human-GEM/data/metabolicTasks/metabolicTasks_Essential.txt';
essentialTasks = parseTaskList('path/Human-GEM/data/metabolicTasks/metabolicTasks_Essential.txt');

% 
 modelF = addBoundaryMets(modelF)
[outModel, modelChanges]=gapfill4EssentialTasks(modelF,ihuman,false); % 

tic
[taskReport, essentialGenes, taskStructure]=checkTasksGenes(outModel,[],true,true,true,essentialTasks)
toc
[rowIndices, colIndices] = find(essentialGenes == 1);
%% Cobra Function
modelCobra=outModel;
modelCobra=simplifyModel(modelCobra)
% modelCobra = ravenCobraWrapper(modelCobra);
% [~, indexes] = ismember(modelCobra.rxns, ihuman.rxns);
%  modelCobra.rxnReferences=ihuman.rxnReferences(indexes);

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
directory2 = 'path/';

% Create the folder
splitStrings = split(modelname, ' ');
folderName = strcat(splitStrings{3},'_',splitStrings{5});
folderPath = fullfile(directory2, folderName);
if ~isfolder(folderPath)
    mkdir(folderPath);
end

save(strcat(directory2,folderName,'/taskReport.mat'), 'taskReport');
save(strcat(directory2,folderName,'/essentialGenes.mat'), 'essentialGenes');
save(strcat(directory2,folderName,'/Res_SGD_C.mat'), 'Res_SGD_C');

% Find row and column indices with value 1
[rowIndices, colIndices] = find(essentialGenes == 1);
rowNamess=outModel.genes(rowIndices);
uniqueRowNames = genvarname(rowNamess);
VariableNamess=taskReport.description(colIndices);
uniqueVariableNamess = genvarname(VariableNamess);
TT = array2table(essentialGenes(rowIndices, colIndices), 'RowNames', uniqueRowNames, 'VariableNames', uniqueVariableNamess);

% Specify the Excel file path
excelFilePath = fullfile(folderPath, 'RPMI_1640_FluxEffect.xlsx');

% % Save the variables as an Excel file
% writeDataToExcel(excelFilePath, RPMI_1640_FluxEffect);

tableFilePath = fullfile(folderPath, 'RPMI_1640_FluxEffect.xlsx');
tableFilePath1 = fullfile(folderPath, 'essentialGenes.xlsx');

% Save the table as an Excel file
writetable(T, tableFilePath,'WriteRowNames', true);
writetable(TT, tableFilePath1,'WriteRowNames', true);
end
%%
filename2 = 'output_SGD.xlsx';  % Replace with the desired Excel file name
for i=1:6
 model=Models_T{i,1};
 modelname=model.description;
splitStrings = split(modelname, ' ');
folderName = strcat(splitStrings{3},'_',splitStrings{5});

% Specify the directory and file name
directort0='/Users/alikaynar/Desktop/UpgrdPostDoc/GEM/GEM_R_Second_Project';
directory3 = strcat(directort0,'/',folderName);  % Replace with the actual directory path
filename = 'Res_SGD_C.mat';  % Replace with the actual file name

% Construct the full file path
fullFilePath = fullfile(directory3, filename);

% Load the .mat file
load(fullFilePath);
modelCobra=model;
%[rowIndices, colIndices] = find(Res_SGD_C.grRateKO == 0);
[rowIndices, colIndices] = find(Res_SGD_C.grRatio < 0.2);
EssentialGenes_SGD_C=modelCobra.genes(rowIndices);
[i,length(rowIndices)]
% Write variable1 to Sheet1
outputSheet1 = fullfile(directort0, filename2);
writecell(EssentialGenes_SGD_C, outputSheet1, 'Sheet', folderName)
end


