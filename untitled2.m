clc
initCobraToolbox;
gtex_data = readtable('/ezgi/Documents/github/Human1_Publication_Data_scripts/tINIT_GEMs/data/gtex_median_tissue_tpm.txt');
gtex_data(1:5,1:5);
% extract the tissue and gene names
data_struct.tissues = gtex_data.Properties.VariableNames(2:end)';  % sample (tissue) names
data_struct.genes = gtex_data.genes;  % gene names
data_struct.levels = table2array(gtex_data(:, 2:end));  % gene TPM values
data_struct.threshold = 1;

load('Human-GEM.mat');
ihuman = addBoundaryMets(ihuman);
model = ravenCobraWrapper(ihuman);
essentialTasks = parseTaskList('/Users/ezgi/Documents/github/Human-GEM/data/metabolicTasks/metabolicTasks_Essential.xlsx');
checkTasks(ihuman, [], true, false, false, essentialTasks);

refModel = ihuman;  % the reference model from which the GEM will be extracted
tissue = 'brain';  % must match the tissue name in data_struct.tissues
celltype = [];  % used if tissues are subdivided into cell type, which is not the case here
hpaData = [];  % data structure containing protein abundance information (not used here)
arrayData = data_struct;  % data structure with gene (RNA) abundance information
metabolomicsData = [];  % list of metabolite names if metabolomic data is available
removeGenes = true;  % (default) remove lowly/non-expressed genes from the extracted GEM
taskFile = [];  % we already loaded the task file, so this input is not required
useScoresForTasks = true;  % (default) use expression data to decide which reactions to keep
printReport = true;  % (default) print status/completion report to screen
taskStructure = essentialTasks;  % metabolic task structure (used instead "taskFile")
params = [];  % additional optimization parameters for the INIT algorithm
paramsFT = [];  % additional optimization parameters for the fit-tasks algorithm

% brainGEM = getINITModel2(refModel, tissue, celltype, hpaData, arrayData, metabolomicsData, removeGenes, taskFile, useScoresForTasks, printReport, taskStructure, params, paramsFT);
% brainGEM.id = 'brain';

% pred_Grates    = [];
% pred_EC_Grates = [];
% meanError      = [];
% mean_EC_error  = [];

load('brainGEM.mat');
model = brainGEM;
model = ravenCobraWrapper(model);
exchModel = setHamsMedium(model); % adapt the function
cd ..
cd GECKO
cd geckomat
cd utilities
brain_ecModel = getSubset_ecModel(model,refModel);
% cd ..
% cd ..
% cd ecModels
% cd ecHumanGEM
% cd model
% load('ecHumanGEM_batch.mat');
% load('ecHumanGEM.mat');
% refModelEC = ecModel;
% ecBrainGEM = getINITModel2(refModelEC, tissue, celltype, hpaData, arrayData, metabolomicsData, removeGenes, taskFile, useScoresForTasks, printReport, taskStructure, params, paramsFT);
clc
cd ..
cd ..
cd ..
cd brainGEM
initCobraToolbox;
model = changeObjective(model,'biomass_human');
cobra_fba = optimizeCbModel(model,'max','one', (1)); % run FBA
FBAvectors_brain_wloop = cobra_fba.x;
save('FBAvectors_brain_wloop');
cobra_fba = optimizeCbModel(model,'max','one', 0); % run FBA
FBAvectors_brain_woloop = cobra_fba.x;
save('FBAvectors_brain_woloop');
cobra_fba_ec = optimizeCbModel(brain_ecModel,'max','one', (1)); % run FBA
FBAvectors_brain_ec_wloop = cobra_fba_ec.x;
save('FBAvectors_brain_ec_wloop');
cobra_fba_ec = optimizeCbModel(brain_ecModel,'max','one', 0); % run FBA
FBAvectors_brain_ec_woloop = cobra_fba_ec.x;
save('FBAvectors_brain_ec_woloop');

