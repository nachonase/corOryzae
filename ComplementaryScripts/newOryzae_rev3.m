%% Reconstruction process
% Purposed: Here is a complete script used for reconstructing a genome-scale metabolic model of 
% a cordycepin-producing Aspergillus oryzae strain.  
% 
% Written by Nachon Raethong, 30-AUG-2023
% Updated by Nachon Raethong, 05-Jan-2024
% Updated by Nachon Raethong, 01-Feb-2024

%% WORKSPACE
cd '/Users/nachonase/Documents/GitHub/BCC7051-GSMM';
%% STEP 1: Get model from Cordyceps militaris (iNR1329)
% Load model template, iNR1329
load 'ComplementaryData/iNR1329.mat';
cordModel = model;

% Identification of orthologous proteins between A. oryzae BCC7051 and iNR1329
blastedBCC7051Cordyceps=getBlast('BCC7051','ComplementaryData/BCC7051.fasta',...
    'iNR1329','ComplementaryData/Cmilitaris_protein.fasta');
%save('ComplementaryData/blastedBCC7051Cordyceps.mat','blastedBCC7051Cordyceps');
%load 'ComplementaryData/blastedBCC7051Cordyceps.mat';
 
% Creation of draft networks from protein orthology inference
bccDraftFromCordyceps = getModelFromHomology({cordModel},...
    blastedBCC7051Cordyceps,'BCC7051',{},1,false,10^-50,200,40);
bccDraftFromCordyceps.id = 'bccDraftFromCordyceps';

% Add Non-gene-associated reactions, including spontaneous, 
% transport and exchange reactions from iNR1329

rxnToAdd.rxns =setdiff(cordModel.rxns,bccDraftFromCordyceps.rxns);
bccDraftFromCordycepsPlusiNR1329 = addRxnsGenesMets(bccDraftFromCordyceps,...
    cordModel,rxnToAdd.rxns,false,...
    'additional rxns from iNR1329 model',2);

% Check flux to Biomass
finalValidateModel = bccDraftFromCordycepsPlusiNR1329;
finalValidateModel = setParam(finalValidateModel,'eq',{'matp'},1);
finalValidateModel = setParam(finalValidateModel,'obj',{'bmOUT'},1);
model=setParam(finalValidateModel,'ub',{'glcIN'},0.1448);
sol = solveLP(model,1);
fprintf(['umax = ' num2str(sol.f*-1) ' per hour' '\n']);
run = finalValidateModel;
run.description = 'BCC7051 from iNR1329';
run.id = 'modelBCC7051fromiNR1329';

% Note: modelBCC7051fromiNR1329 uses biomass reactions of Cordyceps and includes
% many gap reactions from iNR1329.

%% STEP 2: Improve modelBCC7051fromiNR1329 by iWV1346 
% Load iWV1346 model
aoryModel = importModel('ComplementaryData/iWV1346.xml');
% To avoid dupplicated metabolites, we did metabolite naming.
reduced = aoryModel;
[~, textData1]=xlsread('ComplementaryData/supplementary.xlsx','Table S2');
metNames = struct();
metNames.old = textData1(4:end,2);
metNames.new = textData1(4:end,3);
[a, b]=ismember(reduced.metNames,metNames.old);
I=find(a);
reduced.metNames(I)=metNames.new(b(I));
reduced=contractModel(reduced);
% Check flux to A. oryzae biomass
reduced = setParam(reduced,'ub',{'r2205'},[1]);
reduced = setParam(reduced,'obj',{'r2359'},1);
sol = solveLP(reduced);
printFluxes(reduced, sol.x, true);

% Identification of orthologous proteins A. oryzae BCC7051 and iWV1346
blastedBCC7051Aoryzae=getBlast('BCC7051','ComplementaryData/BCC7051.fasta',...
    'iWV1346','ComplementaryData/A_oryzae.faa');
%save('ComplementaryData/blastedBCC7051Aoryzae.mat','blastedBCC7051Aoryzae');
%load 'ComplementaryData/blastedBCC7051Aoryzae.mat';

% Creation of draft networks from protein orthology inference
bccDraftFromAoryzae = getModelFromHomology({reduced},...
    blastedBCC7051Aoryzae,'BCC7051',{},1,false,10^-50,200,40);
bccDraftFromAoryzae.id = 'bccDraftFromAoryzae';
bccDraftFromAoryzae.description = 'BCC7051 from iWV1346';
%% Merge run and Aor
model=mergeModels({run,bccDraftFromAoryzae});
model.description = 'BCC7051 from mergeModels';
reduced2=expandModel(model);
reduced3=contractModel(reduced2);
reduced3.grRules=standardizeGrRules(reduced3);
reduced5 = setParam(reduced3,'obj','bmOUT',1);
reduced5 = setParam(reduced5,'eq','matp',1);
requiredRxns = {'o2IN' 'piIN' 'slfIN' 'nh3IN'};
reduced5 = setParam(reduced5,'ub',requiredRxns,...
    [1000,1000,1000,1000]);
reduced5 = setParam(reduced5,'ub','glcIN',0.2965);
sol = solveLP(reduced5);
printFluxes(reduced5, sol.x);
reduced5.description = 'BCC7051 from mergeModels but not BCC7051 biomass';
reduced5.id = 'reduced5';


%% Import bcc7051 biomass
BCC7051biomass = importExcelModel('ComplementaryData/BCC7051_biomass3.xlsx');
bcc7051BiomassModel = addRxnsGenesMets(reduced5,...
    BCC7051biomass,BCC7051biomass.rxns,false,...
    'import bcc7051 biomass',2);
bcc7051BiomassModel.description = 'model with BCC7051 biomass';
bcc7051BiomassModel.id = 'bcc7051BiomassModel';
sol = solveLP(bcc7051BiomassModel);
printFluxes(bcc7051BiomassModel, sol.x);
%% Remove duplicates and bad reactions
badrxnsToRemove = {'r0069', 'r0135', 'r0136', 'r0140', 'r0141', 'r0206', 'r0207',...
    'r0209', 'r0215', 'r0424', 'r0563', 'r0572', 'r0607', 'r0625', 'r0638', 'r0639',...
    'r0693', 'r0705', 'r0708', 'r0722', 'r0726', 'r0759', 'r1125', 'r1150', 'r0143',...
    'r1124', 'r975', 'r1588', 'R00078_c', 'R00132_c', 'R00317_c', 'R00533_c', 'R00565_c',...
    'R01352_c', 'R01416_c', 'R01451_c', 'R01645_c', 'R01663_c', 'R01845_c', 'R02025_c',...
    'R02208_c', 'R02422_c', 'R02556_c', 'R02964_c', 'R03308_c', 'R03409_c', 'R03435_c',...
    'R03478_c', 'R04452_c', 'R04734_c', 'R05202_c', 'R06264_c', 'R06513_c', 'R06518_c',...
    'R06601_c', 'R07495_c', 'R07766_c', 'R08266_c', 'R09030_c', 'R09087_c', 'R09106_c',...
    'R09395_c', 'R09597_c', 'R10231_c', 'R10309_c', 'R10722_c', 'R10815_c', 'R11218_c',...
    'R11308_c', 'R11861_c', 'r0113', 'r0117', 'r0119', 'r0120', 'r0122', 'r0152',...
    'r0205', 'r0219', 'r0223', 'r0297', 'r0298', 'r0312', 'r0408', 'r0441', 'r0442',...
    'r0467', 'r0477', 'r0492', 'r0505', 'r0518', 'r0528', 'r0543', 'r0609', 'r0610', 'r0613',...
    'r0716', 'r0753', 'r1122', 'r1154', 'r1162', 'r1526', 'wv_R00014', 'wv_R00066',...
    'wv_R00079', 'wv_R00110', 'wv_R00230', 'wv_R00238', 'wv_R00241', 'wv_R00261', ...
    'wv_R00269', 'wv_R00320', 'wv_R00367', 'wv_R00373', 'wv_R00497', 'wv_R00520',...
    'wv_R00523', 'wv_R00524', 'wv_R00577', 'wv_R00645', 'wv_R00653', 'wv_R00655',...
    'wv_R00656', 'wv_R00691', 'wv_R00788', 'wv_R00791', 'wv_R00796', 'wv_R01102', 'wv_R01129',...
    'wv_R01130', 'wv_R01139', 'wv_R01174', 'wv_R01177', 'wv_R01184', 'wv_R01222', 'wv_R01246',...
    'wv_R01258', 'wv_R01264', 'r103', 'r158', 'r159', 'r191', 'r236', 'r274', 'r351', 'r353', ...
    'r458', 'r464', 'r502_bccDraftFromAoryzae', 'r507_bccDraftFromAoryzae',...
    'r508_bccDraftFromAoryzae', 'r509_bccDraftFromAoryzae', 'r559', 'r564', 'r570', ...
    'r571', 'r592', 'r623', 'r685', 'r686', 'r687', 'r728', 'r749', 'r754', 'r769', ...
    'r771', 'r773', 'r798', 'r799', 'r800', 'r811', 'r816', 'r833', 'r835', 'r855', ...
    'r868', 'r878', 'r888', 'r897', 'r899', 'r900', 'r930', 'r939', 'r964', 'r966', ...
    'r969', 'r971', 'r1030_bccDraftFromAoryzae', 'r1033_bccDraftFromAoryzae', ...
    'r1035_bccDraftFromAoryzae', 'r1040_bccDraftFromAoryzae', 'r1046_bccDraftFromAoryzae',...
    'r1049_bccDraftFromAoryzae', 'r1055_bccDraftFromAoryzae', 'r1071_bccDraftFromAoryzae',...
    'r1072_bccDraftFromAoryzae', 'r1073_bccDraftFromAoryzae', 'r1079_bccDraftFromAoryzae', ...
    'r1109_bccDraftFromAoryzae', 'r1110_bccDraftFromAoryzae', 'r1127_bccDraftFromAoryzae', ...
    'r1154_bccDraftFromAoryzae', 'r1155_bccDraftFromAoryzae', 'r1156', 'r1157',...
    'r1158_bccDraftFromAoryzae', 'r1159', 'r1160_bccDraftFromAoryzae', 'r0127',...
    'r1161_bccDraftFromAoryzae', 'r1162_bccDraftFromAoryzae', 'r1163', ...
    'r1164_bccDraftFromAoryzae', 'r1165_bccDraftFromAoryzae', 'r1166_bccDraftFromAoryzae',...
    'r1167_bccDraftFromAoryzae', 'r1168_bccDraftFromAoryzae', 'r1169_bccDraftFromAoryzae', ...
    'r1170_bccDraftFromAoryzae', 'r1171_bccDraftFromAoryzae', 'r1172_bccDraftFromAoryzae', ...
    'r1288_bccDraftFromAoryzae', 'r1293_bccDraftFromAoryzae', 'r1298_bccDraftFromAoryzae', ...
    'r1303_bccDraftFromAoryzae', 'r1308_bccDraftFromAoryzae', 'r1313', ...
    'r1318_bccDraftFromAoryzae', 'r1323_bccDraftFromAoryzae', 'r1328', ...
    'r1333_bccDraftFromAoryzae', 'r1338_bccDraftFromAoryzae', 'r1343_bccDraftFromAoryzae', ...
    'r1348_bccDraftFromAoryzae', 'r1357_bccDraftFromAoryzae', 'r1358_bccDraftFromAoryzae', ...
    'r1363_bccDraftFromAoryzae', 'r1364_bccDraftFromAoryzae', 'r1369_bccDraftFromAoryzae', ...
    'r1370_bccDraftFromAoryzae', 'r1375_bccDraftFromAoryzae', 'r1376_bccDraftFromAoryzae', ...
    'r1381_bccDraftFromAoryzae', 'r1382_bccDraftFromAoryzae', 'r1387_bccDraftFromAoryzae', ...
    'r1388_bccDraftFromAoryzae', 'r1392_bccDraftFromAoryzae', 'r1393_bccDraftFromAoryzae', ...
    'r1417_bccDraftFromAoryzae', 'r1429_bccDraftFromAoryzae', 'r1438_bccDraftFromAoryzae', ...
    'r1443_bccDraftFromAoryzae', 'r1451', 'r1461', 'r1466', 'r1528_bccDraftFromAoryzae', ...
    'r1533', 'r1535', 'r1583', 'r1589_bccDraftFromAoryzae', 'r1685', 'r1686', 'r1689', ...
    'r1734', 'r1737', 'r1739', 'r1743', 'r1749', 'r1780', 'r1798', 'wv_R00572', ...
    'wv_R00611'};
    
goodModel = removeReactions(bcc7051BiomassModel,badrxnsToRemove,true,...
            false,true);

sol = solveLP(goodModel);
printFluxes(goodModel, sol.x);

goodModel2 = setParam(goodModel,'eq','biomass',0);
goodModel2 = setParam(goodModel,'eq','lipid',0);
goodModel2 = setParam(goodModel2,'eq','rna',0);

sol = solveLP(goodModel2);
printFluxes(goodModel2, sol.x);



cordycepsbiomassrxnsToRemove = {'biomass' 'carbohydrate' 'dna' 'lipid' 'matp' 'protein' 'rna'}; 
nocordycepsbcc7051BiomassModel = removeReactions(goodModel2,cordycepsbiomassrxnsToRemove,true,...
            false,true);
nocordycepsbcc7051BiomassModel.id = 'nocordycepsbcc7051BiomassModel';
%adjusting rev of some reactions
adjModel = nocordycepsbcc7051BiomassModel;
adjModel.rev(31) = 1;
adjModel.rev(1743) = 1;
adjModel.rev(1674) = 1;
adjModel.rev(1675) = 1;
adjModel=setParam(adjModel,'lb',{'r0041'},-1000);
adjModel=setParam(adjModel,'lb',{'r666'},-1000);
adjModel=setParam(adjModel,'lb',{'r12'},-1000);
adjModel=setParam(adjModel,'lb',{'r13'},-1000);
adjModel.equations = constructEquations(adjModel);
adjModel.id = 'adjModel';

goodModel22=setParam(adjModel,'ub',{'glcIN'},0.2965);
%goodModel22=setParam(adjModel,'ub',{'glcIN'},0);

sol = solveLP(goodModel22);
sol.f

printFluxes(goodModel22, sol.x);

%% Improve grRules 
modelToImprove = nocordycepsbcc7051BiomassModel;
expression = '^or ';
replace = '';
for i = 1:numel(modelToImprove.rxns)
    modelToImprove.grRules(i) = regexprep(modelToImprove.grRules(i),expression,replace);
end

%% Improve eccodes 
modelToImproveeccodes = modelToImprove;
for i = 1:numel(modelToImproveeccodes.rxns)
    modelToImproveeccodes.eccodes(i) = regexprep(modelToImproveeccodes.eccodes(i),...
        ';ec-code/',';');
    modelToImproveeccodes.eccodes(i) = regexprep(modelToImproveeccodes.eccodes(i),...
        '^ec-code/','');
    modelToImproveeccodes.eccodes(i) = regexprep(modelToImproveeccodes.eccodes(i),...
        ':',';');
    modelToImproveeccodes.eccodes(i) = regexprep(modelToImproveeccodes.eccodes(i),...
        '^;','');
end

modelToImproveeccodes.id = 'modelToImproveeccodesValidetedgrowth';
%%
%exportForGit(modelToImproveeccodes,'modelToImproveeccodesValidetedgrowth',...
%    '/Users/nachonase/Documents/GitHub/BCC7051-GSMM/model_rev3', {'xml','xlsx'});

%% Gapfilling

% Polishing the metabolic network coverage of A. oryzae BCC7051
% In this case, iNR1329 is the only template model that has cordycepin pathway, 
% so the draft initiated from Cordyceps militaris model (iNR1329) is 
% promisingly used as a main network for A. oryzae BCC7051.
% In particular other drafts created from A. nidulans (ANID), Penicillium (iAL1006)
% and N. crassa (MODEL1212060001) as well as an earlier genome-scale metabolic model 
% of  A. oryzae RIB40 (iWV1346) were used as alternative networks 
% supported for enhancing metabolic coverage of the main network. 
% For instance, by checking the EC numbers, some new reactions were added from 
% iAL1006 with annotated genes from alternative networks.
% curatedModel=importExcelModel('model/xlsx/curated_bccDraftFromCordycepsPlusiNR1329.xlsx');

%based protein orthology from bccDraftFromNcrassa
%r0128 EC:3.6.3.14     OAory_01022650      
%r0175 EC:1.1.1.138	  OAory_01028300 or OAory_01053340 or OAory_01055840 or OAory_01075090 or OAory_01075690 or OAory_01017830 or OAory_01017860 or OAory_01045670 or OAory_01046790 or OAory_01052770 or OAory_01075370 or OAory_01097490

%based protein orthology from bccDraftFromAnidulans
%r0294 EC:2.7.8.8      OAory_01015910

%based on protein orthology from bccDraftFromMetaCycDatabase
%r0448 EC:3.6.1.23 	OAory_01007340
%r0397 EC:2.7.4.6      OAory_01082240
%r0398 EC:2.7.4.6      OAory_01082240
%r0399 EC:2.7.4.6      OAory_01082240
%r0400 EC:2.7.4.6      OAory_01082240
%r0443 EC:2.7.4.6      OAory_01082240
%r0444 EC:2.7.4.6      OAory_01082240
%r0445 EC:2.7.4.6      OAory_01082240
%r0446 EC:2.7.4.6      OAory_01082240
%r0447 EC:2.7.4.6      OAory_01082240

listgenesToAdd = {'OAory_01028300', 'OAory_01007340', ...
    'OAory_01082240', 'OAory_01053340', 'OAory_01055840',...
    'OAory_01075090', 'OAory_01075690'};
genesToAdd = struct();
genesToAdd.genes = cell(numel(listgenesToAdd),1);
for i = 1:numel(genesToAdd.genes)
    genesToAdd.genes{i} = listgenesToAdd{i}
end
newModel = addGenesRaven(modelToImproveeccodes,genesToAdd);

listpolishRxns = {'r0128','r0175', 'r0294','r0448','r0397','r0398','r0399','r0400',...
    'r0443','r0444','r0445','r0446','r0447'};
listpolishgrRules = {'OAory_01022650','OAory_01028300 or OAory_01053340 or OAory_01055840 or OAory_01075090 or OAory_01075690 or OAory_01017830 or OAory_01017860 or OAory_01045670 or OAory_01046790 or OAory_01052770 or OAory_01075370 or OAory_01097490',...
    'OAory_01015910','OAory_01007340','OAory_01082240','OAory_01082240','OAory_01082240',...
    'OAory_01082240','OAory_01082240','OAory_01082240','OAory_01082240','OAory_01082240',...
    'OAory_01082240'};
polish = struct();
polish.rxns = cell(numel(listpolishRxns),1);
polish.grRules = cell(numel(listpolishRxns),1);
for i = 1:numel(polish.rxns)
    polish.rxns{i} = listpolishRxns{i}
    polish.grRules{i} = listpolishgrRules{i}
end

[a, b]=ismember(newModel.rxns,polish.rxns);
I=find(a);
newModel.grRules(I)=polish.grRules(b(I));

reducedModel=deleteUnusedGenes(newModel);
reducedModel.id = 'newOryzae';
exportForGit(reducedModel,'newOryzae',...
    '/Users/nachonase/Documents/GitHub/BCC7051-GSMM/newOryzae', {'xml','xlsx','mat','txt'});



exportForGit(reducedModel,'newModelPolisheddeleteUnusedGenes',...
    '/Users/nachonase/Documents/GitHub/BCC7051-GSMM/model_rev3', {'xml','xlsx'});

%Check transport rxns
coolModel = importModel('model_rev3/model/xml/newModelPolisheddeleteUnusedGenes.xml');
coolModel.equations = constructEquations(coolModel);
transportRxnsnum = getTransportRxns(coolModel);
tranS = struct();
tranS.rxns = coolModel.rxns(transportRxnsnum);
tranS.equations = coolModel.equations(transportRxnsnum);
tranS.grRules = coolModel.grRules(transportRxnsnum);
tranS.sols = cell(numel(tranS.rxns),1);

for i = 1:numel(tranS.rxns)
    tran = tranS.rxns{i};
    model=setParam(coolModel,'ub',tran,0);
    model=setParam(coolModel,'lb',tran,0);
    sol = solveLP(model,1);
    tranS.sols{i} = sol.f(1);
end

