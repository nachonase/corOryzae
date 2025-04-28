%% Reconstruction process
% Purposed: Here is a complete script used for reconstructing a metabolic model of 
% Aspergillus oryzae BCC7051 for optimizing cordycepin production.  
% The reconstruction process was carried through the following 4 steps. 
% This is an iterative process which ends up with a validated GEM 
% that can represent A. oryzae BCC7051 metabolism.
% 
% Written by Nachon Raethong, 30-AUG-2023

%% WORKSPACE
cd '/Users/nachonase/Documents/GitHub/BCC7051-GSMM';
%% STEP 1: DATA PREPARATION
% 1.1 Data collection
% Whole-proteome sequences of A. oryzae BCC7051 and template species as well as
% the model files for Cordyceps militaris (iNR1329), A. oryzae RIB40 (iWV1346), 
% Penicillium chrysogenum Wisconsin 54-1255 (iAL1006), 
% Aspergillus nidulans FGSC A4 (ANID) and Neurospora crassa (MODEL1212060001) 
% were collected and keep in /ComplementaryData/ folder. 
% Notably, gene identifiers of these proteome sequences must 
% be same with gene identifiers used in template models.

load 'ComplementaryData/iNR1329.mat';
cordModel = model;
aoryModel = importModel('ComplementaryData/iWV1346.xml');
peniModel = importModel('ComplementaryData/iAL1006 v1.00.xml');
anidModel = importModel('ComplementaryData/ANIDULANS.xml');
neurModel = importModel('ComplementaryData/iJDZ808.xml');

% 1.2 Model templates prepraration and metabolite naming
% The metabolite names used by different template models may introduce a lot of 
% duplicated metabolites and reactions during the reconstruction process. 
% So, metabolites names used in these template
% models were reorganized by mapping metabolite identifiers obtained from 
% different databases (e.g. MetaCyc, KEGG or CHEBI) to the same 
% metabolite names thereby resulted in compatible metabolite namespaces. 
% Overall updated metabolite names and their identifiers are listed in Table S2.
% For addition of new metabolites to the templates network, new metabolites 
% listed in Table S1 were introduced to the templates network by addMets function.

[~, newMets]=xlsread('ComplementaryData/supplementary.xlsx','Table S1');
metsToAdd = struct();
metsToAdd.mets = newMets(2:end,1);
metsToAdd.metNames = newMets(2:end,2);
metsToAdd.metFormulas = newMets(2:end,3);
metsToAdd.compartments = 'c';
peniModel_newMets=addMets(peniModel,metsToAdd); 
[~, textData1]=xlsread('ComplementaryData/supplementary.xlsx','Table S2');
metNames = struct();
metNames.old = textData1(4:end,2);
metNames.new = textData1(4:end,3);
[a, b]=ismember(peniModel_newMets.metNames,metNames.old);
I=find(a);
peniModel_newMets.metNames(I)=metNames.new(b(I));
reduced=contractModel(peniModel_newMets);
reduced = setParam(reduced,'ub',{'glcIN' 'etohIN'},[1 0]);
reduced = setParam(reduced,'obj',{'bmOUT'},1);
sol = solveLP(reduced);
printFluxes(reduced, sol.x, true);
idex_cordycepinPermease = find(ismember(reduced.rxns,'r1468'));
reduced.rxns(idex_cordycepinPermease) = {'nr0010'};
reduced.rxnNames(idex_cordycepinPermease) = {'cordycepin permease'};
idex_cordycepinProduction = find(ismember(reduced.rxns,'penartOUT'));
reduced.rxns(idex_cordycepinProduction) = {'cordycepinOUT'};
reduced.rxnNames(idex_cordycepinProduction) = {'production of cordycepin'};
 

%% STEP 2: DRAFT RECONSTRUCTION
% 2.1 Identification of orthologous proteins 
% A set of orthologous proteins between A. oryzae BCC7051 and template model 
% proteins was obtained by sequence-alignment analysis using BLASTp.
%blastedBCC7051Cordyceps=getBlast('BCC7051','ComplementaryData/BCC7051.fasta',...
%    'iNR1329','ComplementaryData/Cmilitaris_protein.fasta');
%save('ComplementaryData/blastedBCC7051Cordyceps.mat','blastedBCC7051Cordyceps');

%blastedBCC7051Aoryzae=getBlast('BCC7051','ComplementaryData/BCC7051.fasta',...
%    'iWV1346','ComplementaryData/A_oryzae.faa');
%save('ComplementaryData/blastedBCC7051Aoryzae.mat','blastedBCC7051Aoryzae');

%blastedBCC7051Anidulans=getBlast('BCC7051','ComplementaryData/BCC7051.fasta',...
%    'ANID','ComplementaryData/Anidulans.fasta');
%save('ComplementaryData/blastedBCC7051Anidulans.mat','blastedBCC7051Anidulans');

%blastedBCC7051Ncrassa=getBlast('BCC7051','ComplementaryData/BCC7051.fasta',...
%    'MODEL1212060001','ComplementaryData/Ncrassa.fasta');
%save('ComplementaryData/blastedBCC7051Ncrassa.mat','blastedBCC7051Ncrassa');

%blastedBCC7051Penicillium=getBlast('BCC7051','ComplementaryData/BCC7051.fasta',...
%    'iAL1006','ComplementaryData/Penicillium.fasta');
%save('ComplementaryData/blastedBCC7051Penicillium.mat','blastedBCC7051Penicillium');

load 'ComplementaryData/blastedBCC7051Cordyceps.mat';
load 'ComplementaryData/blastedBCC7051Aoryzae.mat';
load 'ComplementaryData/blastedBCC7051Anidulans.mat';
load 'ComplementaryData/blastedBCC7051Penicillium.mat';
load 'ComplementaryData/blastedBCC7051Ncrassa.mat';
 
% 2.2 Creation of draft networks from protein orthology inference
% Orthologous proteins were mapped to template models using getModelFromHomology 
% function. At this point, there are 5 draft networks generated for
% BCC7051.
bccDraftFromCordyceps = getModelFromHomology({cordModel},...
    blastedBCC7051Cordyceps,'BCC7051',{},1,false,10^-50,200,40);
bccDraftFromCordyceps.id = 'bccDraftFromCordyceps';

exportForGit(bccDraftFromCordyceps,'bccDraftFromCordyceps',...
    '/Users/nachonase/Documents/GitHub/BCC7051-GSMM/', 'xlsx');

bccDraftFromAoryzae = getModelFromHomology({aoryModel},...
    blastedBCC7051Aoryzae,'BCC7051',{},1,false,10^-50,200,40);
bccDraftFromAoryzae.id = 'bccDraftFromAoryzae';

exportForGit(bccDraftFromAoryzae,'bccDraftFromAoryzae',...
    '/Users/nachonase/Documents/GitHub/BCC7051-GSMM/', 'xlsx');

bccDraftFromAnidulans = getModelFromHomology({anidModel},...
    blastedBCC7051Anidulans,'BCC7051',{},1,false,10^-50,200,40);
bccDraftFromAnidulans.id = 'bccDraftFromAnidulans';

exportForGit(bccDraftFromAnidulans,'bccDraftFromAnidulans',...
    '/Users/nachonase/Documents/GitHub/BCC7051-GSMM/', 'xlsx');

bccDraftFromPenicillium = getModelFromHomology({peniModel},...
    blastedBCC7051Penicillium,'BCC7051',{},1,false,10^-50,200,40);
bccDraftFromPenicillium.id = 'bccDraftFromPenicillium';
exportForGit(bccDraftFromPenicillium,'bccDraftFromPenicillium',...
    '/Users/nachonase/Documents/GitHub/BCC7051-GSMM/', 'xlsx');

bccDraftFromNcrassa = getModelFromHomology({neurModel},...
    blastedBCC7051Ncrassa,'BCC7051',{},1,false,10^-50,200,40);
bccDraftFromNcrassa.id = 'bccDraftFromNcrassa';
exportForGit(bccDraftFromNcrassa,'bccDraftFromNcrassa',...
    '/Users/nachonase/Documents/GitHub/BCC7051-GSMM/', 'xlsx');

% Add Non-gene-associated reactions, including spontaneous, 
% transport and exchange reactions from iNR1329

rxnToAdd.rxns =setdiff(cordModel.rxns,bccDraftFromCordyceps.rxns);
bccDraftFromCordycepsPlusiNR1329 = addRxnsGenesMets(bccDraftFromCordyceps,...
    cordModel,rxnToAdd.rxns,false,...
    'additional rxns from iNR1329 model',2);

finalValidateModel = bccDraftFromCordycepsPlusiNR1329;
finalValidateModel = setParam(finalValidateModel,'eq',{'matp'},1);
finalValidateModel = setParam(finalValidateModel,'obj',{'bmOUT'},1);
model=setParam(finalValidateModel,'ub',{'glcIN'},0.1448);
sol = solveLP(model,1);
fprintf(['umax = ' num2str(sol.f*-1) ' per hour' '\n']);

exportForGit(bccDraftFromCordycepsPlusiNR1329,'bccDraftFromCordycepsPlusiNR1329',...
    '/Users/nachonase/Documents/GitHub/BCC7051-GSMM/', 'xlsx');

% 2.3 Polishing the metabolic network coverage of A. oryzae BCC7051
% In this case, iNR1329 is the only template model that has cordycepin pathway, 
% so the draft initiated from Cordyceps militaris model (iNR1329) is 
% promisingly used as a main network for A. oryzae BCC7051.
% In particular other drafts created from A. nidulans (ANID), Penicillium (iAL1006)
% and N. crassa (MODEL1212060001) as well as an earlier genome-scale metabolic model 
% of  A. oryzae RIB40 (iWV1346) were used as alternative networks 
% supported for enhancing metabolic coverage of the main network. 
% For instance, by checking the EC numbers, some new reactions were added from 
% iAL1006 with annotated genes from alternative networks.
  

curatedModel=importExcelModel('model/xlsx/curated_bccDraftFromCordycepsPlusiNR1329.xlsx');


%% get new gsmm for BCC7051

%% add BCC7051 biomass
[~, SheetS]=xlsread('ComplementaryData/supplementary.xlsx','Table S4');
Biomass = struct();
Biomass.rxns = SheetS(2:end,1);
Biomass.rxnNames = SheetS(2:end,2);
Biomass.equations = SheetS(2:end,3);
Biomass.grRules = SheetS(2:end,4);

optimizedBiomass = addRxns(DraftBCC7051tran,Biomass,2,'c',true,true);
optimizedBiomass = setParam(optimizedBiomass,'lb',Biomass.rxns,...
    [0, 0, 0, 0, 0, 0, 0]);
optimizedBiomass = setParam(optimizedBiomass,'ub',Biomass.rxns,...
    [1000, 1000, 1000, 1000, 1000, 1000, 1000]);

%% add exchange rxns
exchange = getExchangeRxns(cordModel);
DraftBCC7051Plus = addRxnsGenesMets(optimizedBiomass,...
    cordModel,exchange,...
    false,'exchange',2);

%% add required rxns
needToAdd = {'r0001';'r0002';'r0037';...
    'r0130';'r0133';'r0134';'r0145';'r0175';'r0181';'r0186';'r0191';...
    'r0205';'r0219';'r0223';'r0294';'r0297';'r0298';'r0305';'r0460';...
    'r0462';'r0467';'r0468';'r0469';'r0477';'r0489';'r0492';'r0505';...
    'r0507';'r0509';'r0589';'r0591';'r0609';'r0649';'r0694';'r0716';...
    'r1126';'r1154';'r477'}
DraftBCC7051PlusPlus = addRxnsGenesMets(DraftBCC7051Plus,...
    cordModel,needToAdd,...
    false,'needToAdd',2);

DraftBCC7051PlusPlus = setParam(DraftBCC7051PlusPlus,'obj','bmOUT',1);
DraftBCC7051PlusPlus = setParam(DraftBCC7051PlusPlus,'eq','matp',1);
requiredRxns = {'o2IN' 'piIN' 'slfIN' 'nh3IN'};
DraftBCC7051PlusPlus = setParam(DraftBCC7051PlusPlus,'ub',requiredRxns,...
    [1000,1000,1000,1000]);
DraftBCC7051PlusglcIN = setParam(DraftBCC7051PlusPlus,'ub','glcIN',0.2965);

sol = solveLP(DraftBCC7051PlusglcIN);
printFluxes(DraftBCC7051PlusglcIN, sol.x);
%%

DraftBCC7051FromiWV1346 = getModelFromHomology({reduced},...
    blastBCC7051RIB40,'BCC7051',{},2,false,10^-5,100);

genesToAdd=struct();
genesToAdd.genes = setdiff(DraftBCC7051FromiWV1346.genes,DraftBCC7051PlusglcIN.genes);
DraftBCC7051Plusd= addGenesRaven(DraftBCC7051Plus, genesToAdd);

model = addRxnsGenesMets(DraftBCC7051Plusd,DraftBCC7051FromiWV1346,...
    {'r1126';'r1154';'r477'},true,'');

new = contractModel(model);
new=deleteUnusedGenes(new)
newDraftBCC7051PlusglcIN = setParam(new,'ub','glcIN',0.2965);

sol = solveLP(newDraftBCC7051PlusglcIN);
printFluxes(newDraftBCC7051PlusglcIN, sol.x);

exportForGit(newDraftBCC7051PlusglcIN,'newDraftBCC7051PlusglcIN',...
    '',{'xlsx','xml'});

