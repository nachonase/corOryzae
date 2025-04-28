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
aoryModel = importModel('ComplementaryData/iWV1346.xml');
% 1.2 Model templates prepraration and metabolite naming
[~, textData1]=xlsread('ComplementaryData/supplementary.xlsx','Table S2');
metNames = struct();
metNames.old = textData1(4:end,2);
metNames.new = textData1(4:end,3);
[a, b]=ismember(aoryModel.metNames,metNames.old);
I=find(a);
aoryModel.metNames(I)=metNames.new(b(I));
reduced=contractModel(aoryModel);
reduced = setParam(reduced,'ub',{'r2205'},[1]);
reduced = setParam(reduced,'obj',{'r2359'},1);
sol = solveLP(reduced);
printFluxes(reduced, sol.x, true);

%% STEP 2: DRAFT RECONSTRUCTION
% 2.1 Identification of orthologous proteins 
% A set of orthologous proteins between A. oryzae BCC7051 and template model 
% proteins was obtained by sequence-alignment analysis using BLASTp.
%blastedBCC7051Aoryzae=getBlast('BCC7051','ComplementaryData/BCC7051.fasta',...
%    'iWV1346','ComplementaryData/A_oryzae.faa');
%save('ComplementaryData/blastedBCC7051Aoryzae.mat','blastedBCC7051Aoryzae');

load 'ComplementaryData/blastedBCC7051Aoryzae.mat';
% 2.2 Creation of draft networks from protein orthology inference
% Orthologous proteins were mapped to template models using getModelFromHomology 
% function. At this point, there are 5 draft networks generated for
% BCC7051.
bccDraftFromAoryzae = getModelFromHomology({reduced},...
    blastedBCC7051Aoryzae,'BCC7051',{},1,false,10^-50,200,40);
bccDraftFromAoryzae.id = 'bccDraftFromAoryzae';

run =  importModel('model/xml/newDraftBCC7051PlusglcIN.xml')
model=mergeModels({run,bccDraftFromAoryzae})

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

exportForGit(reduced5,'reduced5',...
    '/Users/nachonase/Documents/GitHub/BCC7051-GSMM/', 'xlsx');

