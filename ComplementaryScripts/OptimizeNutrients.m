%% load model
modelFileName = 'C:\Users\Lenovo\Documents\GitHub\BCC7051-GSMM\model\xml\newOryzae.xml';
model = readCbModel(modelFileName);
coolModel = ravenCobraWrapper(model);
coolModel = importModel('C:\Users\Lenovo\Documents\GitHub\BCC7051-GSMM\model\xml\newOryzae.xml');
model = coolModel;

model = setParam(model,'ub',{'cordycepinOUT'},1000);
model = setParam(model,'ub',{'bmOUT'},1000);

%% constraint
uptake = {};

% Loop through each reaction ID
for i = 1:length(model.rxns)
    rxnID = model.rxns{i};
    if endsWith(rxnID, 'IN', 'IgnoreCase', false)  % case sensitive
        uptake{end+1} = rxnID; %#ok<SAGROW>
    end
end

%% Optimize requirement nutrients for growth

for 

modelSenerioIcI=setParam(coolModel,'ub',{'glcIN'},1000);
modelSenerioIcI=setParam(modelSenerioIcI,'eq',{'R04591_c'},0);
modelSenerioIcI = setParam(modelSenerioIcI,'eq',{'bmOUT'},1);
modelSenerioIcI=setParam(modelSenerioIcI,'ub',{'cordycepinOUT'},1000);
modelSenerioIcI = setParam(modelSenerioIcI,'obj',{'cordycepinOUT'},1);
solmodelSenerioIcI = solveLP(modelSenerioIcI,1);
fprintf(['umax = ' num2str(solmodelSenerioIcI.f*-1) ' per hour' '\n']);
printFluxes(modelSenerioIcI, solmodelSenerioIcI.x, true)


targets=FSEOF(modelSenerioII,"bmOUT","cordycepinOUT",100)


modelC = ravenCobraWrapper(modelSenerioIcI);
[grRatio, grRateKO, grRateWT, hasEffect, delRxn, fluxSolution] = singleRxnDeletion(modelC);



coolModel = importModel('C:\Users\Lenovo\Documents\GitHub\BCC7051-GSMM\model\xml\newOryzae.xml');
modelSenerioIcI=setParam(coolModel,'ub',{'glcIN'},1000);
modelSenerioIcI=setParam(modelSenerioIcI,'ub',{'cordycepinOUT'},1000);
modelSenerioIcI = setParam(modelSenerioIcI,'obj',{'bmOUT'},1);
solmodelSenerioIcI = solveLP(modelSenerioIcI,1);
fprintf(['umax = ' num2str(solmodelSenerioIcI.f*-1) ' per hour' '\n']);
printFluxes(modelSenerioIcI, solmodelSenerioIcI.x, true)
targets=FSEOF(modelSenerioIcI,"bmOUT","cordycepinOUT",100)


modelC = ravenCobraWrapper(modelSenerioIcI);
[grRatio, grRateKO, grRateWT, hasEffect, delRxn, fluxSolution] = singleRxnDeletion(modelC);



%% Experimental senerio fitting
% load model
cd 'C:\Users\Lenovo\Documents\GitHub\BCC7051-GSMM';
%
coolModel = importModel('model/xml/newOryzae.xml',false);
coolModel.equations = constructEquations(coolModel);
%% Senerio I: wild-type strain grew in glucose
modelSenerioI=setParam(coolModel,'ub',{'glcIN'},0.632);
modelSenerioI = setParam(modelSenerioI,'obj',{'bmOUT'},1);
solmodelSenerioI = solveLP(modelSenerioI,1);
fprintf(['umax = ' num2str(solmodelSenerioI.f*-1) ' per hour' '\n']);
printFluxes(modelSenerioI, solmodelSenerioI.x, true)
%Experiment = 0.0317 ± 0.0036 h-1

targets=FSEOF(modelSenerioI,"bmOUT","cordycepinOUT",100)

modelC = ravenCobraWrapper(modelSenerioI);
[grRatio, grRateKO, grRateWT, hasEffect, delRxn, fluxSolution] = singleRxnDeletion(modelC);

%% Senerio II: cordyceping-producing strain grew in glucose
modelSenerioII=setParam(coolModel,'ub',{'glcIN'},0.5487);
modelSenerioII=setParam(modelSenerioII,'eq',{'cordycepinOUT'},0.0126);
modelSenerioII = setParam(modelSenerioII,'obj',{'bmOUT'},1);
solmodelSenerioII = solveLP(modelSenerioII,1);
fprintf(['umax = ' num2str(solmodelSenerioII.f*-1) ' per hour' '\n']);
printFluxes(modelSenerioII, solmodelSenerioII.x, true)
%Experiment = 0.0253 ± 0.0019 h-1


modelSenerioIcI=setParam(modelSenerioIcI,'eq',{'R04591_c'},0);


targets=FSEOF(modelSenerioII,"bmOUT","cordycepinOUT",100)

modelC = ravenCobraWrapper(modelSenerioII);
[grRatio, grRateKO, grRateWT, hasEffect, delRxn, fluxSolution] = singleRxnDeletion(modelC);


%% Analysis different fluxes for cordycepin production by A.oryzae strains
followChanged(coolModel,solmodelSenerioII.x,solmodelSenerioI.x, 95)



modelSenerioII7=setParam(coolModel,'ub',{'glcIN'},0.2547);
modelSenerioII7 = setParam(modelSenerioII7,'obj',{'cordycepinOUT'},1);
solmodelSenerioII7 = solveLP(modelSenerioII7,1);
fprintf(['cordycepin = ' num2str(solmodelSenerioII7.f*-1) ' per hour' '\n']);


%% Vary nutrients
modelBigger = coolModel;
sol = solveLP(modelBigger,1);
printFluxes(modelBigger,sol.x);

%% variedResultfoRGROWTH
modelRaven = modelBigger;
uptakeRxns = modelRaven.rxns(strncmp(modelRaven.rxnNames,'uptake ',7));
uptakeRxnNames = modelRaven.rxnNames(strncmp(modelRaven.rxnNames,'uptake ',7));
requiredRxns = {'o2IN' 'piIN' 'slfIN'};

modelRaven = setParam(modelRaven,'ub',uptakeRxns,0); % block all uptake
modelRaven = setParam(modelRaven,'ub',requiredRxns,1000); % allow only required metabolites
modelRaven = setParam(modelRaven,'lb',{'bmOUT' 'cordycepinOUT'},0);
modelRaven = setParam(modelRaven,'ub',{'bmOUT' 'cordycepinOUT'},1000);
modelRaven = setParam(modelRaven,'obj',{'bmOUT'},1);

%modelRaven = setParam(modelRaven,'ub',{'glcIN','nh3IN'},1000);
%sol = solveLP(modelRaven,1);
%printFluxes(modelRaven,sol.x);

mVersatile_I = cell(numel(uptakeRxns),1);
solVersatile_I= cell(numel(uptakeRxns),1);
model_mVersatile_I = modelRaven; % allow only required metabolites 
for i = 1:numel(uptakeRxns)
    model = model_mVersatile_I;
    model_i = setParam(model,'ub',uptakeRxns{i},1000);
    sol = solveLP(model_i,1);
    mVersatile_I{i} = (sol.f*-1);
    solVersatile_I{i} = sol;
end

solVersatile_II= cell(numel(uptakeRxns),1);
mVersatile_II = cell(numel(uptakeRxns),1);
model_mVersatile_II = setParam(modelRaven,'ub','glcIN',1); % allow only required metabolites and glucose
for i = 1:numel(uptakeRxns)
    model = model_mVersatile_II;
    model_i = setParam(model,'ub',uptakeRxns{i},1000);
    sol = solveLP(model_i,1);
    mVersatile_II{i} = (sol.f*-1);
    solVersatile_II{i} = sol;
end


solVersatile_III = cell(numel(uptakeRxns),1);
mVersatile_III = cell(numel(uptakeRxns),1);
model_mVersatile_III = setParam(modelRaven,'ub','nh3IN',1); % allow only required metabolites and ammonia
for i = 1:numel(uptakeRxns)
    model = model_mVersatile_III;
    model_i = setParam(model,'ub',uptakeRxns{i},1000);
    sol = solveLP(model_i,1);
    mVersatile_III{i} = (sol.f*-1);
    solVersatile_III{i} = sol;
end

variedResultfoRGROWTH = table(uptakeRxns,uptakeRxnNames,mVersatile_I,mVersatile_II,mVersatile_III);

%%variedResultfoRcorDYCEPING

modelRaven = setParam(modelRaven,'eq','bmOUT',0.05);
modelRaven = setParam(modelRaven,'obj',{'cordycepinOUT'},1);

mVersatile_I = cell(numel(uptakeRxns),1);
solVersatile_I= cell(numel(uptakeRxns),1);
model_mVersatile_I = modelRaven; % allow only required metabolites 
for i = 1:numel(uptakeRxns)
    model = model_mVersatile_I;
    model_i = setParam(model,'ub',uptakeRxns{i},1000);
    sol = solveLP(model_i,1);
    mVersatile_I{i} = (sol.f*-1);
    solVersatile_I{i} = sol;
end

solVersatile_II= cell(numel(uptakeRxns),1);
mVersatile_II = cell(numel(uptakeRxns),1);
model_mVersatile_II = setParam(modelRaven,'ub','glcIN',1); % allow only required metabolites and glucose
for i = 1:numel(uptakeRxns)
    model = model_mVersatile_II;
    model_i = setParam(model,'ub',uptakeRxns{i},1000);
    sol = solveLP(model_i,1);
    mVersatile_II{i} = (sol.f*-1);
    solVersatile_II{i} = sol;
end


solVersatile_III = cell(numel(uptakeRxns),1);
mVersatile_III = cell(numel(uptakeRxns),1);
model_mVersatile_III = setParam(modelRaven,'ub','nh3IN',1); % allow only required metabolites and ammonia
for i = 1:numel(uptakeRxns)
    model = model_mVersatile_III;
    model_i = setParam(model,'ub',uptakeRxns{i},1000);
    sol = solveLP(model_i,1);
    mVersatile_III{i} = (sol.f*-1);
    solVersatile_III{i} = sol;
end

variedResultfoRcorDYCEPING = table(uptakeRxns,uptakeRxnNames,mVersatile_I,mVersatile_II,mVersatile_III);


