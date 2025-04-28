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
%% Senerio II: cordyceping-producing strain grew in glucose
modelSenerioII=setParam(coolModel,'ub',{'glcIN'},0.5487);
modelSenerioII=setParam(modelSenerioII,'eq',{'cordycepinOUT'},0.0126);
modelSenerioII = setParam(modelSenerioII,'obj',{'bmOUT'},1);
solmodelSenerioII = solveLP(modelSenerioII,1);
fprintf(['umax = ' num2str(solmodelSenerioII.f*-1) ' per hour' '\n']);
printFluxes(modelSenerioII, solmodelSenerioII.x, true)
%Experiment = 0.0253 ± 0.0019 h-1
%% Analysis different fluxes for cordycepin production by A.oryzae strains
followChanged(coolModel,solmodelSenerioII.x,solmodelSenerioI.x, 95)



modelSenerioII7=setParam(coolModel,'ub',{'glcIN'},0.2547);
modelSenerioII7 = setParam(modelSenerioII7,'obj',{'cordycepinOUT'},1);
solmodelSenerioII7 = solveLP(modelSenerioII7,1);
fprintf(['cordycepin = ' num2str(solmodelSenerioII7.f*-1) ' per hour' '\n']);


