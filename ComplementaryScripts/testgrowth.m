
testgrowth = struct();
testgrowth.umax = cell(numel(goodModel22.rxns),1);
testgrowth.rxns = cell(numel(goodModel22.rxns),1);
for i = 1:numel(goodModel22.rxns)
    modelS = setParam(goodModel22,'lb',goodModel22.rxns{i},-1000);
    modelS = setParam(modelS,'obj','bmOUT',1);
    solS = solveLP(modelS,1);
    testgrowth.umax{i} = solS.f*-1;
    testgrowth.rxns{i} = goodModel22.rxns{i};
end
