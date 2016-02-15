function Ind = VarIndx(Values)
global TrialInfo Analyzer  NEV NS6 NS3
Parmass = Analyzer.L.param;
Ind = 1;
for i = 1:length(Values)
    VarValDims = 1;
    for h = 1:i-1
            VarValDims = VarValDims*(length(eval(Parmass{1,h}{2})));
    end
    Ind = Ind + (find(eval(Parmass{1,i}{2}) == Values(i))-1)*VarValDims;
end