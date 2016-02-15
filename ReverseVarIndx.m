function Inds = ReverseVarIndx(Values)
global TrialInfo Analyzer  NEV NS6 NS3 Ind
Parmass = Analyzer.L.param;
In = Values-1;
for n = 1:length(Parmass(1,:))
    i = length(Parmass(1,:))+1-n;
    VarValDims = 1;
    for h = 1:i-1
            VarValDims = VarValDims*(length(eval(Parmass{1,h}{2})));
    end
    Ind = floor(In/VarValDims);
    In = In - Ind*VarValDims;
    Inds(i) = Ind;
end
Inds = Inds+1;
for i = 1:length(Inds)
    tmp = eval(Parmass{1,i}{2});
    Inds(i) = tmp(Inds(i));
end