function [row] = compareGroups(tripletsTable,cond,proximityThersh,ignoreDierection)

d1 = cond.dierection1;
d2 = cond.dierection2;
c1 = cond.class1;
c2 = cond.class2;
cM = cond.classMain;

if isempty(ignoreDierection)
    ignoreDierection = false;
end
    
if ignoreDierection==false
    condInds1 = tripletsTable.N1Class==c1 & tripletsTable.N2Class==c2 &...
        tripletsTable.synapseType1==d1 & tripletsTable.synapseType2==d2 & tripletsTable.mainNeuronClass == cM;
    %reversing also works, a->c, b->c in equal to b->c, a->c. a->c, b<-c is
    %equal to b<-c, a->c
    condInds2 = tripletsTable.N1Class==c2 & tripletsTable.N2Class==c1 &...
        tripletsTable.synapseType1==d2 & tripletsTable.synapseType2==d1 & tripletsTable.mainNeuronClass == cM;
    condInds = condInds1 | condInds2;

    
    
    
else
    condInds = tripletsTable.N1Class==c1 & tripletsTable.N2Class==c2 &...
    tripletsTable.mainNeuronClass == cM;
end

thisIsSig = and((tripletsTable.probShuffleA<proximityThersh),tripletsTable.probShuffleB<proximityThersh);


thisType = zeros(length(thisIsSig),1);
thisType(condInds) = 1;
[tbl,~,Pval] = crosstab(thisType & ~tripletsTable.isPair,thisIsSig & ~tripletsTable.isPair);


if isnan(Pval)
    Ncond = nan;
    Nothers = nan;
    fCond = nan;
    fOthers = nan;
    pvLogistic = nan;
    pvNormMean =nan;
    condNormMean = nan;
    othersNormMean= nan;
    condNormSTD= nan;
    logisticCoeff = nan;
    
    pvPercentile = nan;
    condPercentileMean = nan;
    condPercentileSTD = nan;
    othersPercentileMean = nan;
    coeffPercentile = nan;
    logisticCoeffSE = nan;
    coeffPercentileSE = nan;
else
    Ncond = sum(tbl(2,:));
    Nothers = sum(tbl(1,:));
    fCond = tbl(2,2)/Ncond;
    fOthers = tbl(1,2)/Nothers;
    yCond = tbl(2,2);
    yOthers = tbl(1,2);
    %Pval =  proportionTest(yCond,Ncond,yOthers,Nothers,'right');

    isTo1 = tripletsTable.synapseType1=='toNeuron';
    isTo1 = (isTo1 - mean(isTo1))/std(isTo1);
    isTo2 = tripletsTable.synapseType2=='toNeuron';
    isTo2 = (isTo2 - mean(isTo2))/std(isTo2);
    isInter = tripletsTable.mainNeuronClass=='interneuron';
    isInter = (isInter-mean(isInter))/std(isInter);
    isMotor = tripletsTable.mainNeuronClass=='motorneuron';
    isMotor = (isMotor-mean(isMotor))/std(isMotor);
    isSensory = tripletsTable.mainNeuronClass=='sensory';
    X = [thisType tripletsTable.nSynapses1 tripletsTable.nSynapses2 ...
       tripletsTable.commonPartnersNum isTo1 isTo2 isInter isMotor ];

     
    y = thisIsSig;
    
    [B,dev,stats] = mnrfit(X,categorical(~y));
    pvLogistic = stats.p(2);
    logisticCoeff = B(2);
    logisticCoeffSE = stats.se(2);
    %multi linear fit
    mdl = fitlm( X,tripletsTable.normMeanDist);    
    pvNormMean = mdl.Coefficients.pValue(2);
    condNormMean = nanmean(tripletsTable.normMeanDist(thisType==1));
    condNormSTD = nanstd(tripletsTable.normMeanDist(thisType==1));
    othersNormMean = nanmean(tripletsTable.normMeanDist(thisType==0));
    
    
    mdl = fitlm( X,tripletsTable.precentilesMeanDist);    
    pvPercentile = mdl.Coefficients.pValue(2);
    coeffPercentile = mdl.Coefficients.Estimate(2);
    coeffPercentileSE = mdl.Coefficients.SE(2);
    condPercentileMean = nanmean(tripletsTable.precentilesMeanDist(thisType==1));
    if condPercentileMean==1
        disp(cond)
    end
    condPercentileSTD = nanstd(tripletsTable.precentilesMeanDist(thisType==1));
    othersPercentileMean = nanmean(tripletsTable.precentilesMeanDist(thisType==0));
    
end




row = table(cond.classMain,cond.class1, cond.class2, cond.dierection1,cond.dierection2,fCond,fOthers,Pval,Ncond,Nothers,pvLogistic,logisticCoeff,pvNormMean,condNormMean,othersNormMean,condNormSTD,...
    pvPercentile, condPercentileMean,condPercentileSTD, othersPercentileMean,coeffPercentile,logisticCoeffSE,coeffPercentileSE, 'VariableNames',...
    {'mainClass','n1Class','n2Class','n1Dierection','n2Dierection','fCond','fOthers','Pval','Ncond','Nothers','pvLogistic','logisticCoeff','pvNormMean','condNormMean','othersNormMean','condNormSTD'...
    ,'pvPercentile', 'condPercentileMean','condPercentileSTD', 'othersPercentileMean','coeffPercentile','logisticSE','coeffPercentileSE'});




end

