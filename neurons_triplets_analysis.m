 clear
load('Supplementary_File1all_neurons.mat')
load('Supplementary_File2_tripletsTable.mat')
neuronNames = categorical(extractfield(all_neurons,'NameStr'));

tripletsTableOriginal = tripletsTable;

%%


global proximityThersh; 
global useMean;
%The fraction of shuffels the distance between synapses of two neurons must
%be smaller than to consider the synapses of the two neurons close on the
%common neuron
proximityThersh = 0.05; 

%Use mean distace and not minimal distace between synapses
useMean = true;
%The folder in which figures should be saved
figsSaveFold = 'C:\figs\';
%%
%plot hist of all triplets
figure('units','normalized','outerposition',[0 0 1 1])
fSize = 30;
subplot(1,2,1)
nBins = 20;
histogram(max([tripletsTable.probShuffleB tripletsTable.probShuffleA]'),nBins,'Normalization','probability')
hold on
binCounts = histcounts(max([tripletsTable.probShuffleB tripletsTable.probShuffleA]'),20);
binCounts = binCounts/sum(binCounts);
bar(1/nBins/2,binCounts(1),1/nBins)
xlim([0 1])
hold on

xlabel('Bootstap P-value')
set(gca,'XTick',[0:0.2:1]);
ylabel('Fraction of triplets')
figName = 'allTripletsShufflePV';
set(gca,'FontSize',fSize)


%%
%add normalized mean distance and precentile mean distance to all triplets
mainNsNames = unique(tripletsTable.commonNeuron);
% precentilesMeanDist is place of the mean distace between two neuorons on a
% common neuron compared to all pairwise distances between neurons on the common neuron
% (divided by the number of pairwise distances)
precentilesMeanDist = nan(height(tripletsTable),1);

%normalizedMeanDist is the mean disntace divided by the maximal distance
%between and pair on neurons on the common neuron
normalizedMeanDist = nan(height(tripletsTable),1);
for i=1:length(mainNsNames)
    thisSub = tripletsTable(tripletsTable.commonNeuron==mainNsNames(i),:);
    thisMeanDists = sort(thisSub.meanDist);
    thisInds = find(tripletsTable.commonNeuron==mainNsNames(i));
    nDists = length(thisMeanDists);
    thisN = all_neurons(neuronNames==mainNsNames(i));
    thisMDist = max(thisN.dist_map(:)); 
    for k=1:length(thisInds)
        precentilesMeanDist(thisInds(k)) = find(thisMeanDists==tripletsTable.meanDist(thisInds(k)),1)/length(thisMeanDists);
        normalizedMeanDist(thisInds(k)) = tripletsTable.meanDist(thisInds(k))/thisMDist;
    end
    
end
tripletsTable.precentilesMeanDist =precentilesMeanDist;
tripletsTable.normMeanDist = normalizedMeanDist;
%%
%add distance between cell bodies of n1 and n2 where both exist
cellBodyDists = zeros(height(tripletsTable),1);
for i=1:height(tripletsTable)
    ind1 = find(neuronNames==tripletsTable.N1name(i));
    ind2 = find(neuronNames==tripletsTable.N2name(i));
    
    x1 = all_neurons(ind1).cellBodyX;
    y1 = all_neurons(ind1).cellBodyY;
    z1 = all_neurons(ind1).cellBodyZ;
    
    x2 = all_neurons(ind2).cellBodyX;
    y2 = all_neurons(ind2).cellBodyY;
    z2 = all_neurons(ind2).cellBodyZ;
    
    cellBodyDists(i) = sqrt((x1-x2)^2+(y1-y2)^2+(z1-z2)^2);
end
tripletsTable.cellBodyDistsN1N2 = cellBodyDists;


%%
% right left comparison - find all entries in the table in which n1 and
% n2 are a L/R pair.
nName1 = tripletsTable.N1name;
nName2 = tripletsTable.N2name;
leftRightMeans = [];
otherMean = [];
isPair = zeros(height(tripletsTable),1);
for i=1:height(tripletsTable)
    n1_name = string(tripletsTable.N1name(i));
    n2_name = string(tripletsTable.N2name(i));
    n1Start = extractBefore(n1_name,strlength(n1_name)); 
    n2Start = extractBefore(n2_name,strlength(n2_name));
    if n1Start==n2Start
        n1End = extractAfter(n1_name,strlength(n1_name)-1);
        n2End = extractAfter(n2_name,strlength(n2_name)-1);
        if or(n1End=='R' && n2End=='L',n1End=='L' && n2End=='R')
           isPair(i)=1;
           
        end
       
    end
end
tripletsTable.isPair = isPair;
%%

%add common partners num to table
commonNeighborTableAll = getCommonNeighborTable(tripletsTable);
commonPartnersNum = zeros(height(tripletsTable),1);
for i=1:height(commonNeighborTableAll)
    commonPartnersNum(tripletsTable.N1name==commonNeighborTableAll.n1Name(i) &...
        tripletsTable.N2name==commonNeighborTableAll.n2Name(i)) = commonNeighborTableAll.nNeighbors(i);
    commonPartnersNum(tripletsTable.N2name==commonNeighborTableAll.n1Name(i) &...
        tripletsTable.N1name==commonNeighborTableAll.n2Name(i)) = commonNeighborTableAll.nNeighbors(i);
end
sum(isnan(commonPartnersNum))
tripletsTable.commonPartnersNum = commonPartnersNum;
[pPairNPartners ~] = ranksum(commonNeighborTableAll.nNeighbors(commonNeighborTableAll.isPair==1),commonNeighborTableAll.nNeighbors(commonNeighborTableAll.isPair==0))
disp(['for pairs, the mean number of common partners is ' num2str(mean(commonNeighborTableAll.nNeighbors(commonNeighborTableAll.isPair==1)))])
disp(['for others, the mean number of common partners is ' num2str(mean(commonNeighborTableAll.nNeighbors(commonNeighborTableAll.isPair==0)))])



%%
%Plot correlations (common neighbors plot)
mSize = 30;
fSize = 30;
meanDists = zeros(1,7);
pClose =  zeros(1,7);
allDists = {};
for i=1:6
    thisRows  = tripletsTable(tripletsTable.commonPartnersNum==i,:);
    thisMeans = thisRows.meanDist;
    meanDists(i) = nanmean(thisMeans);
    allDists{i} = thisMeans;
    pClose(i) = sum(thisRows.probShuffleA<proximityThersh & thisRows.probShuffleB<proximityThersh)/height(thisRows);
    
end
thisRows  = tripletsTable(tripletsTable.commonPartnersNum>6,:);
thisMeans = thisRows.meanDist;
meanDists(7) = nanmean(thisMeans);
allDists{7} = thisMeans;
pClose(7) = sum(thisRows.probShuffleA<proximityThersh & thisRows.probShuffleB<proximityThersh)/height(thisRows);
figure('units','normalized','outerposition',[0 0 1 1])

subplot(1,2,1)
plot(1:7,pClose,'.','MarkerSize',mSize)
xlabel('Number of common neighbors')
set(gca,'xtick',[1:7])
set(gca,'xticklabels',{'1','2','3','4','5','6','7+'})
ylabel('Fraction of adjacent pairs')
set(gca,'FontSize',fSize)
xlim([0 8])
pbaspect([1 1 1])
subplot(1,2,2)
plot(1:7,meanDists,'.','MarkerSize',mSize)
set(gca,'xtick',[1:7])
set(gca,'xticklabels',{'1','2','3','4','5','6','7+'})
ylabel('Mean distance between synapses [\mum]')
xlabel('Number of common neighbors')
xlim([0 8])
set(gca,'FontSize',fSize)
pbaspect([1 1 1])


figure('units','normalized','outerposition',[0 0 1 1])
subplot(1,2,1)
plot(1:7,pClose,'.','MarkerSize',mSize)
xlabel('Number of common neighbors')
set(gca,'xtick',[1:7])
set(gca,'xticklabels',{'1','2','3','4','5','6','7+'})
ylabel('Fraction of adjacent pairs')
set(gca,'FontSize',fSize)
xlim([0 8])
subplot(1,2,2)
edges = [0:3:30 1000];
nBins = length(edges)-1;
plotHist = zeros(7,nBins);

for i=1:7
    [N,edges] = histcounts(allDists{i},edges,'Normalization','probability');
    plotHist(i,:) = N;
end
imagesc(plotHist')
set(gca,'YDir','normal')

ticks = {};
for i=1:length(edges)-1
    ticks{i} = num2str(edges(i));
end
ticks{end} = [ticks{end} '+'];
set(gca,'YTick',1:length(ticks))
set(gca,'YTickLabel',ticks)
set(gca,'XTick',1:7)
set(gca,'xticklabels',{'1','2','3','4','5','6','7+'})

xlabel('Number of common neighbors')
ylabel('Mean distance between synapses [\mum]')
set(gca,'FontSize',fSize)

cbh = colorbar;
colormap jet
cbh.Title.String = "Fraction of pairs";


[corrMeanNCommonNeighbors PCorr] = corrcoef(tripletsTable.commonPartnersNum,tripletsTable.meanDist)
%%

% compare all types and dierections to find groups in which the fraction of
% close synapses is significantly higher than expected.
tripletsTableMod = tripletsTable;
tripletsTableMod.mainNeuronClass(tripletsTableMod.mainNeuronClass=='command') = 'interneuron';
tripletsTableMod.N1Class(tripletsTableMod.N1Class=='command') = 'interneuron';
tripletsTableMod.N2Class(tripletsTableMod.N2Class=='command') = 'interneuron';

pvalsTable = table();
types = categorical({'sensory', 'interneuron', 'motorneuron'});
directions = categorical({'fromNeuron', 'toNeuron'});
for iMain = 1:length(types)
    cond.classMain = types(iMain);
    for iN1 =  1:length(types)
        cond.class1 = types(iN1);
        for iN2 = 1:length(types)
            cond.class2 = types(iN2);
            for dN1 = 1:length(directions)
                 cond.dierection1 = directions(dN1);
                for dN2 = 1:length(directions)
                    cond.dierection2 = directions(dN2);
                    if height(pvalsTable)>0
                        %check if the symmetric condition was already
                        %checked
                        if (any((pvalsTable.mainClass==cond.classMain) & (pvalsTable.n1Class==cond.class2) & (pvalsTable.n2Class==cond.class1) & (pvalsTable.n1Dierection==cond.dierection2) & (pvalsTable.n2Dierection==cond.dierection1)))
                        continue
                        end
                    end
                    thisRow = compareGroups(tripletsTableMod,cond,proximityThersh,false);
                    pvalsTable = [pvalsTable; thisRow];
                    
                end
            end
        end
    end
end


%%

%find significant groups
pvThresh = 0.0001;
pvalsTable(isnan(pvalsTable.Pval),:)=[];
comparisonsNum =  height(pvalsTable);
%add bonferroni correction, set max pval to one.
pvalsTable.PvalCorrected = min(pvalsTable.Pval*comparisonsNum,ones(length(pvalsTable.Pval),1));
pvalsTable.PvalCorrectedLogistic = min(pvalsTable.pvLogistic*comparisonsNum,ones(length(pvalsTable.Pval),1));
pvalsTable.PvalCorrectedNormMean = min(pvalsTable.pvNormMean*comparisonsNum,ones(length(pvalsTable.Pval),1));
pvalsTable.PvalCorrectedPercentile = min(pvalsTable.pvPercentile*comparisonsNum,ones(length(pvalsTable.Pval),1));

pvalsTable.pvFinal = pvalsTable.PvalCorrectedLogistic;
 
pvalsTable.pvFinal(pvalsTable.logisticCoeff<0)=1;

figure('units','normalized','outerposition',[0 0 1 1])
semilogy((pvalsTable.pvFinal),'.','markerSize',15);
hold on
sigClose = pvalsTable(pvalsTable.pvFinal<pvThresh,:);
sigCloseInds = find(pvalsTable.pvFinal<pvThresh);
semilogy(sigCloseInds,(pvalsTable.pvFinal(pvalsTable.pvFinal<pvThresh)),'.r','markerSize',15);
mainNSig = cellstr(sigClose.mainClass);
N1Sig = cellstr(sigClose.n1Class);
N2Sig = cellstr(sigClose.n2Class);
sigCloseMod = sigClose;
sigCloseMod.n1Dierection(sigCloseMod.n1Dierection=='fromNeuron') = 'to';
sigCloseMod.n1Dierection(sigCloseMod.n1Dierection=='toNeuron') = 'from';
sigCloseMod.n2Dierection(sigCloseMod.n2Dierection=='fromNeuron') = 'to';
sigCloseMod.n2Dierection(sigCloseMod.n2Dierection=='toNeuron') = 'from';
N1dierection = cellstr(sigCloseMod.n1Dierection);
N2dierection = cellstr(sigCloseMod.n2Dierection);

for i=1:length(mainNSig)
    thisStr = [N1dierection{i} ' ' N1Sig{i} newline ...
        N2dierection{i} ' ' N2Sig{i} newline ...
        'on ' mainNSig{i}];
        text(sigCloseInds(i),(sigClose.pvFinal(i)),thisStr,'FontSize',24)
end
ylim([ min(pvalsTable.pvFinal)/100 1])
xlim([0 length(pvalsTable.pvFinal)+1])
ylabel('Bonferroni corrected p value')
set(gca,'FontSize',24)
set(gca,'YTick',[1e-7 1e-6 1e-5 1e-4 1e-3 1e-2 1e-1 1])
pbaspect([1 1 1])

%%


%%
%plot effect size and pVal

figure('units','normalized','outerposition',[0 0 1 1])
subplot(1,3,3)
errs = sqrt(sigClose.fCond.*(1-sigClose.fCond)./sigClose.Ncond);
errorbar(1:length(errs),sigClose.fCond,errs,'o','MarkerSize',10,'LineWidth',3,...
    'color','black','MarkerFaceColor','black','CapSize',30)
set(gca,'FontSize',30)
hold on
plot([0,length(errs)+1],[mean(pvalsTable.fOthers), mean(pvalsTable.fOthers)],'r','LineWidth',3)
xlim([0, length(errs)+1])
ymax = max(errs+sigClose.fCond)*1.1;
ylim([0,ymax])
xlabel('Condition')
set(gca,'XTickLabel',{})
ylabel('Fraction close')
for i=1:height(sigClose)
    thisStr = [N1dierection{i} ' ' N1Sig{i} newline ...
        N2dierection{i} ' ' N2Sig{i} newline ...
        'on ' mainNSig{i}];
        text(i+0.1,sigClose.fCond(i),thisStr,'FontSize',30)
end
%pbaspect([1 1 1])

subplot(1,3,1:2)
semilogy((pvalsTable.pvFinal),'black.','markerSize',30);
hold on
sigClose = pvalsTable(pvalsTable.pvFinal<pvThresh,:);
sigCloseInds = find(pvalsTable.pvFinal<pvThresh);
semilogy(sigCloseInds,(pvalsTable.pvFinal(pvalsTable.pvFinal<pvThresh)),'.r','markerSize',30);
mainNSig = cellstr(sigClose.mainClass);
N1Sig = cellstr(sigClose.n1Class);
N2Sig = cellstr(sigClose.n2Class);
sigCloseMod = sigClose;
sigCloseMod.n1Dierection(sigCloseMod.n1Dierection=='fromNeuron') = 'to';
sigCloseMod.n1Dierection(sigCloseMod.n1Dierection=='toNeuron') = 'from';
sigCloseMod.n2Dierection(sigCloseMod.n2Dierection=='fromNeuron') = 'to';
sigCloseMod.n2Dierection(sigCloseMod.n2Dierection=='toNeuron') = 'from';
N1dierection = cellstr(sigCloseMod.n1Dierection);
N2dierection = cellstr(sigCloseMod.n2Dierection);

for i=1:length(mainNSig)
    thisStr = [N1dierection{i} ' ' N1Sig{i} newline ...
        N2dierection{i} ' ' N2Sig{i} newline ...
        'on ' mainNSig{i}];
        text(sigCloseInds(i),(sigClose.pvFinal(i)),thisStr,'FontSize',30)
end
ylim([ min(pvalsTable.pvFinal)/100 1])
xlim([0 length(pvalsTable.pvFinal)+1])
ylabel('Bonferroni corrected p value')
set(gca,'FontSize',30)
set(gca,'YTick',[1e-7 1e-6 1e-5 1e-4 1e-3 1e-2 1e-1 1])
pbaspect([1 1 1])


% 



%%
%look for correlation between number of common neighbors and distance between
%synapses
     
% set the threshold for neighbors before binning the rest
maxN=7;
allNs = 1:max(commonNeighborTableAll.nNeighbors);  
allNsMeanAndProbsAll = [];
allNsMeanDists = [];
numberOfMeasurmentsAll = [];
allNsSigA = [];
allNsMeanAndProbsInIn= [];
numberOfMeasurmentsInIn = [];

allNsMeanAndProbsOutOut = [];
numberOfMeasurmentsOutOut = [];

allNsMeanAndProbsInOut = [];
numberOfMeasurmentsInOut = [];

allMeanDists = [];

for i=1:length(allNs)
    allNsSigA = [allNsSigA; mean(commonNeighborTableAll.aFrac(commonNeighborTableAll.nNeighbors==allNs(i)))];    
    allNsMeanDists = [allNsMeanDists; mean(commonNeighborTableAll.meanDist(commonNeighborTableAll.nNeighbors==allNs(i)))];

    allNsMeanAndProbsAll = [allNsMeanAndProbsAll; mean(commonNeighborTableAll.AandBFrac(commonNeighborTableAll.nNeighbors==allNs(i)))];
    numberOfMeasurmentsAll = [numberOfMeasurmentsAll; sum(commonNeighborTableAll.nNeighbors==allNs(i))];
    
    allMeanDists = [allMeanDists; mean(commonNeighborTableAll.meanDist(commonNeighborTableAll.nNeighbors==allNs(i)))];
end
allNsMeanAndProbsAllMod = [];
allNsMeanAndProbsAllMod(1:maxN-1) = allNsMeanAndProbsAll(1:maxN-1);
allNsMeanAndProbsAllMod(end+1) = nansum(allNsMeanAndProbsAll(maxN:end).*numberOfMeasurmentsAll(maxN:end))/nansum(numberOfMeasurmentsAll(maxN:end));

allMeanDistsMod = [];
allMeanDistsMod(1:maxN-1) = allMeanDists(1:maxN-1);
allMeanDistsMod(end+1) = nansum(allMeanDists(maxN:end).*numberOfMeasurmentsAll(maxN:end))/nansum(numberOfMeasurmentsAll(maxN:end));

lastPoint = maxN;

allNsMeanAndProbsAll(end) = sum(allNsMeanAndProbsAll(lastPoint:end).*numberOfMeasurmentsAll(lastPoint:end))/sum(numberOfMeasurmentsAll(lastPoint:end));

allNsMeanAndProbsAll(end) = sum(allNsMeanAndProbsAll(lastPoint:end).*numberOfMeasurmentsAll(lastPoint:end))/sum(numberOfMeasurmentsAll(lastPoint:end));

figure('units','normalized','outerposition',[0 0 1 1])
subplot(1,2,1)
plot(1:maxN,allNsMeanAndProbsAllMod,'.','MarkerSize',25)
hold all

populationFracAandB = sum(tripletsTable.probShuffleA<proximityThersh & tripletsTable.probShuffleB<proximityThersh)/height(tripletsTable);
%plot([1 maxN],[populationFracAandB; populationFracAandB],'--')
pbaspect([1 1 1])
xlabel('Number of common neighbors')
ylabel('Fraction of partners that are close')
xlim([0,maxN+1])
set(gca,'XTick',1:maxN)
set(gca,'XTickLabel',[split(num2str(1:maxN-1)); '>=' num2str(maxN)])
set(gca,'FontSize',24)
%legend('','Population Mean','Location','northwest')
[correlationsFrac PValsFrac] = corrcoef(1:maxN,allNsMeanAndProbsAllMod)
subplot(1,2,2)
plot(1:maxN,allMeanDistsMod,'.','MarkerSize',25)
pbaspect([1 1 1])

xlabel('Number of common neighbors')

ylabel('Mean min distance between pairs [um]')

if useMean==true
    ylabel('Mean distance between pairs [um]')
end

xlim([0,maxN+1])
set(gca,'XTick',1:maxN)

set(gca,'XTickLabel',[split(num2str(1:maxN-1)); '>=' num2str(maxN)])
set(gca,'FontSize',24)
[correlationsDist PValsDist] = corrcoef(1:maxN,allMeanDistsMod)




%%

function  commonNeighborTable= getCommonNeighborTable(tripletsTable)
% returns a table 

%get all the neurons names
names = unique(tripletsTable.N1name);
%global variable that sets the threshold for a pair to be considered close
global proximityThersh;
commonNeighborTable = table();
% go over all the neuron names
for n1=1:length(names)
    % for n2, go over all the pairs n1,2 not reviewed yet
    for  n2 = n1+1:length(names)
        %make sure these are different neurons
        if names(n1)~=names(n2)
            % look at all rows of neurons n1,n2
            thisSubTable = tripletsTable(or(tripletsTable.N1name==names(n1)&tripletsTable.N2name==names(n2),...
                tripletsTable.N2name==names(n1)&tripletsTable.N1name==names(n2)),:);
            if height(thisSubTable)>0
                %count the number of common neighbors
                thisN = length(unique(thisSubTable.commonNeuron));
                %             thisAandBSigProb =  sum(thisSubTable.probShuffleA<proximityThersh)/height(thisSubTable);
                %             thisBSigProb =  sum(thisSubTable.probShuffleB<proximityThersh)/height(thisSubTable);
                %find the proportion of close pairs
                thisAandBSigProb =  sum(thisSubTable.probShuffleB<proximityThersh & thisSubTable.probShuffleA<proximityThersh...
                    )/height(thisSubTable);
                thisASigProb = sum(thisSubTable.probShuffleA<proximityThersh)...
                    /height(thisSubTable);
           
                
                %get the mean distance
                thisMeanDist = mean([thisSubTable.meanDist ; thisSubTable.meanDist]);
                isPair = mean(thisSubTable.isPair);
                
                thisRow = table(thisN,thisAandBSigProb,thisASigProb,thisMeanDist,...
                    names(n1),names(n2),n1,n2,isPair,'VariableNames',{'nNeighbors','AandBFrac','aFrac','meanDist'...
                    ,'n1Name','n2Name','n1Ind','n2Ind','isPair'});
                commonNeighborTable = [commonNeighborTable; thisRow];
            end
        end
    end
end
end

        

