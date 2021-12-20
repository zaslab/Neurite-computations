% This script computes the number of close sections and the number of
% synapses between all neurons and compares them
%It then computes for each pair the aggregation of synapses

load('Supplementary_File1all_neurons.mat')
dz = 0.08; %distance between EM images is 80 nm


%%
allNSynapses = all_neurons(1).all_synapses;
for i=2:length(all_neurons)
    allNSynapses = [allNSynapses; all_neurons(i).all_synapses];
end
%%

allUniqueSynapses = unique(allNSynapses.idx);


dLsChem = nan(height(allNSynapses),1);
count =1;
sectionsNumChem =  nan(height(allNSynapses),1);
for ii = 1:length(allUniqueSynapses)
    %For each, find all instances
    thisSynapse = allNSynapses(allNSynapses.idx==allUniqueSynapses(ii),:);
    %Ignore electical synapses
    if thisSynapse.dierection(1)~='electrical'
        %Seperate pre from post parts
        thisPre = thisSynapse(thisSynapse.dierection=='fromNeuron',:);
        thisPre = thisPre(1,:);
        thisPreXYZ = [thisPre.x; thisPre.y; thisPre.z];
        thisPosts =  thisSynapse(thisSynapse.dierection=='toNeuron',:);
        %Compute the distance between the pre part and each of the post
        %parts
        for jj = 1:height(thisPosts)
            thisPostXYZ = [thisPosts.x(jj); thisPosts.y(jj); thisPosts.z(jj)];
            dLsChem(count) = sqrt(sum((thisPreXYZ - thisPostXYZ).^2));
            %Also save the size of the synapse
            sectionsNumChem(count) = thisPre.sections;
            count = count+1;
        end
    end
    
end
dLsChem = dLsChem(~isnan(dLsChem));
sectionsNumChem = sectionsNumChem(~isnan(sectionsNumChem));

%Repeat for electrical syns only
%%
dLsElec = nan(height(allNSynapses),1);
count =1;
sectionsNumElec =  nan(height(allNSynapses),1);
for ii = 1:length(allUniqueSynapses)
    %For each, find all instances
    thisSynapse = allNSynapses(allNSynapses.idx==allUniqueSynapses(ii),:);
    
    %Ignore electical synapses
    if thisSynapse.dierection(1)=='electrical' & height(thisSynapse)>1
        %Seperate pre from post parts
        thisPre = thisSynapse(1,:);
        
        thisPreXYZ = [thisPre.x; thisPre.y; thisPre.z];
        thisPosts =  thisSynapse(2:end,:);
        %Compute the distance between the pre part and each of the post
        %parts
        for jj = 1:height(thisPosts)
            thisPostXYZ = [thisPosts.x(jj); thisPosts.y(jj); thisPosts.z(jj)];
            dLsElec(count) = sqrt(sum((thisPreXYZ - thisPostXYZ).^2));
            %Also save the size of the synapse
            sectionsNumElec(count) = thisPre.sections;
            count = count+1;
        end
    end
    
end

dLsElec = dLsElec(~isnan(dLsElec));
sectionsNumElec = sectionsNumElec(~isnan(sectionsNumElec));

%%
figure('units','normalized','outerposition',[0 0 1 1])
trimThresh = 1.5;
minSynDist = median(dLsChem);
dLThresh = quantile(dLsChem,0.5);
dLsChemTrim = dLsChem;
dLsChemTrim(dLsChemTrim>trimThresh) = trimThresh;
histogram(dLsChemTrim,[0:0.05:trimThresh],'Normalization','probability')
xlabel('Distance between synapse ends [um]')
ylabel('Fraction of synapses')
set(gca,'FontSize',28)
xlim([-0 trimThresh])
hold on
yl = ylim;
h1 = plot([dLThresh, dLThresh],yl,'--r','LineWidth',4);
pbaspect([1 1 1])
legend(h1,'Median')
set(gca,'YTick',0:0.04:0.2)

set(gca,'XTick',0:0.5:trimThresh)
set(gca,'XTickLabel',{'0','0.5','1',[num2str(trimThresh) '+']})



%%
neuronsPairsSum = table();

for n1 = 1:length(all_neurons)
    thisN1 = all_neurons(n1);
    for n2 = (n1+1):length(all_neurons)
        if ~(n1==n2)
           thisN2 = all_neurons(n2);
           neuronsPairsSum = [neuronsPairsSum; computePair(thisN1, thisN2,dLThresh)];
        end
    end
    disp(n1);
end

%%
closeSectionsN = unique(neuronsPairsSum.nCloseSections);
nSynsMean = [];
for i=1:length(closeSectionsN)
    nSynsMean(i) = mean(neuronsPairsSum.nSynapses(neuronsPairsSum.nCloseSections==closeSectionsN(i)));
end
normalSectionsNumber = quantile(neuronsPairsSum.nCloseSections(neuronsPairsSum.nSynapses>0),0.5)
isClose = neuronsPairsSum.nCloseSections>normalSectionsNumber;
isCon = neuronsPairsSum.nSynapses>0;
conAndClose = sum(isClose&isCon)
closeAndNotCon = sum(isClose&~isCon)

crosstab(isClose,isCon)
plot(closeSectionsN,nSynsMean,'.')
neuronsPairsSumOrig = neuronsPairsSum;
%%

%%
%bin by the number of close sections.
binSize = 100;
dotSize= 60;
neuronsPairsSum = neuronsPairsSumOrig;
%Sort by the number of close sections
neuronsPairsSum = neuronsPairsSum(~isnan(neuronsPairsSum.cellBodyDist),:);
[~, sortedInds] = sort(neuronsPairsSum.nCloseSections);
neuronsPairsSumAll = neuronsPairsSumOrig;
[~, sortedIndsAll] = sort(neuronsPairsSumAll.nCloseSections);

counter = 0; %count for the number of samples in each bin 
groupSyns = []; %temp for saving the number of syns in each bin
groupSections = []; % temp for saving the number of sections in each bin
groupDists =[]; % temp for the distances in each bin

pSynGroupSec = []; %synapse probability per bin
nSynGroupSec = []; %synapses mean number per bin
meanDistsGroupSec = []; %mean dist per bin
nSectionsGroupSec = [];
groupCounter = 0;
for i=1:length(sortedInds)
    counter = counter + 1;
    groupSyns = [groupSyns neuronsPairsSum.nSynapses(sortedInds(i))];
    groupSections = [groupSections neuronsPairsSum.nCloseSections(sortedInds(i))];
    groupDists = [groupDists neuronsPairsSum.cellBodyDist(sortedInds(i))];
    if counter>=binSize
        counter = 0;
        groupCounter = groupCounter+1;
        
        pSynGroupSec = [pSynGroupSec sum(groupSyns>0)/length(groupSyns)];
        nSynGroupSec = [nSynGroupSec mean(groupSyns)];
        meanDistsGroupSec = [meanDistsGroupSec mean(groupDists)];
        nSectionsGroupSec = [nSectionsGroupSec mean(groupSections)];
        groupSyns = [];
        groupSections = [];
        groupDists =[];
    end
end
pSynGroupSec(end) = (pSynGroupSec(end)*binSize+sum(groupSyns>0))/(length(groupSyns)+binSize);
nSynGroupSec(end) = (nSynGroupSec(end)*binSize + sum(groupSyns))/(length(groupSyns)+binSize);
meanDistsGroupSec(end) = (meanDistsGroupSec(end)*binSize + sum(groupDists))/(length(groupDists)+binSize);

%Repeat for all synapses
counter = 0; %count for the number of samples in each bin 
groupSynsAll = []; %temp for saving the number of syns in each bin
groupSectionsAll = []; % temp for saving the number of sections in each bin
groupDistsAll =[]; % temp for the distances in each bin

pSynGroupSecAll = []; %synapse probability per bin
nSynGroupSecAll = []; %synapses mean number per bin
meanDistsGroupSecAll = []; %mean dist per bin
nSectionsGroupSecAll = [];
groupCounterAll = 0;
for i=1:length(sortedIndsAll)
    counter = counter + 1;
    groupSynsAll = [groupSynsAll neuronsPairsSumAll.nSynapses(sortedIndsAll(i))];
    groupSectionsAll = [groupSectionsAll neuronsPairsSumAll.nCloseSections(sortedIndsAll(i))];
    groupDistsAll = [groupDistsAll neuronsPairsSumAll.cellBodyDist(sortedIndsAll(i))];
    if counter>=binSize
        counter = 0;
        groupCounterAll = groupCounterAll+1;
        
        pSynGroupSecAll = [pSynGroupSecAll sum(groupSynsAll>0)/length(groupSynsAll)];
        nSynGroupSecAll = [nSynGroupSecAll mean(groupSynsAll)];
        meanDistsGroupSecAll = [meanDistsGroupSecAll mean(groupDistsAll)];
        nSectionsGroupSecAll = [nSectionsGroupSecAll mean(groupSectionsAll)];
        groupSynsAll = [];
        groupSectionsAll = [];
        groupDistsAll =[];
    end
end
pSynGroupSecAll(end) = (pSynGroupSecAll(end)*binSize+sum(groupSynsAll>0))/(length(groupSynsAll)+binSize);
nSynGroupSecAll(end) = (nSynGroupSecAll(end)*binSize + sum(groupSynsAll))/(length(groupSynsAll)+binSize);
meanDistsGroupSecAll(end) = (meanDistsGroupSecAll(end)*binSize + sum(groupDistsAll))/(length(groupDistsAll)+binSize);



%Group by the cell body distance
[temp, sortedInds] = sort(neuronsPairsSum.cellBodyDist);
sortedInds = sortedInds(~isnan(temp));
groupsSections = nan(length(sortedInds),1);
counter = 0;
groupCounter = 1;
groupSyns = [];
groupSections = [];
groupDists =[];

pSynGroupDist = [];
nSynGroupDist = [];
meanDistsGroupDist = [];

for i=1:length(sortedInds)
    counter = counter + 1;
    
    groupSyns = [groupSyns neuronsPairsSum.nSynapses(sortedInds(i))];
    groupSections = [groupSections neuronsPairsSum.nCloseSections(sortedInds(i))];
    groupDists = [groupDists neuronsPairsSum.cellBodyDist(sortedInds(i))];
    if counter>=binSize
        counter = 0;
        groupCounter = groupCounter+1;
        
        pSynGroupDist = [pSynGroupDist sum(groupSyns>0)/length(groupSyns)];
        nSynGroupDist = [nSynGroupDist mean(groupSyns)];
        meanDistsGroupDist = [meanDistsGroupDist mean(groupDists)];
        
        groupSyns = [];
        groupSections = [];
        groupDists =[];
    end
end

pSynGroupDist(end) = (pSynGroupDist(end)*binSize+sum(groupSyns>0))/(length(groupSyns)+binSize);
nSynGroupDist(end) = (nSynGroupDist(end)*binSize + sum(groupSyns))/(length(groupSyns)+binSize);
meanDistsGroupDist(end) = (meanDistsGroupDist(end)*binSize + sum(groupDists))/(length(groupDists)+binSize);
fSize = 24

figure('units','normalized','outerposition',[0 0 1 1])
subplot(2,2,3)
scatter(nSectionsGroupSecAll.*dz,pSynGroupSecAll,dotSize,'filled','MarkerFaceColor',[1 0.8 0.8])
hold on
scatter(nSectionsGroupSec.*dz,pSynGroupSec,dotSize,'filled')
xlabel('Mean length on which neurites are proximal [\mum]')
ylabel('Fraction of pairs that share a synapse')



set(gca,'FontSize',fSize)
ylimP = ylim();
subplot(2,2,4)
scatter(nSectionsGroupSecAll.*dz,nSynGroupSecAll,dotSize,'filled','MarkerFaceColor',[1 0.8 0.8])
hold on
scatter(nSectionsGroupSec.*dz,nSynGroupSec,dotSize,'filled')

legend('All Pairs','Pairs with registered cell body postion')
xlabel('Mean length on which neurites are proximal [\mum]')
ylabel('Mean number of synapses')
set(gca,'FontSize',fSize)
ylimN = ylim();
subplot(2,2,1)
scatter(meanDistsGroupDist,pSynGroupDist,dotSize,'filled')
xlabel('Mean distance between cell bodies [\mum]')
ylabel('Probability to have a synapse')
set(gca,'FontSize',fSize)
ylim(ylimP)
subplot(2,2,2)
scatter(meanDistsGroupDist,nSynGroupDist,dotSize,'filled')
xlabel('Mean distance between cell bodies [\mum]')
ylabel('Mean number of synapses')
set(gca,'FontSize',fSize)
ylim(ylimN)


modelTable = neuronsPairsSum;
modelTable(isnan(modelTable.cellBodyDist),:) = [];

Y = nan(length(modelTable.nSynapses),1);
Y(modelTable.nSynapses==0) = 1;
Y(modelTable.nSynapses>0) = 2;
X = table((modelTable.cellBodyDist-nanmean(modelTable.cellBodyDist))/nanstd(modelTable.cellBodyDist), ...
    (modelTable.nCloseSections-nanmean(modelTable.nCloseSections))/nanstd(modelTable.nCloseSections),'VariableNames',{'Dist','CloseSectionsNum'});
[BLogi,~,statsLogi] = mnrfit(table2array(X),Y);
XY = table(X.Dist,X.CloseSectionsNum,Y,'VariableNames',{'Dist','CloseSectionsNum','Nsyns'});
mdlLin = fitlm(XY);

disp('after normalizing by mean and std:') 
disp(['Losgistic regression Vals for synapse probability: ' num2str(BLogi(2)) ' (cellBodyDist) and ' num2str(BLogi(3)) '(number of proximal sections)'])
disp(['Pvals are accordigly: ' num2str(statsLogi.p(2)) ' and ' num2str(statsLogi.p(3))])


disp('Linear regression Vals for synapses number: ')
mdlLin.Coefficients

corrNsynnCloseSections = corrcoef(modelTable.nCloseSections,modelTable.nSynapses);
corrNsynnCellDist = corrcoef(modelTable.cellBodyDist,modelTable.nSynapses);
corrSynProbCellDist = corrcoef(modelTable.cellBodyDist,modelTable.nSynapses);

disp(['correlation of nCloseSections and synapsesNum: is ' ...
    num2str(corrNsynnCloseSections(1,2))])

corrcoef(modelTable.nCloseSections,Y)
corrcoef(modelTable.cellBodyDist,double(modelTable.nSynapses>0))

%%

allProbsChem= nan(length(all_neurons),length(all_neurons));
dummyProbsChem= nan(length(all_neurons),length(all_neurons));

allProbsElec= nan(length(all_neurons),length(all_neurons));
dummyProbsElec= nan(length(all_neurons),length(all_neurons));
realMeans = nan(length(all_neurons),length(all_neurons));
randMeans = nan(length(all_neurons),length(all_neurons));

min_dL = dLThresh;
doPlot = 0;

for ii=1:length(all_neurons)
    disp(ii)
    thisN1 = all_neurons(ii);
    n1Syns = thisN1.all_synapses;
    n1Name= thisN1.NameStr;
    for jj = 1:length(all_neurons)
        thisN2 = all_neurons(jj);
        n2Name= thisN2.NameStr;
        synNumTo = sum(and(n1Syns.pre==n2Name,n1Syns.dierection=='toNeuron'));
        synNumFrom = sum(and(and(categorical(n1Syns.post1)==n2Name,categorical(n1Syns.pre)==n1Name),n1Syns.dierection=='fromNeuron'));
        if and(synNumFrom>=2,jj~=ii)
            [allProbsChem(ii,jj), dummyProbsChem(ii,jj), allProbsElec(ii,jj), dummyProbsElec(ii,jj),  realMeans(ii,jj), randMeans(ii,jj)]...
                = estimateDistancesProb(thisN1,thisN2,min_dL,doPlot,'from');
            
        end
        
    end 
end
save('pairWiseResultsFrom.mat')
%%
load('pairWiseResultsFrom.mat')
%%
%analyze the results (chem synapses)
figure('units','normalized','outerposition',[0 0 1 1])
allProbsChem(allProbsChem<0) = nan;
%dummyProbsChem(isnan(allProbsChem)) = nan;
histogram(allProbsChem(~isnan(allProbsChem)),20,'Normalization','probability')
hold on
%histogram(dummyProbsChem(~isnan(dummyProbsChem)),20,'Normalization','probability')
[N,edges] = histcounts(dummyProbsChem(~isnan(dummyProbsChem)),20);
plot(edges(1:end-1)+0.025,N/sum(N),'r--','LineWidth',2)
legend('Real synapses', 'Random Synapses')
ylabel('Fraction of pairs')
xlabel('Bootstrap P-value')
set(gca,'FontSize',30)
pbaspect([1 1 1])
correctedProbs = allProbsChem;
%correctedProbs(correctedProbs>1) = 1;
%correctedProbs = reshape(mafdr(correctedProbs(:),'BHFDR',true),size(correctedProbs));

topM = triu(correctedProbs)';
botM = tril(correctedProbs);
topM(topM==0) = nan;
botM(botM==0) = nan;

%num of pairs
sum(correctedProbs(:)<0.05);

%both ij and ji are significanly close
sum(and(topM(:)<0.05,botM(:)<0.05));

disp(['Out of ' num2str(sum(correctedProbs(:)<0.05)) ' significant pair, in ' num2str(sum(and(topM(:)<0.05,botM(:)<0.05)*2)...
) ' pairs are significant on both ends'])
[a b] = kstest2(allProbsChem(:),dummyProbsChem(:))
disp([num2str(sum(~isnan(correctedProbs(:)))) ' were tested'])
disp([num2str(sum(correctedProbs(:)<0.05)/sum(~isnan(correctedProbs(:)))) ' were significant'])
saveas(gcf,'C:\Users\Owner\Google Drive\nir-rotem_dentritic_computation\codes\newCodes\figs\SynapsesAggregationHistFrom.png')
saveas(gcf,'C:\Users\Owner\Google Drive\nir-rotem_dentritic_computation\codes\newCodes\figs\SynapsesAggregationHistFrom.emf')
savefig(gcf,'C:\Users\Owner\Google Drive\nir-rotem_dentritic_computation\codes\newCodes\figs\SynapsesAggregationHistFrom.fig')

realMeans(realMeans==0) = nan;
randMeans(randMeans==0)= nan;
disp(['mean distance of real syns is ' num2str(nanmean(realMeans(:)))])
disp(['mean distance of rand syns is ' num2str(nanmean(randMeans(:)))])

randMeans = tril(randMeans);
realMeans = tril(realMeans);
realMeans(realMeans==0) = nan;
randMeans(randMeans==0)= nan;

ranksum(realMeans(:),randMeans(:))
disp(['n= ' num2str(sum(~isnan(randMeans(:))))])
%%
%Both cond
allProbsChem= nan(length(all_neurons),length(all_neurons));
dummyProbsChem= nan(length(all_neurons),length(all_neurons));

allProbsElec= nan(length(all_neurons),length(all_neurons));
dummyProbsElec= nan(length(all_neurons),length(all_neurons));
realMeans = nan(length(all_neurons),length(all_neurons));
randMeans = nan(length(all_neurons),length(all_neurons));

min_dL = dLThresh;
doPlot = false;

for ii=1:length(all_neurons)
    thisN1 = all_neurons(ii);
    n1Syns = thisN1.all_synapses;
    n1Name= thisN1.NameStr;
    for jj = 1:length(all_neurons)
        thisN2 = all_neurons(jj);
        n2Name= thisN2.NameStr;
        synNumTo = sum(and(n1Syns.pre==n2Name,n1Syns.dierection=='toNeuron'));
        synNumFrom = sum(and(and(categorical(n1Syns.post1)==n2Name,categorical(n1Syns.pre)==n1Name),n1Syns.dierection=='fromNeuron'));
        if synNumFrom>=1 && synNumTo>=1 && jj~=ii 
            [allProbsChem(ii,jj), dummyProbsChem(ii,jj), allProbsElec(ii,jj), dummyProbsElec(ii,jj),  realMeans(ii,jj), randMeans(ii,jj)]...
                = estimateDistancesProb(thisN1,thisN2,min_dL,doPlot,'both');
            
        end
        
    end 
end
save('pairWiseResultsBoth.mat')
%%
load('pairWiseResultsBoth.mat')
%%
%analyze the results (chem synapses)
figure('units','normalized','outerposition',[0 0 1 1])
subplot(1,2,1)
allProbsChem(allProbsChem<0) = nan;
%dummyProbsChem(isnan(allProbsChem)) = nan;
histogram(allProbsChem(~isnan(allProbsChem)),20,'Normalization','probability')
hold on
%histogram(dummyProbsChem(~isnan(dummyProbsChem)),20,'Normalization','probability')
[N,edges] = histcounts(dummyProbsChem(~isnan(dummyProbsChem)),20);
plot(edges(1:end-1)+0.025,N/sum(N),'r--','LineWidth',2)
legend('Registered synapses', 'Random Synapses')
ylabel('Fraction of pairs')
xlabel('Bootstrap P-value')
set(gca,'FontSize',30)
%pbaspect([1 1 1])
correctedProbs = allProbsChem;
%correctedProbs(correctedProbs>1) = 1;
%correctedProbs = reshape(mafdr(correctedProbs(:),'BHFDR',true),size(correctedProbs));

topM = triu(correctedProbs)';
botM = tril(correctedProbs);
topM(topM==0) = nan;
botM(botM==0) = nan;

%num of pairs
sum(correctedProbs(:)<0.05);

%both ij and ji are significanly close
sum(and(topM(:)<0.05,botM(:)<0.05));

disp(['Out of ' num2str(sum(correctedProbs(:)<0.05)) ' significant pair, in ' num2str(sum(and(topM(:)<0.05,botM(:)<0.05)*2)...
) ' pairs are significant on both ends'])
[a b] = kstest2(allProbsChem(:),dummyProbsChem(:))
disp([num2str(sum(~isnan(correctedProbs(:)))) ' were tested'])
disp([num2str(sum(correctedProbs(:)<0.05)/sum(~isnan(correctedProbs(:)))) ' were significant'])

randMeans = tril(randMeans);
realMeans = tril(realMeans);
realMeans(realMeans==0) = nan;
randMeans(randMeans==0)= nan;

disp(['mean distance of real syns is ' num2str(nanmean(realMeans(:)))])
disp(['mean distance of rand syns is ' num2str(nanmean(randMeans(:)))])

ranksum(realMeans(:),randMeans(:))

disp(['n= ' num2str(sum(~isnan(randMeans(:))))])

%%


%%

function [row] = computePair(neuron1, neuron2, min_dL)

nerveRingZ = 20;
n1Name = categorical(cellstr(neuron1.NameStr));
n2Name = categorical(cellstr(neuron2.NameStr));

n1Locs = table(neuron1.rows.objName1, neuron1.rows.x1, neuron1.rows.y1, neuron1.rows.z1,'VariableNames',{'obj', 'x','y','z'});
n1Locs = [n1Locs ; table(neuron1.rows.objName2,neuron1.rows.x2, neuron1.rows.y2, neuron1.rows.z2,'VariableNames',{'obj','x','y','z'})];
n1Locs = unique(n1Locs);
n1Locs = n1Locs(n1Locs.z<nerveRingZ,:);


n2Locs = table(neuron2.rows.objName1, neuron2.rows.x1, neuron2.rows.y1, neuron2.rows.z1,'VariableNames',{'obj', 'x','y','z'});
n2Locs = [n2Locs ; table(neuron2.rows.objName2,neuron2.rows.x2, neuron2.rows.y2, neuron2.rows.z2,'VariableNames',{'obj','x','y','z'})];
n2Locs = unique(n2Locs);
n2Locs = n2Locs(n2Locs.z<nerveRingZ,:);

objsDists = pdist2(table2array(n1Locs(:,2:4)),table2array(n2Locs(:,2:4)));

minD12 = min(objsDists,[],2);

closeLocs1 = n1Locs(minD12<=min_dL,:);
%closeLocs1Idx = arrayfun(@(x)find(neuron1.uniqueObjs==x,1),closeLocs1.obj);


realSyns = neuron1.all_synapses;
realSynsInds1 = find(or(realSyns.pre==n2Name,and(realSyns.post1==n2Name,realSyns.pre==n1Name)));

n1X = neuron1.cellBodyX;
n1Y = neuron1.cellBodyY;
n1Z = neuron1.cellBodyZ;

n2X = neuron2.cellBodyX;
n2Y = neuron2.cellBodyY;
n2Z = neuron2.cellBodyZ;
if ~isnan(n1X) && ~isnan(n2X)
    cellBodyDist = sqrt((n1X-n2X)^2+(n1Y-n2Y)^2+(n1Z-n2Z)^2);
else
    cellBodyDist = nan;
end

%realSynsDistmat = neuron1.dist_map(realSynsInds1,realSynsInds1);
synnum = length(realSynsInds1);




%realSynsDistmatOthers = reshape(realSynsDistmat(~eye(size(realSynsDistmat))), size(realSynsDistmat, 2)-1, []);
%meanMin1Dists = mean(min(realSynsDistmatOthers,[],2));
% n1AdjMat = neuron1.adj_mat;
% n1Objs = neuron1.uniqueObjs;
% if synnum>0
%     for i=1:synnum
%         thisSynObj = realSyns(realSynsInds1(i),:);
%         thisSynDists = n1AdjMat(thisSynObj.adjMatInd,:);
row = table(n1Name,n2Name,synnum,height(closeLocs1),cellBodyDist,'VariableNames',{'neuron1', 'neuron2','nSynapses','nCloseSections','cellBodyDist'});
end


