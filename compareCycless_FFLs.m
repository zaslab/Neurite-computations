clear
load('Supplementary_File1all_neurons.mat')
load('Supplementary_File2_tripletsTable.mat')
figsSaveFold = 'C:\Users\Owner\Google Drive\nir-rotem_dentritic_computation\codes\newCodes\figs\';


%%
%find triplets of FFL and Cycles
Cycles = findCycle(all_neurons,true);
FFLs = findFFL(all_neurons,true);
neuronsNames = categorical(extractfield(all_neurons,'NameStr'));


%%
%get Inds of C from table (not all exist due to the limit on the
%number of permutations)
cyclesBInds = [];

for i=1:height(Cycles)
        cyclesBInds = [cyclesBInds; find(tripletsTable.commonNeuron==Cycles.B(i) & ...
        tripletsTable.N1name==Cycles.A(i) & tripletsTable.synapseType1=='toNeuron'...
        & tripletsTable.N2name==Cycles.C(i) & tripletsTable.synapseType2=='fromNeuron')];
end
%%
%get Inds of FFLs from grand table (not all exist due to the limit on the
%number of permutations)
fflBInds = [];
for i=1:height(FFLs)
    fflBInds = [fflBInds; find(tripletsTable.commonNeuron==FFLs.B(i) & ...
        tripletsTable.N1name==FFLs.A(i) & tripletsTable.synapseType1=='toNeuron'...
        & tripletsTable.N2name==FFLs.C(i) & tripletsTable.synapseType2=='fromNeuron')];
end

%%
%Count the number of synapses in cycles and FFLs
cyclesSynsAC = [];
cyclesSynsAB = [];
cyclesSynsBC = [];

for i=1:height(Cycles)
    
    onASyns = find(tripletsTable.commonNeuron==Cycles.A(i) & ...
        tripletsTable.N1name==Cycles.C(i) & tripletsTable.synapseType1=='toNeuron'...
        & tripletsTable.N2name==Cycles.B(i) & tripletsTable.synapseType2=='fromNeuron');
    if ~isempty(onASyns)
        cyclesSynsAC(i) = tripletsTable.nSynapses1(onASyns);
        cyclesSynsAB(i) = tripletsTable.nSynapses2(onASyns);
    else
        cyclesSynsAC(i) = nan;
        cyclesSynsAB(i) = nan;
    end
    
    onCSyns = find(tripletsTable.commonNeuron==Cycles.C(i) & ...
        tripletsTable.N1name==Cycles.B(i) & tripletsTable.synapseType1=='toNeuron'...
        & tripletsTable.N2name==Cycles.A(i) & tripletsTable.synapseType2=='fromNeuron');
    if ~isempty(onCSyns)
        cyclesSynsBC(i) = tripletsTable.nSynapses1(onCSyns);
    else
        cyclesSynsBC(i) = nan;
    end
    
    
    
        
end


fflSynsAC = [];
fflSynsAB = [];
fflSynsBC = [];

for i=1:height(FFLs)
    
        onASyns = find(tripletsTable.commonNeuron==FFLs.A(i) & ...
        tripletsTable.N1name==FFLs.C(i) & tripletsTable.synapseType1=='fromNeuron'...
        & tripletsTable.N2name==FFLs.B(i) & tripletsTable.synapseType2=='fromNeuron');
    if ~isempty(onASyns)
        fflSynsAC(i) = tripletsTable.nSynapses1(onASyns);
        fflSynsAB(i) = tripletsTable.nSynapses2(onASyns);
    else
        fflSynsAC(i) = nan;
        fflSynsAB(i) = nan;
    end
    
    onCSyns = find(tripletsTable.commonNeuron==FFLs.C(i) & ...
        tripletsTable.N1name==FFLs.B(i) & tripletsTable.synapseType1=='toNeuron'...
        & tripletsTable.N2name==FFLs.A(i) & tripletsTable.synapseType2=='toNeuron');
    if ~isempty(onCSyns)
        fflSynsBC(i) = tripletsTable.nSynapses1(onCSyns);
    else
        fflSynsBC(i) = nan;
    end
    
    
    
        
end
disp(['for FFL, the mean number of synapses is ' num2str(nanmean(fflSynsAB)) ' (AB) '  num2str(nanmedian(fflSynsAC)) ' (AC) '...
    num2str(nanmean(fflSynsBC)) ' (BC)']);
disp(['total mean in FFLS is:' num2str(nanmean([fflSynsAB fflSynsAC fflSynsBC]))])
disp(['for cycles, the mean number of synapses is the same for all due to symmetry- ' num2str(nanmean(cyclesSynsAB))])

disp(['The mean number of synapses is different (PV = ' num2str(ranksum(cyclesSynsAB,[fflSynsAB fflSynsAC fflSynsBC])) ')']) 


%%
figure('units','normalized','outerposition',[0 0 1 1])

combinedFFL = [fflSynsAB fflSynsAC fflSynsBC];
combinedFFLH = combinedFFL;

CycleSynsH = cyclesSynsAB;
CycleSynsH = CycleSynsH(~isnan(CycleSynsH));

[barsCycles] = histc(CycleSynsH,[.5:9.5 1000])/length(CycleSynsH);

[barsFFL] = histc(combinedFFLH,[.5:9.5 1000])/length(combinedFFLH);
barsFFL = barsFFL(1:10);
barsCycles = barsCycles(1:10);

subplot(2,1,1)
bar([barsFFL; barsCycles]')
set(gca,'FontSize',24);
set(gca,'XTickLabel',{'1','2','3','4','5','6','7','8','9','>=10'})
ylabel('Fraction')
xlabel('# of synapses between neurons')
legend('FFLs','Cycles')

figName = 'Histogram_synapsesNum_FL_Cycles';
saveas(gcf,[figsSaveFold figName '.png'])
saveas(gcf,[figsSaveFold figName '.emf'])
savefig(gcf,[figsSaveFold figName '.fig'])



%%
%Find packed FFLs from all the FFLs list
minL = 4;
packedFFLs = [];
allFFLsSynInds = table();
FFLcounts = 0;
allPackedCountsFFL = nan(height(FFLs),1);
allSynCounts = nan(height(FFLs),1);
packedFFL_triplets = 0;
fflLargestDists = [];
fflMeanDistPerTriplet = [];
%go over all FFLs
for fflInd=1:height(FFLs)
    thisTripletDists = [];
    thisNeuronA = all_neurons(neuronsNames==cellstr(FFLs.A(fflInd)));
    thisNeuronB = all_neurons(neuronsNames==cellstr(FFLs.B(fflInd)));
    thisNeuronC = all_neurons(neuronsNames==cellstr(FFLs.C(fflInd)));
    %look on all A->B and A->C synapses on A
    [BConA,synsAB,synsAC] = getSubDistMat(thisNeuronA,thisNeuronB.NameStr,'fromNeuron',thisNeuronC.NameStr,'fromNeuron');
    %for each A->B syn, find the closest A->C syn.
    [minBCdists minBCInds] = min(BConA,[],2);
    thisFFLsynInds = table();
    packedCount = 0;
    thisFFL_triplet_packed = 0;
    %go over all the A->B syns
    for synInd = 1:height(synsAB)
        isPacked = true;
        thisLargestDist = minBCdists(synInd);
        %Check if the closest A->C syn is closer than minL
        if ~(minBCdists(synInd)<minL)
            isPacked = false;
        end
        %Follow the A->B synapse
        thisSynAB = synsAB(synInd,:);
        %and the A->C synapse
        thisSynAC =  synsAC(minBCInds(synInd),:);
        %on neuron B, look for all the A->B and B->C syns.
        [AConB,synsBA,synsBC] = getSubDistMat(thisNeuronB,thisNeuronA.NameStr,'toNeuron',thisNeuronC.NameStr,'fromNeuron');
        %Find the A->B synapse on neuron B by the synapse index.
        thisSynABonBInd = find(thisSynAB.idx==synsBA.idx,1);
        %Look for the closest B->C synapse
        [thisClosestCBsyn thisClosestCBsynInd] = min(AConB(thisSynABonBInd,:));
        thisLargestDist = max([thisLargestDist thisClosestCBsyn]);
        if isempty(thisClosestCBsyn)
            isPacked = false;
        else
            %is it closer than minL?
            if ~(thisClosestCBsyn<minL)
                isPacked = false;
            end
            thisBCsyn = synsBC(thisClosestCBsynInd,:);
            %Find all the B->C and A->C distances synapses on C
            [ABonC,synsCA,synsCB] = getSubDistMat(thisNeuronC,thisNeuronA.NameStr,'toNeuron',thisNeuronB.NameStr,'toNeuron');
            %Find the A->C and B->C synapses on C and find the disntace
            %between them
            thisBonCInd = find(thisBCsyn.idx==synsCB.idx,1);
            thisAonCInd = find(thisSynAC.idx==synsCA.idx,1);
            thisDistOnC = ABonC(thisAonCInd,thisBonCInd);
            %If the distance is also smaller than minL, we found a
            %packed FFL.
            thisLargestDist = max([thisLargestDist thisDistOnC]);
            if ~(thisDistOnC<minL)
                isPacked = false;
            end
            %Save the synapses and add a counter to the packed
            %counter
            if isPacked==true
                thisSyns = table(FFLs.A(fflInd),FFLs.B(fflInd),FFLs.C(fflInd),thisSynAB.idx,thisSynAC.idx,...
                    thisBCsyn.idx,fflInd,'VariableNames',{'A','B','C','AtoBInd','AtoCInd','BtoCInd','fflInd'});
                thisFFLsynInds = [thisFFLsynInds;thisSyns];
                packedCount = packedCount+1;
                thisFFL_triplet_packed = 1;
            end      
            %Count the total number or AB synapses for an estimate of the
            %number of FFL circuits (a triplet can have many circuits)
            FFLcounts = FFLcounts+1;
            fflLargestDists(end+1) = thisLargestDist;
            thisTripletDists(end+1) = thisLargestDist;
        end     
    end
    allFFLsSynInds = [allFFLsSynInds; thisFFLsynInds];
    fflMeanDistPerTriplet(fflInd) = nanmean(thisTripletDists);
    packedFFL_triplets = packedFFL_triplets+thisFFL_triplet_packed;
    allPackedCountsFFL(fflInd) = packedCount;
    allSynCounts(fflInd) = height(synsAB);
    
end
fractionPackedFFLs = height(allFFLsSynInds)/FFLcounts;
%%
packedCycles = [];
allCyclesSynInds = table();
FLcounts = 0;
allPackedCountsFL = nan(height(Cycles),1);
allSynCounts = nan(height(Cycles),1);
packedFL_triplets = 0;
flLargestDists = [];
flMeanDistPerTriplet = [];
for flInd=1:height(Cycles)
    thisTripletDists = [];
    thisNeuronA = all_neurons(neuronsNames==cellstr(Cycles.A(flInd)));
    thisNeuronB = all_neurons(neuronsNames==cellstr(Cycles.B(flInd)));
    thisNeuronC = all_neurons(neuronsNames==cellstr(Cycles.C(flInd)));
    [BConA,synsAB,synsAC] = getSubDistMat(thisNeuronA,thisNeuronB.NameStr,'fromNeuron',thisNeuronC.NameStr,'toNeuron');
    [minBCdists minBCInds] = min(BConA,[],2);
    thisCyclesynInds = table();
    packedCount = 0;
    thisFL_packed = 0;
    for synInd = 1:height(synsAB)
        isPacked = true;
        thisLargestDist = minBCdists(synInd);
        
        if ~(minBCdists(synInd)<minL)
            isPacked = false;
        end
        thisSynAB = synsAB(synInd,:);
        thisSynAC =  synsAC(minBCInds(synInd),:);
        [AConB,synsBA,synsBC] = getSubDistMat(thisNeuronB,thisNeuronA.NameStr,'toNeuron',thisNeuronC.NameStr,'fromNeuron');
        thisSynABonBInd = find(thisSynAB.idx==synsBA.idx,1);
        [thisClosestCBsyn thisClosestCBsynInd] = min(AConB(thisSynABonBInd,:));
        thisLargestDist = max([thisLargestDist thisClosestCBsyn]);
        if isempty(thisClosestCBsyn)
            thisClosestCBsyn = 10000;
            isPacked = false;
        else
            
            if  ~(thisClosestCBsyn<minL)
                isPacked = false;
            end
            thisBCsyn = synsBC(thisClosestCBsynInd,:);
            [ABonC,synsCA,synsCB] = getSubDistMat(thisNeuronC,thisNeuronA.NameStr,'fromNeuron',thisNeuronB.NameStr,'toNeuron');
            thisBonCInd = find(thisBCsyn.idx==synsCB.idx,1);
            thisAonCInd = find(thisSynAC.idx==synsCA.idx,1);
            thisDistOnC = ABonC(thisAonCInd,thisBonCInd);
            thisLargestDist = max([thisLargestDist thisDistOnC]);
            if ~(thisDistOnC<minL)
                isPacked = false;
            end
            
            if isPacked==true
                thisSyns = table(Cycles.A(flInd),Cycles.B(flInd),Cycles.C(flInd),thisSynAB.idx,thisSynAC.idx,...
                    thisBCsyn.idx,'VariableNames',{'A','B','C','AtoBInd','BtoCInd','CtoAInd'});
                thisCyclesynInds = [thisCyclesynInds;thisSyns];
                packedCount = packedCount+1;
                thisFL_packed = 1;
                
            end
            FLcounts = FLcounts+1;
            flLargestDists(end+1) = thisLargestDist;
            thisTripletDists(end+1) = thisLargestDist;
        end
        
    end
    allCyclesSynInds = [allCyclesSynInds; thisCyclesynInds];
    flMeanDistPerTriplet(flInd) = mean(thisTripletDists);
    packedFL_triplets = packedFL_triplets+thisFL_packed;
    allPackedCountsFL(flInd) = packedCount;
    allSynCounts(flInd) = height(synsAB);
end
fractionPackedCycles = height(allCyclesSynInds)/FLcounts;
%%

packedFractrionPval = proportionTest(height(allCyclesSynInds), FLcounts,height(allFFLsSynInds), FFLcounts);
disp(['From ' num2str(FLcounts) ' FL circuits ' num2str(height(allCyclesSynInds)) ' were packed (fraction of ' num2str(fractionPackedCycles) ')'])
disp(['From ' num2str(FFLcounts) ' FFL circuits ' num2str(height(allFFLsSynInds)) ' were packed (fraction of ' ...
    num2str(height(allFFLsSynInds)/FFLcounts) ')'])
disp(['two tailed proportion test is:' num2str(packedFractrionPval)])

fractionPackedFL_triplets = packedFL_triplets/height(Cycles);
fractionPackedFFL_triplets = packedFFL_triplets/height(FFLs);
disp(['From ' num2str(height(Cycles)) ' FL triplets ' num2str(packedFL_triplets) ' were packed (fraction of ' num2str(fractionPackedFL_triplets) ')'])
disp(['From ' num2str(height(FFLs)) ' FFL triplets ' num2str(packedFFL_triplets) ' were packed (fraction of ' num2str(fractionPackedFFL_triplets) ')'])
disp(['two tailed proportion test is:' num2str(proportionTest(packedFL_triplets ,height(Cycles),packedFFL_triplets,height(FFLs)))])


%%


figure('units','normalized','outerposition',[0 0 1 1])

maxBin = max([flMeanDistPerTriplet fflMeanDistPerTriplet]);
edges = [0:4:maxBin+1];
histogram(fflMeanDistPerTriplet,edges,'Normalization','probability')
hold on
histogram(flMeanDistPerTriplet,edges,'Normalization','probability')
legend('FFLs','Cycles')
xlabel('Average circuit length [\mum]')
ylabel('Fraction')
set(gca,'FontSize',24)
figName = 'packed_Circiuits_histogram_mean_per_triplet';
pbaspect([1 1 1])

saveas(gcf,[figsSaveFold figName '.png'])
saveas(gcf,[figsSaveFold figName '.emf'])
savefig(gcf,[figsSaveFold figName '.fig'])

ranksum(fflMeanDistPerTriplet,flMeanDistPerTriplet)


