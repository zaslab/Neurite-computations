%comparing chem and electrical synapses
load('Supplementary_File1all_neurons.mat')
neuronsNames = cellfun(@(x) categorical(cellstr(x)),(extractfield(all_neurons,'NameStr')));
totalDistsChem = [];
totalDistsChemNorm = [];
totalDistsElec = [];
totalDistsElecNorm = [];
totalPvalsRSE_C = [];
totalPvalsRS_L = [];
totalPvalsRS_R = [];

A_ec = [];
B_ec = [];

donePairsA =categorical(); 
donePairsB = categorical();
doneCount = 0;
doPlot = true;
for i = 1:length(all_neurons)
    thisN = all_neurons(i);
    thisNname = categorical(cellstr(thisN.NameStr));
    thisDistMat = thisN.dist_map;
    thisSynListAll = thisN.all_synapses;
    thisPartnersList = unique([thisSynListAll.post1 thisSynListAll.pre]);
    thisPartnersList = thisPartnersList(thisPartnersList~=thisNname);
    for j=1:length(thisPartnersList)
        thisPartner = thisPartnersList(j);
        thisSynList = thisSynListAll(thisSynListAll.pre==thisPartner | thisSynListAll.post1==thisPartner,:);
        
        
        electricalSynsInds = find(thisSynList.type=='electrical');
        uInds = unique(thisSynList.idx(thisSynList.type=='chemical'));
        [~, inds] = intersect(thisSynList.idx,uInds);
        synIndsChem = unique(inds);
        
        maxDist = max(thisDistMat(:));
        allChemElecDists = thisDistMat(electricalSynsInds,synIndsChem);
        allChemDistsTemplate = tril(ones(size(allChemElecDists)),-1);
        allChemElecDists= allChemElecDists(allChemDistsTemplate>0);
        allChemElecDists = allChemElecDists(allChemElecDists>0);
        
   
        allElecDists = thisDistMat(electricalSynsInds,electricalSynsInds);
        allElecDistsTemplate =  tril(ones(size(allElecDists)),-1);
        allElecDists = allElecDists(allElecDistsTemplate>0);
        allElecDists = allElecDists(allElecDists>0);
        
        allChemDists = thisDistMat(synIndsChem,synIndsChem);
        allChemDistsTemplate =  tril(ones(size(synIndsChem)),-1);
        allChemDists = allChemDists(allChemDistsTemplate>0);
        allChemDists = allChemDists(allChemDists>0);

        indA = find(thisPartner==donePairsA);
        indB  = find(thisNname==donePairsB);
        isNew = isempty(intersect(indA,indB));
        if ~isNew
            
            doneCount = doneCount+1;
        end
        
        if (length(allElecDists(:))>1) && (length(allChemElecDists(:))>1)...
            && isNew
        
            donePairsA = [donePairsA thisNname];
            donePairsB = [donePairsB thisPartner];
            thisPV = ranksum([allElecDists; allChemDists],allChemElecDists,'tail','left');
            totalPvalsRSE_C = [totalPvalsRSE_C; thisPV];
             A_ec = [A_ec; thisNname];
             B_ec = [B_ec; thisPartner];
                 
        end
    end
end

LPV = -2*sum(log(totalPvalsRSE_C));
totPval =1-chi2cdf(LPV,2*length(totalPvalsRSE_C))
totalPvalsRS_chem_gap = totalPvalsRSE_C;


%%
%comparing chem dierctions

totalDistsIn = [];
totalDistsInNorm = [];
totalDistsOut = [];
totalDistsOutNorm = [];
totalPvalsRS_InOut = [];
compNum = 213;
A = [];
B = [];
donePairsA =categorical(); 
donePairsB = categorical();
doPlot = false;
doneCount= 0;
for i = 1:length(all_neurons)
    thisN = all_neurons(i);
    thisNname = categorical(cellstr(thisN.NameStr));
    thisDistMat = thisN.dist_map;
    thisSynListAll = thisN.all_synapses;
    thisPartnersList = unique([thisSynListAll.post1 thisSynListAll.pre]);
    thisPartnersList = thisPartnersList(thisPartnersList~=thisNname);
   
    for j=1:length(thisPartnersList)
        thisPartner = thisPartnersList(j);
        thisSynList = thisSynListAll(thisSynListAll.pre==thisPartner | thisSynListAll.post1==thisPartner,:);
       
        
        uIndsIn = find(thisSynList.type=='chemical' & thisSynList.pre==thisPartner & thisSynList.dierection=='toNeuron');
        uIndsOut = find(thisSynList.type=='chemical' & thisSynList.post1==thisPartner & thisSynList.dierection=='fromNeuron');      
        maxDist = max(thisDistMat(:));
       
        allInInDists = thisDistMat([uIndsIn , uIndsIn]);
        allInDistsTemplate = tril(ones(size(allInInDists)),-1);
        allInInDists= allInInDists(allInDistsTemplate>0);
        allInInDists = allInInDists(allInInDists>0);
        
        allOutOutDists = thisDistMat(uIndsOut,uIndsOut);
        allOutDistsTemplate = tril(ones(size(allOutOutDists)),-1);
        allOutOutDists= allOutOutDists(allOutDistsTemplate>0);
        allOutOutDists = allOutOutDists(allOutOutDists>0);
        
        allInOutDists = thisDistMat(uIndsOut,uIndsIn);
        allInOutDists = allInOutDists(allInOutDists>0);
        
        indA = find(thisPartner==donePairsA);
        indB  = find(thisNname==donePairsB);
        isNew = isempty(intersect(indA,indB));
        if ~isNew
            
            doneCount = doneCount+1;
        end
        
        if (length(allInInDists(:))>1) && (length(allOutOutDists)>1) && isNew
            donePairsA = [donePairsA thisNname];
            donePairsB = [donePairsB thisPartner];
            
            totalDistsIn = [totalDistsIn ; allInInDists(:)];
            totalDistsOut = [totalDistsOut; allOutOutDists(:)];
            
            totalDistsInNorm = [totalDistsInNorm; totalDistsIn(:)/maxDist];
            totalDistsOutNorm = [totalDistsOutNorm; totalDistsOut(:)/maxDist];

            thisPV = ranksum([allOutOutDists(:);allInInDists(:)],allInOutDists,'tail','left');
            totalPvalsRS_InOut = [totalPvalsRS_InOut; thisPV];
          
        end
    end
end


LPV = -2*sum(log(totalPvalsRS_InOut));
PV_in_out =1-chi2cdf(LPV,2*length(totalPvalsRS_InOut))
disp(['analyzed ' num2str(length(totalPvalsRS_InOut)) ' pairs'])
%%
%plotting both histograms
figure('units','normalized','outerposition',[0 0 1 1])
subplot(1,2,1)
histogram(totalPvalsRS_chem_gap,20)
xlabel('p value')
ylabel('Number of pairs')
fSize = 30;
set(gca,'FontSize',fSize)
set(gca,'YTick',1:9)
title('gap  and chemical ')


subplot(1,2,2)
histogram(totalPvalsRS_InOut,20,'Normalization','probability')
xlabel('p value')
ylabel('Fraction of pairs')

set(gca,'FontSize',fSize)
title(' dierection')


