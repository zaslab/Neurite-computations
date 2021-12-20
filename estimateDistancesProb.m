function [shuffleProbChem, dummyProbChem, shuffleProbElec,dummyProbElec, realMean, randMean] = estimateDistancesProb(neuron1,neuron2,min_dL,doPlot,D1)

if ~exist('doPlot')
    doPlot = 0;
end

shuffleProbChem = nan;
dummyProbChem = nan;
dummyProbElec = nan;
realMean = nan;
randMean = nan;

nerveRingZ = 20;
repeats = 200;
minRepeats = 20;
n1Name = neuron1.NameStr;
n2Name = neuron2.NameStr;

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
farLocs1Inds = find(minD12>min_dL); % used for plot
closeLocs1Idx = arrayfun(@(x)find(neuron1.uniqueObjs==x,1),closeLocs1.obj);


realSyns = neuron1.all_synapses;

fromIndsChem = find(and(categorical(realSyns.post1)==n2Name, realSyns.dierection == 'fromNeuron'));
toIndsChem = find(and(categorical(realSyns.pre)==n2Name, realSyns.dierection == 'toNeuron'));
if strcmp(D1,'from')
    allIndsChem = fromIndsChem;
elseif strcmp(D1,'to')
     allIndsChem = toIndsChem;

elseif strcmp(D1,'both')
     allIndsChem = unique([toIndsChem; fromIndsChem]);
     nFrom = length(fromIndsChem);
     nTo = length(toIndsChem);
end

fromIndsElec = find(and(categorical(realSyns.post1)==n2Name, realSyns.dierection == 'electrical'));
toIndsElec = find(and(categorical(realSyns.pre)==n2Name, realSyns.dierection == 'electrical'));

allIndsElec = unique([fromIndsElec; toIndsElec]);


allAdjMatIndsChem = realSyns.adjMatInd(allIndsChem);
allAdjMatIndsElec = realSyns.adjMatInd(allIndsElec);



adjMat = neuron1.adj_mat;
%First - compute the distance between all the close objects
closeLocs1Idx = unique([closeLocs1Idx; allAdjMatIndsChem; allAdjMatIndsElec]);

closeObjsDistMat = nan(height(closeLocs1),height(closeLocs1));
for ii = 1:length(closeLocs1Idx)
    thisObj1Ind = closeLocs1Idx(ii);
    for jj = 1:length(closeLocs1Idx)
        thisObj2Ind = closeLocs1Idx(jj);
        closeObjsDistMat(ii,jj) = graphshortestpath(adjMat,thisObj1Ind,thisObj2Ind,'Directed',false);
    end
end
if strcmp(D1,'both')
    realSynsDistmatChem = neuron1.dist_map(toIndsChem,fromIndsChem);
    meanDistsChem = mean(realSynsDistmatChem(:));
    synNumChem = length(allIndsChem);
else   
    realSynsDistmatChem = neuron1.dist_map(allIndsChem,allIndsChem);
    synNumChem = length(allIndsChem);
    realSynsDistmatnoDiagChem = reshape(realSynsDistmatChem(~eye(size(realSynsDistmatChem))), size(realSynsDistmatChem, 2)-1, [])';
    meanDistsChem = mean(realSynsDistmatnoDiagChem(:));

    
end

%Shuffle all Chemical Syns, measure all Chemical synapses
useAllperm = false;
if length(closeLocs1Idx)<= synNumChem
    shuffleProbChem = -1.1;
    
elseif nchoosek(length(closeLocs1Idx),synNumChem)<repeats && nchoosek(length(closeLocs1Idx),synNumChem)>minRepeats
    allPerms = nchoosek(1:length(closeLocs1Idx),synNumChem);
    useAllperm = true;
    repeats = size(allPerms,1);
    
elseif nchoosek(length(closeLocs1Idx),synNumChem)<minRepeats
    shuffleProbChem = -1.2;
end

allMeanDists1 = nan(repeats,1);
run = false;
if isnan(shuffleProbChem)
    run = true;
    for jj = 1:repeats
        %thisSynapsesInds = closeLocs1Idx(allsynapsesInds(:,jj));
        if useAllperm
            thisSynapsesInds = allPerms(jj,1:synNumChem);
        else
            thisSynapsesInds = randperm(length(closeLocs1Idx),synNumChem);
        end
        if strcmp(D1,'both')
            thisDistMat = closeObjsDistMat(thisSynapsesInds(1:nFrom),thisSynapsesInds(nFrom+1:end));
            allMeanDists1(jj) =  mean(thisDistMat(:));
        else
            
            thisDistMat = closeObjsDistMat(thisSynapsesInds,thisSynapsesInds);
            
            thisSynsDistmatnoDiag = reshape(thisDistMat(~eye(size(thisDistMat))), size(thisDistMat, 2)-1, [])';
            allMeanDists1(jj) =  mean(thisSynsDistmatnoDiag(:));

        end
    end
end

if run==true
    shuffleProbChem = (sum(allMeanDists1<meanDistsChem)+1)/(sum(~isnan(allMeanDists1))+1);
    dummyProbChem = (sum(allMeanDists1<allMeanDists1(randi(length(allMeanDists1))))+1)/(sum(~isnan(allMeanDists1))+1);
    realMean = meanDistsChem;
    randMean = mean(allMeanDists1);
   
else
    dummyProbChem = nan;
    shuffleProbChem = nan;
end

realSynsDistmatElec = neuron1.dist_map(allIndsElec,allIndsElec);
synNumElec = length(allIndsElec);
if synNumElec>=2
    shuffleProbElec = -2;
    realSynsDistmatnoDiagElec = reshape(realSynsDistmatElec(~eye(size(realSynsDistmatElec))), size(realSynsDistmatElec, 2)-1, [])';
    meanDistsElec= mean(realSynsDistmatnoDiagElec(:));

    useAllperm = false;
    if length(closeLocs1Idx)<= synNumElec
        shuffleProbElec = -1.1;
        
    elseif nchoosek(length(closeLocs1Idx),synNumElec)<repeats && nchoosek(length(closeLocs1Idx),synNumElec)>minRepeats
        allPerms = nchoosek(1:length(closeLocs1Idx),synNumElec);
        useAllperm = true;
        repeats = size(allPerms,1);
        
    elseif nchoosek(length(closeLocs1Idx),synNumElec)<minRepeats
        shuffleProbElec=-1.2;
    end
    
    allMeanDists1 = nan(repeats,1);
    
    if shuffleProbElec == -2
        for jj = 1:repeats
            %thisSynapsesInds = closeLocs1Idx(allsynapsesInds(:,jj));
            if useAllperm
                thisSynapsesInds = allPerms(jj,1:synNumElec);
            else
                thisSynapsesInds = randperm(length(closeLocs1Idx),synNumElec);
            end
            thisDistMat = closeObjsDistMat(thisSynapsesInds,thisSynapsesInds);
            
            thisSynsDistmatnoDiag = reshape(thisDistMat(~eye(size(thisDistMat))), size(thisDistMat, 2)-1, [])';

            allMeanDists1(jj) =  mean(thisSynsDistmatnoDiag(:));

            
        end
    end
    
    
    shuffleProbElec = (sum(allMeanDists1<meanDistsElec)+1)/(repeats+1);
    dummyProbElec = (sum(allMeanDists1<=allMeanDists1(randi(length(allMeanDists1)))))/(repeats);

    
    
else
    shuffleProbElec = nan;
end

close all

if and(shuffleProbChem<0.05,doPlot==1)
    pallet = jet(100);
   
    %normD = (minD12-min(minD12))/(max(minD12)-min(minD12));

    normD = minD12;
    colorRange = min_dL-min(minD12);
 
    normD = ((normD-min(minD12))/colorRange);
    normD = round(normD*99)+1;
    normD(normD>100)=100;
    
    fromSyns = realSyns(fromIndsChem,:);
    toSyns = realSyns(toIndsChem,:);
    figure('units','normalized','outerposition',[0 0 1 1])
    hold on
    %scatter3(n1Locs.z, -1*n1Locs.y,n1Locs.x,'black')
    scatter3(n1Locs.z, -1*n1Locs.y,n1Locs.x,60,pallet(normD,:),'filled')
    %scatter3(n2Locs.z, -1*n2Locs.y,n2Locs.x,'cyan')
    %scatter3(closeLocs1.z,-1*closeLocs1.y,closeLocs1.x,50,'MarkerEdgeColor','r','LineWidth',2)
    scatter3(fromSyns.z,-1*fromSyns.y,fromSyns.x,80,[0.4 0.4 0.4],'LineWidth',4)
    scatter3(toSyns.z,-1*toSyns.y,toSyns.x,80,[0.8 0.8 0.8],'LineWidth',4)
    daspect([1 1 1])
    caxis([min(minD12)*1000 min_dL*1000])
    
    hcb = colorbar;
    title(hcb,'Distance [nm]')
    set(gca,'FontSize',15)
    colormap jet
    title([neuron1.NameStr ' ' neuron2.NameStr ' p =' num2str(shuffleProbChem)...
        'Snyanpses num is ' num2str(synNumChem)]);
    legend('Non close Sections', [neuron1.NameStr '->' neuron2.NameStr]...
        ,[neuron2.NameStr '->' neuron1.NameStr])

end

if and(shuffleProbChem<0.05,doPlot==2)
    %normD = (minD12-min(minD12))/(max(minD12)-min(minD12));
    scatterSize = 60;
    colN1 = [240;5 ; 5]'/255;
    colN2 = [5; 40; 220]'/255;
    colS1 = [153; 12; 162]'/255;
    colS2 = [32; 147; 24]'/255;
    
    n1LocsPlot = n1Locs(farLocs1Inds,:);
   
    fromSyns = realSyns(fromIndsChem,:);
    toSyns = realSyns(toIndsChem,:);
    figure('units','normalized','outerposition',[0 0 1 1])
    hold on
    scatter3(n1LocsPlot.z, -1*n1LocsPlot.y,n1LocsPlot.x,scatterSize,'MarkerEdgeColor',colN1,'LineWidth',1)
    scatter3(n2Locs.z, -1*n2Locs.y,n2Locs.x,scatterSize,'MarkerEdgeColor',colN2,'LineWidth',1)
    scatter3(closeLocs1.z,-1*closeLocs1.y,closeLocs1.x,scatterSize,'filled','MarkerFaceColor','red','MarkerEdgeColor','black')
    scatter3(fromSyns.z,-1*fromSyns.y,fromSyns.x,scatterSize,'filled','MarkerFaceColor',colS1)
    scatter3(toSyns.z,-1*toSyns.y,toSyns.x,scatterSize,'filled','MarkerFaceColor',colS2)
    daspect([1 1 1])
  
    
    set(gca,'FontSize',15)
    title([neuron1.NameStr ' ' neuron2.NameStr ' p =' num2str(shuffleProbChem)...
        'Snyanpses num is ' num2str(synNumChem)]);
    legend(neuron1.NameStr,neuron2.NameStr,'Close sections', [neuron1.NameStr '->' neuron2.NameStr]...
        ,[neuron2.NameStr '->' neuron1.NameStr])
    
end
