load('Supplementary_File1all_neurons.mat')
%%

allDistsOnNeuron = [];
allDistsOclid = [];
for i=1:length(all_neurons)
    thisSyns = all_neurons(i).synapsesObjs;
    points = [thisSyns.x1 thisSyns.y1 thisSyns.z1];
    D = pdist(points);
    
    thisDists = all_neurons(i).dist_map;
    thisDists(thisDists==0)=-1;
    thisDists = tril(thisDists,-1);
    thisDists(thisDists==0)=nan;
    thisDists(thisDists==-1) = 0;
    allDistsOnNeuron = [allDistsOnNeuron; thisDists(~isnan(thisDists))];
    allDistsOclid = [allDistsOclid; D'];
   
end
allDistsOclidOrig = allDistsOclid;
allDistsOnNeuronOrig = allDistsOnNeuron;
%%
figure('units','normalized','outerposition',[0 0 1 1])

%allDistsOnNeuronCap = allDistsOnNeuron;
%allDistsOnNeuronCap(allDistsOnNeuronCap>30) = 30;
maxOc = 20;
maxOnNeuron=28;
allDistsOclid = allDistsOclidOrig;
allDistsOnNeuron = allDistsOnNeuronOrig;
allDistsOclid(allDistsOclid>maxOc) = maxOc;
allDistsOnNeuron(allDistsOnNeuron>maxOnNeuron) = maxOnNeuron;

histogram2(allDistsOclid,allDistsOnNeuron,[0:2:maxOc],[0:2:maxOnNeuron],'DisplayStyle','tile','Normalization','probability')
hcb = colorbar;
title(hcb,'Fraction')
set(gca,'ColorScale','log')

set(gca,'YTick',[0 4 8 12 16 20 24 28])
set(gca,'YTickLabel',{'0', '4', '8', '12', '16' ,'20','24', '28+'})
set(gca,'XTick',[0 4 8 12 16 20])
set(gca,'XTickLabel',{'0', '4', '8', '12', '16', '20+'})

xlabel('Euclidian distance [um]')
ylabel('Distance along neurite [um]')
set(gca,'FontSize',28)
hold on
plot([0 20], [0 20],'r','LineWidth',4)
%xlim([0 12])
daspect([1 1 1])
%%

%%

cols = jet(length(all_neurons));
figure('units','normalized','outerposition',[0 0 1 1])
hold on
zLim = 14;
nTypes = categorical(extractfield(all_neurons,'type'));
for i=1:length(all_neurons)
    thisRows = all_neurons(i).rows(all_neurons(i).rows.z1<zLim,:);
    if nTypes(i)=='interneuron'
        thisCol = [117,112,179]/255;
    elseif nTypes(i)=='sensory'
        thisCol = [27,158,119]/255;
        
    else
        thisCol = [217,95,2]/255;
        
    end
    
    for ii=1:height(thisRows)
        
    %plot3([thisRows.z1(ii) thisRows.z2(ii)], -1*[thisRows.y1(ii) thisRows.y2(ii)]...
    %    ,[thisRows.x1(ii) thisRows.x2(ii)],'Color',cols(i,:));  
    if (thisRows.cellbody1(ii)==1 || thisRows.cellbody2(ii)==1) && abs(thisRows.z1(ii) - thisRows.z2(ii))<0.5
        lw = 5;
        lineCol = min([[0.3 0.3 0.3]+thisCol; [1 1 1]]);
        
    else
        lw = 1;
        lineCol = thisCol;
    end
        
    plot3([thisRows.z1(ii) thisRows.z2(ii)], -1*[thisRows.y1(ii) thisRows.y2(ii)]...
        ,[thisRows.x1(ii) thisRows.x2(ii)],'Color',lineCol,'LineWidth',lw);  
    end
    
   
end

plot3(allNSynapses.z(allNSynapses.z<zLim & allNSynapses.type=='chemical'),-1*allNSynapses.y(allNSynapses.z<zLim & allNSynapses.type=='chemical'),allNSynapses.x(allNSynapses.z<zLim & allNSynapses.type=='chemical')...
    ,'.','color','black','MarkerSize',10)
xlim([0 zLim])
ylim([-8 0])
daspect([1 1 1])

%%
figure('units','normalized','outerposition',[0 0 1 1])
hold on
zLim = 14;
pairNames = categorical(unique(extractfield(all_neurons,'nameStart')));
colsPairs = jet(length(pairNames));
for i=1:length(all_neurons)
    thisRows = all_neurons(i).rows(all_neurons(i).rows.z1<zLim,:);
    thisCol = colsPairs(find(pairNames==all_neurons(i).nameStart),:);
    for ii=1:height(thisRows)
    plot3([thisRows.x1(ii) thisRows.x2(ii)], -1*[thisRows.y1(ii) thisRows.y2(ii)]...
        ,[thisRows.z1(ii) thisRows.z2(ii)],'Color',thisCol);  
    end
end

plot3(allNSynapses.x(allNSynapses.z<zLim),-1*allNSynapses.y(allNSynapses.z<zLim),allNSynapses.z(allNSynapses.z<zLim)...
    ,'black.','MarkerSize',8)

daspect([1 1 1])
ylim([-8 0])
saveas(gcf,[folderName 'allNsFront.png'])
savefig([folderName 'allNsFront.fig'])


