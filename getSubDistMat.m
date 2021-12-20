function [subDistMat,rowsA,colsB] = getSubDistMat(neuron,nameA,dierectionA,nameB,dierectionB)

if strcmp(dierectionA,'fromNeuron')
    relevantIndsA = find(and(neuron.all_synapses.post1 == nameA, neuron.all_synapses.dierection == dierectionA));
elseif strcmp(dierectionA, 'toNeuron')
    relevantIndsA = find(and(neuron.all_synapses.pre == nameA, neuron.all_synapses.dierection == dierectionA));
end

if strcmp(dierectionB, 'fromNeuron')
    relevantIndsB = find(and(neuron.all_synapses.post1 == nameB, neuron.all_synapses.dierection == dierectionB));
elseif strcmp(dierectionB, 'toNeuron')
    relevantIndsB = find(and(neuron.all_synapses.pre == nameB, neuron.all_synapses.dierection == dierectionB));
end
subDistMat = neuron.dist_map(relevantIndsA,relevantIndsB);
rowsA = neuron.all_synapses(relevantIndsA,:);
colsB = neuron.all_synapses(relevantIndsB,:);

end

