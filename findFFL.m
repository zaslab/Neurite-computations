function [ FFLs ] = findFFL(all_neurons,isOnly)
%This fuction receives the structured dataset of all the neurons and
%returns a table of all the Feed Forward loops within the dataset.
%notation: the FFL is composed of A,B,C neurons where A->B->C and A->C.
%isOnly - search for circuits that only have the connections of an FFL (if
%true) or at least the connections of an FFL (if false)
if ~exist('isOnly','var')
    isOnly = true;
end
neuronsNum = length(all_neurons);
%The final talble
FFLs = table();
neuronsNames = categorical(extractfield(all_neurons,'NameStr'));

%Go over all the neurons
for ii = 1:neuronsNum
    
    thisA = all_neurons(ii);
    disp(ii);
    %For each A neuron save the Z->A and A->Z partenrs on A.
    
    thisAouts = unique(thisA.all_synapses.post1(thisA.all_synapses.dierection=='fromNeuron'));
    thisAins = unique(thisA.all_synapses.pre(thisA.all_synapses.dierection=='toNeuron'));
    
    %Go over all the postsynaptic partners
    for jj = 1:length(thisAouts)
        
        %Find the record in the dateset for each posysynaptic neuron
        thisB = all_neurons(thisAouts(jj)==neuronsNames);
        
        if ~isempty(thisB)

            %Get all Z->B partners on B
            BInSyns = categorical(unique(thisB.all_synapses.pre(thisB.all_synapses.dierection=='toNeuron')));
            
            %Make sure the A->B connection is also marked on the B neuron
            %and not only on the A neuron.
            if ~isempty(find(BInSyns==thisA.NameStr,1))
                %find all possilble C neurons (which are the B->Z partners
                %on B)
                thisCs = categorical(unique(thisB.all_synapses.post1(thisB.all_synapses.dierection=='fromNeuron')));
                thisBouts = thisCs;
                %For the case of the only FFL circuit, make sure there is
                %no B->A conncetion
                if or(isempty(find(thisBouts==thisA.NameStr, 1)),~isOnly)
                    %Go over all possible C partners
                    for kk = 1:length(thisCs)
                        thisC = thisCs(kk);
                        
                        thisCN = all_neurons(thisCs(kk)==neuronsNames);
                        if ~isempty(thisCN)
                            %For each possible C neuron in the Dataset,
                            %save the Z->C neurons on C.
                            thisCIns = categorical(unique(thisCN.all_synapses.pre(thisCN.all_synapses.dierection=='toNeuron')));
                            
                            %make sure there is a B->C synapse and a A->C
                            %synapse on C.
                            if and(~isempty(find(thisB.NameStr==thisCIns, 1)),...
                                    ~isempty(find(thisA.NameStr==thisCIns,1)))
                                %For the case of the only FFL circuit, make
                                %sure there is no C->A and C->B synapses.
                                if or(isempty(find(thisAins==thisC, 1)) && isempty(find(BInSyns==thisC,1)),~isOnly)
                                    %make sure there is also a A->C
                                    %connection on A.
                                    if ~isempty(find(thisAouts==thisC,1))
                                        thisFFL = struct();
                                        thisFFL.A = categorical(cellstr(thisA.NameStr));
                                        thisFFL.B = categorical(cellstr(thisB.NameStr));
                                        thisFFL.C = thisC;
                                        FFLs = [FFLs ; struct2table(thisFFL)];
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end
