function [ Cycles ] = findCycle(all_neurons,isOnly)


if ~exist('isOnly','var')
    isOnly = false;
end

neuronsNum = length(all_neurons);
Cycles = table();
neuronsNames = categorical(extractfield(all_neurons,'NameStr'));

for ii = 1:neuronsNum
    thisA = all_neurons(ii);
    disp(ii);
    thisAouts = unique(thisA.all_synapses.post1(thisA.all_synapses.dierection=='fromNeuron'));
    thisAins = unique(thisA.all_synapses.pre(thisA.all_synapses.dierection=='toNeuron'));
    for jj = 1:length(thisAouts)
        thisB = all_neurons(thisAouts(jj)==neuronsNames);
        if ~isempty(thisB)
            thisBIns = thisB.all_synapses.pre(thisB.all_synapses.dierection=='toNeuron');
            thisBouts = categorical(unique(thisB.all_synapses.post1(thisB.all_synapses.dierection=='fromNeuron')));
            thisCs = thisBouts;
            if and(or(~isOnly, isempty(find(thisCs==thisA.NameStr,1))),~isempty(find(thisBIns==thisA.NameStr,1)))
                for kk = 1:length(thisCs)
                    
                    thisC = all_neurons(thisBouts(kk)==neuronsNames);
                    if ~isempty(thisC)
                        thisCname = thisCs(kk);
                        thisCIns = thisC.all_synapses.pre(thisC.all_synapses.dierection=='toNeuron');
                        thisCouts = thisC.all_synapses.post1(thisC.all_synapses.dierection=='fromNeuron');
                        if (~isempty(find(thisCname==thisAins, 1))) && (~isempty(find(thisCIns==thisB.NameStr,1)))...
                                && (~isempty(find(thisCouts==thisA.NameStr,1)))
                            if or(~isOnly,isempty(find([thisBIns; thisAouts]==thisCname,1)))
                                
                                thisCyc = struct();
                                thisCyc.A = categorical(cellstr(thisA.NameStr));
                                thisCyc.B = categorical(cellstr(thisB.NameStr));
                                thisCyc.C = categorical(cellstr(thisCname));
                                Cycles = [Cycles ; struct2table(thisCyc)];
                            end
                        end
                    end
                end
            end
        end
        
    end
end
