function [equivalence, equivalenceSplits] = getEquivalenceMatches(prediction, response)
nMethods = size(prediction, 2);
pairXOR = zeros(size(prediction, 1));

for method = 2:nMethods
    if method == 2
        pairXOR = xor(prediction(:, method-1), prediction(:, method)); %1st pair
    else
        pairXOR = xor(pairXOR, prediction(:, method)); %Keeps updating with subsequent pairs
    end
end
equivalence = ~pairXOR;

equivalenceSplits = zeros(2, 2);
for sample = 1:size(equivalence, 1)
    if equivalence(sample) %Match
        if response(sample) %Time Cell
            equivalenceSplits(2, 1) = equivalenceSplits(2, 1) + 1;
        else %Other Cell
            equivalenceSplits(2, 2) = equivalenceSplits(2, 2) + 1;
        end
    else %No Match
        if response(sample) %Time Cell
            equivalenceSplits(1, 1) = equivalenceSplits(1, 1) + 1;
        else %Other Cell
            equivalenceSplits(1, 2) = equivalenceSplits(1, 2) + 1;
        end
    end
end