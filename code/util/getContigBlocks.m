function [blockSizes, blockStarts, blockSymbols] = getContigBlocks( ids )
%  Given vector "ids" of discrete ids (integers)
%     find all contiguous blocks within it

blockSizes = [];
blockStarts = [];
blockSymbols = [];
aa = 1;
while aa <= length( ids )
   
    blockStarts(end+1) = aa;
    blockSymbols(end+1) = ids(aa);
    
    if aa == length( ids )
        % Reached the end of the input vector
        blockSizes(end+1) = 1;
        break;
    else
        % Search forwards
        didFindMismatch = 0;
        for bb = aa+1:length(ids)
            if ids(bb) ~= ids(aa);
                didFindMismatch = 1;
                bb = bb - 1;
                break;
            end
        end
        if didFindMismatch
            blockSizes(end+1) = bb - aa + 1;
        else
            blockSizes(end+1) = length(ids) - aa + 1;
            break;
        end
        
    end
        
    aa = bb+1;
    
end
    