function [blockSizes, blockStarts, blockSymbols] = getContigBlocks( ids )
%  Given a sequence "ids" of discrete positive integers
%     find all contiguous blocks within it
%USAGE:
%  [S, L, A] = getContigBlocks( [1 1 1 2 2 2 2 1 3 3] );
%  will return [B,L,S], where
%      S : [3 4 1 2] is vector of block sizes
%      L : [1 4 8 9] is vector of starting locations 
%      A : [1 2 1 3] is vector of block labels (original vals in "ids")

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
    
