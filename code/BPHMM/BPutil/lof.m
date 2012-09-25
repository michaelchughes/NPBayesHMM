function [F_sorted sort_ind] = lof(F, recurN, doVerbose)
% Map binary matrices to Left-Ordered Form (lof)
%   by ordering columns in descending order from left-to-right
%      according to magnitude of binary number expressed by each column.
% Empty columns will be *removed*, so size(F_sorted,2) <= size(F,2)
% Assumes first row represents the *most* significant bit.
% Sorting performed using Matlab's built-in  bin2dec function,
%   which we use to map each column to a unique decimal integer.                                                   
% bin2dec.m can only handle binary strings of length ~50 or less (R2010a)
%     so we must perform recursive sorting if F has more than 50 rows.
% EXAMPLE
%   F = [0 0 0 1 1 0; 1 0 1 0 1 0; 0 1 1 1  0 0 ; 0 1 0 0 0 0];
%   F_sorted = lof( F );
%       F                      F_sorted
%   0 0 0 1 1 0       lof     1 1 0 0 0    (Note empty col. was removed)
%   1 0 1 0 1 0   ------- >   1 0 1 1 0
%   0 1 1 1 0 0               0 1 1 0 1
%   0 1 0 0 0 0               0 0 0 0 1
% COMMENTS
%  Most users should happily omit all input args but the first one.
%  The others are only for students that seek a deeper understanding of
%    different possible algorithms for obtaining left-ordered form
%  In practice, the default options should always yield the fastest code.
% INPUT (* indicates optional input which can be omitted )
%   F          := binary matrix (entries are either zero or one )
%   *recurN    := integer indicator for which recursion type to perform
%                           1 : sort first row, then recurse on F(2:end,:)
%                           2 : sort first 2 rows,   recurse on F(3:end,:)
%                (default) 50 : sort first 50 rows,  recurse on F(50:end,:)
%   *doVerbose := optional flag to print debugging info.
%                (default) 0  : no progress messages printed to stdout
% OUTPUT
%   F_sorted  :=  matrix F in Left-Ordered Form
%   sort_ind  :=  vector gives indices for transformation of F to F_sorted
%                      F_sorted = F( sort_ind, : )
% REFERENCES
%   Griffiths and Ghahramani. 
%     "Infinite Latent Feature Models and the Indian Buffet Process"
%   In particular: Fig. 2 and section 2.2 "Equivalence classes"

if ~exist( 'doVerbose', 'var' )
    doVerbose = 0;
end

if ~exist( 'recurN', 'var' )
    recurN = 50;
end

% Sort columns according to Left-Ordered Form
if doVerbose
    initBufStr = ' ';
else
    initBufStr = '';
end

switch recurN
    case 1
        sort_ind = recursiveLoF( F , initBufStr);
    case 2
        sort_ind = recursiveLoF2( F , initBufStr);
    case 50
        sort_ind = recursiveLoF50( F , initBufStr);
end

F_sorted = F(:,sort_ind);

% Remove Empty Columns
posColIDs = sum(F_sorted,1) > 0;
F_sorted = F_sorted(:, posColIDs);
sort_ind = sort_ind( posColIDs );

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [sort_ind] = recursiveLoF50( F, bufStr )
    if ~isempty( bufStr )
      fprintf( '%s Reforming matrix of size %d x %d\n', bufStr, size(F,1), size(F,2) );
      bufStr = [bufStr ' '];
    end
    
    MM = 50;
    
    if size( F , 2) == 0
        sort_ind = [];
    elseif size(F, 2) == 1
        sort_ind = 1;
    elseif size( F, 1 ) <= MM
        sort_ind = sortColsByBin2Dec( F );
    else
       
        % Sort the first fifty rows (which is max allowed by bin2dec)
        [sort_ind, sort_vals] = sortColsByBin2Dec( F( 1:MM , : )  );
        
        % Figure out where duplicates exist, and sort remaining rows within
        unique_vals = unique( sort_vals );
        dupHist = histc( sort_vals, unique_vals );
        
        dupIDs = find( dupHist > 1 );
        
        for dupID = dupIDs
        
            curIDs = find( sort_vals == unique_vals(dupID) );
            
            recurF = F( MM+1:end,  sort_ind( curIDs )  );
            
            subSortIDs = recursiveLoF50(  recurF,  bufStr );
            
            sort_ind( curIDs ) = sort_ind( curIDs(subSortIDs) );
            
        end
        
    end
    
    % make sure sort_ind has all unique entries
    assert( length( sort_ind) == length( unique(sort_ind) ), 'ERROR: Left-Ordered Form col. swap will fail.' );
    
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [sort_ind] = recursiveLoF2( F, bufStr )
    if ~isempty( bufStr )
      fprintf( '%s Reforming matrix of size %d x %d\n', bufStr, size(F,1), size(F,2) );
      bufStr = [bufStr ' '];
    end
    
    if size( F , 2) == 0
        sort_ind = [];
    elseif size(F, 2) == 1
        sort_ind = 1;
    elseif size( F, 1 ) < 50
        sort_ind = sortColsByBin2Dec( F );
    else
       ind_1 = find( F( 1, : ) );
       ind_11 = find( F(2, ind_1) );
       ind_10 = setdiff( 1:length( ind_1 ), ind_11 );
       sort_11 = recursiveLoF( F( 3:end, ind_1( ind_11 ) ), bufStr );
       sort_10 = recursiveLoF( F( 3:end, ind_1( ind_10 ) ), bufStr );

       
       ind_0 = setdiff( 1:size(F,2), ind_1 );       
       ind_01 = find( F(2, ind_0) );
       ind_00 = setdiff( 1:length( ind_0 ), ind_01 );
       sort_01 = recursiveLoF( F( 3:end, ind_0( ind_01 )  ), bufStr );
       sort_00 = recursiveLoF( F( 3:end, ind_0( ind_00 )  ), bufStr );
       
       ind_1 = ind_1( [ind_11( sort_11 )  ind_10( sort_10 ) ]   );
       ind_0 = ind_0( [ind_01( sort_01 )  ind_00( sort_00 )]   );

       sort_ind = [ ind_1  ind_0 ];
    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [sort_ind] = recursiveLoF( F, bufStr )
    if ~isempty( bufStr )
      fprintf( '%s Reforming matrix of size %d x %d\n', bufStr, size(F,1), size(F,2) );
      bufStr = [bufStr ' '];
    end
    
    if size( F , 2) == 0
        sort_ind = [];
    elseif size(F, 2) == 1
        sort_ind = 1;
    elseif size( F, 1 ) < 50
        sort_ind = sortColsByBin2Dec( F );
    else
       ind_1 = find( F( 1, : ) );
       sort_1 = recursiveLoF( F( 2:end, ind_1 ), bufStr );

       ind_0 = setdiff( 1:size(F,2), ind_1 );       
       sort_0 = recursiveLoF( F( 2:end, ind_0 ), bufStr );
       
       sort_ind = [ ind_1( sort_1 )  ind_0( sort_0 ) ];
    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [sort_ind, sort_vals] = sortColsByBin2Dec( F  )

[N,Kz] = size(F);
%val = zeros(1,Kz);

% Rolled my own bin2dec... should be much faster!
val = sum( bsxfun(@times, F, pow2(N-1:-1:0)' ), 1);
%for kk=1:Kz
    %col_kk = num2str(F(:,kk)');
%    colStr = sprintf( '%d', F(:,kk)' );
%    val(kk) = bin2dec(colStr);        
%end

[sort_vals, sort_ind] = sort(val,'descend');

end