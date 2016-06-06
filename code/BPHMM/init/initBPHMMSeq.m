function [Psi, algParams, outParams] = initBPHMMSeq( data, model, initParams, algParams, outParams )
% Initialize BP-HMM hidden variables for MCMC sampling
%  via a **sequential** process
% That first fits a crude model to a subset of M sequences alone,
%  and then iteratively adds sequences ONE AT A TIME
%    trying to create unique features when appropriate
% Be warned: there is no guaranteed number of features that exist
%  after this procedure. It depends not only on number of sequences
%  but also on the complexity of those sequences

M = initParams.nSubsetObj;

% ------------------------------- Remember old state to use again afterward
curStream = RandStream.getGlobalStream();
entryState = curStream.State;

% Reset PRNG state to associate with specific *task*
reset( RandStream.getGlobalStream(), outParams.taskID );

% -------------------------------------------- Init BPHMM for first M sequences
if strcmp( class(data), 'SeqData' )    
    Dsubset = SeqData();
elseif strcmp( class(data), 'ARSeqData' )
    Dsubset = ARSeqData( data.R );
end
subsetObjIDs = randsample( data.N, min(M,data.N)  );
remObjIDs    = setdiff( 1:data.N, subsetObjIDs );
for ii = subsetObjIDs'
   Dsubset = Dsubset.addSeq( data.seq(ii), data.seqNames{ii} ); 
end
outP = outParams;
outP.doPrintHeaderInfo = 0;

fprintf( 'Initializing: %d unique on subset of %d sequences\n', initParams.F.nUniquePerObj, M );

% Initialize model on first M sequences
%  using a few unique featuress per sequence
subPsi = initBPHMMFresh( Dsubset, model, initParams, algParams, outP );

if ~isempty( algParams.doAnneal )
   subPsi.invTemp = 0; % forget reversibility... give me a good init! 
end

% Run through a few iters of full sampler on this small dataset
for warmiter = 1:5
    subPsi = sampleSharedFeats( subPsi, Dsubset ); 
    subPsi = sampleStateSeq( subPsi, Dsubset );
    
    subPsi = sampleSplitMerge_SeqAlloc( subPsi, Dsubset, algParams );

    subPsi.ThetaM = subPsi.ThetaM.sampleAllTheta( Dsubset, subPsi.stateSeq );
    subPsi.TransM = subPsi.TransM.sampleAllEta( subPsi.F, subPsi.stateSeq );
    
    % Try to delete extra features if possible
    algPdeath = algParams;
    algPdeath.Debug.MoveType = 0; %force death move!
    algPdeath.doAvoidCache   = 1;
    if isfield( subPsi, 'cache' )
        subPsi = rmfield( subPsi, 'cache' );
    end
    subPsi = sampleUniqueFeats( subPsi, Dsubset, algPdeath, 0 );
    subPsi = sampleStateSeq( subPsi, Dsubset );

    checkPsiFeatureConsistency( subPsi );
    
end


fprintf( 'Sampling features for remaining sequences\n' );
tic;
for rr = 1:length( remObjIDs )
    Dsubset = Dsubset.addSeq( data.seq( remObjIDs(rr) ), data.seqNames{remObjIDs(rr)} );
    curObjID = Dsubset.N;
    % Assume new object possesses all features
    %   and has prior mean on transition params
    subPsi.F(curObjID,:) = 1;
    subPsi.TransM = subPsi.TransM.getAllEta_PriorMean( 1:curObjID, subPsi.F, []);
    
    %if isfield( subPsi, 'cache' )
    %       subPsi = rmfield( subPsi, 'cache' );
    %end
    %subPsi = sampleStateSeq( subPsi, Dsubset, curObjID );
    %checkPsiFeatureConsistency( subPsi );

    % Perform MH updates to turn shared features on and off
    [subPsi] = sampleSharedFeats( subPsi, Dsubset, curObjID );
    
    % Attempt to do birth move using *entire sequence* data
    algPbirth = algParams;
    algPbirth.Debug.MoveType = 1;   
    algPbirth.Debug.wstart   = 1;
    algPbirth.Debug.wend     = Dsubset.Ts( curObjID );
    [subPsi,~, ~] = sampleUniqueFeats( subPsi, Dsubset, algPbirth, 0, curObjID );
    subPsi.TransM.K = size( subPsi.F, 2 );
    subPsi.TransM.N = size( subPsi.F, 1 );

    % Every so often, resample ALL shared feat assignments + emit params
    if rr < 3 || mod(rr,5)==0 || rr == length(remObjIDs)
        if isfield( subPsi, 'cache' )
            subPsi = rmfield( subPsi, 'cache' );
        end
        subPsi = sampleSharedFeats( subPsi, Dsubset );
        % Try to do feature death moves
        subPsi = sampleUniqueFeats( subPsi, Dsubset, algPdeath, 0 );

        if isfield( subPsi, 'cache' )
            subPsi = rmfield( subPsi, 'cache' );
        end
        subPsi = sampleStateSeq( subPsi, Dsubset );
        checkPsiFeatureConsistency( subPsi );

        subPsi.ThetaM = subPsi.ThetaM.sampleAllTheta( Dsubset, subPsi.stateSeq );
        subPsi.TransM = subPsi.TransM.sampleAllEta( subPsi.F, subPsi.stateSeq );
        
        checkPsiFeatureConsistency( subPsi );
    end
    
    if mod(rr,20)==0
        fprintf( '    seq %4d/%4d  after %4.0f sec | nF= %d\n', rr, length(remObjIDs), toc, size( subPsi.F,2) );
    end
end
 
% Finally, 
%   shuffle so that the data sequence is in the right (original) order
if isfield( subPsi, 'cache' )
    subPsi = rmfield( subPsi, 'cache' );
end
Psi = subPsi;
reorderIDs = [subsetObjIDs(:)' remObjIDs(:)'];
Psi.F = subPsi.F( reorderIDs, : );
Psi.stateSeq = subPsi.stateSeq( reorderIDs );

Psi.TransM = Psi.TransM.sampleAllEta( Psi.F, Psi.stateSeq );
Psi.ThetaM = Psi.ThetaM.sampleAllTheta( data, Psi.stateSeq );

checkPsiFeatureConsistency( Psi );
Psi = reallocateFeatIDs( Psi );
checkPsiFeatureConsistency( Psi );


% ---------------------------------------------------------  Reset stream
curStream = RandStream.getDefaultStream();
curStream.State = entryState;