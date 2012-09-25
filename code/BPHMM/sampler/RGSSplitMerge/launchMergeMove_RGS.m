function [propPsi] = launchMergeMove_RGS( Psi, data_struct, hyperparams, model, featIDs )

propF = Psi.F;

ActiveObjIDs = find( sum(Psi.F(:,featIDs),2) > 0 )';
propF( :, featIDs ) = 0;
propF( ActiveObjIDs, end+1 ) = 1;

kmerge = size( propF, 2);
propFeatIDs = [kmerge kmerge];



propStateSeq = Psi.stateSeq;
for ii = 1:length( propStateSeq )
    if propF(ii,kmerge) == 0
        continue;
    end
    
    aIDs =  propStateSeq(ii).z == featIDs(1);
    bIDs =  propStateSeq(ii).z == featIDs(2);
    propStateSeq(ii).z( aIDs | bIDs ) = kmerge;
end

% ----------------------------------------------- Sample TS and theta    
propTheta = sampleTheta_RGS( propF, propStateSeq, Psi.theta, data_struct, model, propFeatIDs );
propTS    = sampleTransParams_RGS( propF, propStateSeq, Psi.TS, hyperparams, model, ActiveObjIDs );

propPsi.activeFeatIDs = propFeatIDs;
propPsi.F = propF;
propPsi.stateSeq = propStateSeq;
propPsi.theta = propTheta;
propPsi.TS = propTS;