function [propPsi] = launchSplitMove_RGS( Psi, data_struct, hyperparams, model, anchorObjIDs, featIDs )

kold = featIDs(1);

ActiveObjIDs = find( Psi.F( :, kold ) > 0 )';

propF = Psi.F;
knew = size( propF,2) + 1;

propF( :, knew ) = 0;
propF( ActiveObjIDs, knew ) = 1;

propFeatIDs = [kold knew];


propStateSeq = Psi.stateSeq;

for aa = ActiveObjIDs
    oldts = propStateSeq( aa ).z == kold;
    
    if aa == anchorObjIDs(1)
        propStateSeq( aa ).z( oldts ) = kold;
    elseif aa == anchorObjIDs(2)
        propStateSeq( aa ).z( oldts ) = knew;
    else
        if rand > 0.5
            propStateSeq( aa ).z( oldts ) = kold;
        else
            propStateSeq( aa ).z( oldts ) = knew;
        end
    end    
end

propTheta = sampleTheta_RGS( propF, propStateSeq, Psi.theta, data_struct, model, propFeatIDs );
propTS    = sampleTransParams_RGS( propF, propStateSeq, Psi.TS, hyperparams, model, ActiveObjIDs );

propPsi.activeFeatIDs = propFeatIDs;
propPsi.F = propF;
propPsi.stateSeq = propStateSeq;
propPsi.theta = propTheta;
propPsi.TS = propTS;