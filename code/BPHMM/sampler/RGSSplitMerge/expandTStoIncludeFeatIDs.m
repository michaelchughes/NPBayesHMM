function [propTS] = expandTStoIncludeFeatIDs( F, TS, hyperparams, activeObjIDs, featIDs )

alpha0 = hyperparams.alpha;
kappa0 = hyperparams.kappa;

propTS = TS;

for aa = activeObjIDs

    if F( aa, featIDs(1) ) &&  F( aa, featIDs(2) )
        continue;
    end
    
    curFeatIDs = TS.obj(aa).availFeatIDs;
    Kz_ii = length( curFeatIDs );    
    propTS.obj(aa).pi_init = ones( 1, Kz_ii+1 );
    Pz = zeros(Kz_ii+1, Kz_ii+1);
    
    Pz(1:Kz_ii, 1:Kz_ii) = propTS.obj(aa).pi_z;
    Pz(end,:) = randgamma( alpha0*ones(1,Kz_ii+1) );
    Pz(:, end) = randgamma( alpha0*ones(1,Kz_ii+1) );
    Pz(end,end) = randgamma( alpha0 +kappa0 );

    if F( aa, featIDs(1) )
        kadd = featIDs(2);
    elseif F( aa, featIDs(2) )    
        kadd = featIDs(1);
    else
        error( 'BAD stuff happening' );
    end

    [activeFeatIDs, sortIDs] = sort( [TS.obj(aa).availFeatIDs kadd] );
    propTS.obj(aa).availFeatIDs = activeFeatIDs;
    propTS.obj(aa).pi_z = Pz( sortIDs, sortIDs );
end