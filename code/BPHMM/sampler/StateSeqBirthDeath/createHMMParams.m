function [hatTS, hatTheta] = createHMMParams( F, stateSeq, hyperparams, data_struct, model, featIDs, WindowParams )

objIDs = 1:size(F,1);

fIDs = false( 1, size(F,2) );
fIDs( featIDs ) = true;
F(:, ~fIDs ) = 0;

hatTS = getEta_PriorMean( F, [], hyperparams, model, objIDs, featIDs );

if any( featIDs == WindowParams.featID )
   stateSeq(WindowParams.objID).z( WindowParams.zIDs ) = WindowParams.featID; 
end

[Xstats] = getXSuffStats( F, stateSeq, data_struct, model, objIDs, featIDs );
hatTheta = setThetaToPosteriorMean( [], Xstats, model.obsModel, featIDs );
