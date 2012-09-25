restoredefaultpath;
addpath( genpath('.') );
lightspeedDir = getUserSpecifiedPath('LightspeedLibrary');
addpath( genpath( lightspeedDir ) );

try
    randgamma( 100 );
catch e
    fprintf( 'Please double check user specified LightspeedLibrary.path\n' );
    fprintf( 'Current path: %s\n', lightspeedDir );
    error( 'ERROR: Lightspeed toolbox not properly linked.\n' );

end