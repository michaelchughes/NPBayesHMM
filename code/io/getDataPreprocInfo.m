function Preproc = getDataPreprocInfo( dataName, dataParams )
% Obtain struct of "preprocess" parameters 
%  these indicate how to create/transform raw data
%  into desired format to be processed by the given time series model
% When adding a new dataset, be sure to EDIT THIS FILE

MocapChannelNames = {'root.ty', 'lowerback.rx', 'lowerback.ry', 'upperneck.ry', ...
                     'rhumerus.rz', 'rradius.rx','lhumerus.rz', 'lradius.rx', ...
                     'rtibia.rx', 'rfoot.rx', 'ltibia.rx', 'lfoot.rx'...    
                    };
                
KitchenMocapChannelNames = ...
    {'root.ty', 'lowerback.rx', 'lowerback.ry', 'neck.ry', ...
     'ltibia.rx', 'lfoot.rx', 'rtibia.rx', 'rfoot.rx', ...     
     'lhumerus.rx', 'lhumerus.rz', 'lradius.rz', 'lhand.rx', 'lhand.rz', ...
     'rhumerus.rx', 'rhumerus.rz', 'rradius.rz', 'rhand.rx', 'rhand.rz'};
     
               
% =============================================== DEFAULT PARAMS
PP = struct(); 
switch lower(dataName)
    % ----------------------------------------------------- Synth Data
    case {'synth', 'synthmultinomial', 'synthdiscrete'}
        PP.nObj = 50;
        PP.T = 100;
        PP.nStates = 10;
        PP.V = 1000; % # vocabulary symbols
        PP.obsDim = -5;
        PP.pEmitFavor = 0.9; % prob. state emits symbol from its favored set
    case {'synthbernoulli','synthbinary'}
        PP.nObj = 50;
        PP.T = 100;
        PP.nStates = 5;
        PP.nDims = 40;
    case {'synthgaussian', 'synthnormal'}
        PP.nObj = 50;
        PP.T = 100;
        PP.nStates = 8;
        PP.obsDim = 5;
    case {'synthar', 'synthautoreg'}
        PP.nObj = 50;
        PP.T = 100;
        PP.nStates = 8;
        PP.obsDim = 2;
        PP.R = 1;
    % ----------------------------------------------------- Video Clips
    case {'kthfinal'}
        PP.nObj = 0;
        PP.V = 1000;
        PP.detType = 'STIP';
        PP.featureName = 'HOF';
        PP.timeBinSize = 2/25; 
                
    case {'kthmashup'}
        PP.nObj = 0;
        PP.V = 1000;
        PP.detType = 'STIP';
        PP.featureName = 'HOF';
        PP.timeBinSize = 4/25; 
                
    case {'olympicsports', 'olympics8'}
        PP.nObj = 0; % zero means use every sequence available!
        PP.V = 1000;
        PP.detType = 'STIP';
        PP.featureName = 'HOG+HOF';
        PP.timeBinSize = 4/25; %use 4 frames per block, given 25 fps
                
    case {'cmukitchen'}
        PP.nObj = 0; % zero means use every sequence available!
        PP.V = 1000;
        PP.detType = 'STIP';
        PP.featureName = 'HOG+HOF';
        PP.timeBinSize = 0.5;
                
    case {'kitchenbig'}
        PP.nObj = 0; % zero means use every sequence available!
        PP.V = 1000;
        PP.detType = 'STIP';
        PP.featureName = 'HOG+HOF';
        PP.timeBinSize = 1;        
    % ----------------------------------------------------- Motion Capture
    case {'kitchenmocap'}
        PP.channelNames = KitchenMocapChannelNames;
        PP.R = 1;
        PP.windowSize = 12;
        PP.obsDim = length( KitchenMocapChannelNames );
    case {'mocapbig'}
        PP.obsDim =12;
        PP.R = 1;
        PP.windowSize = 12;
        %PP.channelIDs = [2 4 5 7 11 12 15 16 20 21 25 26];
        PP.channelNames = MocapChannelNames;
    case {'mocap6'}
        PP.nObj = 6;
        PP.obsDim = 12;
        PP.R = 1;
        PP.windowSize = 12;
        %PP.channelIDs = [2 4 5 7 11 12 15 16 20 21 25 26];
        PP.channelNames = MocapChannelNames;        
end

% =============================================== UPDATE WITH USER INPUT
Preproc = updateParamsWithUserInput( PP, dataParams );
