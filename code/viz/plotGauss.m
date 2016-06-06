function [] = plotGauss( Mu, Sigma, doInv, MyColors, figH )
% This function plots the contours of a 2-D gaussian, on current figure
% For any Gaussian, we can use the eigenvectors of covariance Sigma
%    to find ellipses which have "constant probability" under the pdf
% INPUT
%    Mu   :   KxD  matrix of means (each row is a mean)
%   Sigma :   DxDxK matrix of covar. matrices (each 3d slice is cov mat)
%   doInv :   set to zero if Sigma is a covariance,
%                      1  if Sigma is a precision matrix (needs inversion)
%   MyColors : [R G B] triple color to use for the contour plot

if ~exist( 'doInv', 'var' )
    doInv =0;
end

if ~exist( 'MyColors', 'var' );
    hold all;
    MyColors = jet( 100 );
    MyColors = MyColors( randsample( 100, size(Mu,1) ) , : );
else
    hold on;
end
if ~exist( 'figH', 'var' )
    figH = gca();
end

[K,D] = size( Mu );
for kk = 1:K
    if doInv
        Sigma(:,:,kk) = Sigma(:,:,kk) \ eye(D);
    end
end
Sigma = Sigma(1:2, 1:2, :);
for kk = 1:K
    taskVis = 'on';
    %for R = linspace(0.1, 2, 4)
    for R = linspace(0.5, 2, 2)
         
        t = -pi:.01:pi;
        k = length(t);
        x = R*sin(t);
        y = R*cos(t);
        
        
        [vv,dd] = eig( Sigma(:,:,kk) );
        A = real((vv*sqrt(dd))');
        z = [x' y']*A;
        
        plot(figH, z(:,1)+Mu(kk,1),z(:,2)+Mu(kk,2), '.', 'MarkerSize', 15, 'Color', MyColors(kk,:), 'HandleVisibility', taskVis);
        taskVis = 'off';
    end
    
end
