function Dpost = preprocessMocapData(D,windowSize)
% Block average sensor channel matrix D
%   using given windowSize

numWindows = floor(size(D,2)/windowSize);
endWindow = rem(size(D,2),numWindows);

Dpost = zeros(size(D,1),numWindows+1);

for ii=1:numWindows
    Dpost(:,ii) = mean(D(:,(ii-1)*windowSize+1:ii*windowSize),2);
end

if endWindow>1
    Dpost(:,numWindows+1) = mean(D(:,numWindows*windowSize+1:end),2);
else
    Dpost(:,numWindows+1) = D(:,end);
end

return;