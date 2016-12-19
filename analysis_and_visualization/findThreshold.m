function [dataThresh,expval] = findThreshold(initData,devCount,minVal)
%%
if isempty(mfilename) || nargin == 0
    initData = ctrlData*inverseVal;
    devCount = 2;
    minVal = 2;
    close all
end
testCt = 1000;
testvals = linspace(0.0005,1,testCt);
ptally = zeros(testCt,1);
for i = 1:testCt
    testData = initData;
%     testData(testData < minVal) = [];
    testData = testData-minVal;
    testData = log1p(testData);
    expval = testvals(i);
    testData = testData.^expval;
%     [~,p] = kstest(testData);
%     p = abs(skewness(testData));
    p = abs(kurtosis(testData));
    ptally(i) = p;
end
if isempty(mfilename) || nargin == 0
    plotData = ptally;
    plot(plotData,'.')
end
%%
[~,maxval] = min(ptally);
i = maxval(1);
testData = initData;
% testData(testData < minVal) = [];
testData = testData-minVal;
testData = log1p(testData);
expval = testvals(i);
testData = testData.^expval;

if isempty(mfilename) || nargin == 0
    plotData = (testData);
    ksdensity(plotData)
    densities = ksdensity(plotData);
    hold all
    plot(plotData,rand(numel(plotData),1)*max(densities),'.')
%     plot(plotData,rand(numel(plotData),1),'.')
    % set(gca,'ylim',[0 0.1])
end
if minVal == 0
    dataThresh = mean(testData)+std(testData)*devCount;
else
    dataThresh = mean(testData)+std(testData)*devCount;
end

