mu = mean(ctrlData);
skew = skewness(ctrlData);
kurt = kurtosis(ctrlData);
m = numel(ctrlData);
n = 1;
ctrl25 = prctile(ctrlData,25);
ctrl75 = prctile(ctrlData,75);
testCt = 1000;
testvals = linspace(0.0005,range(ctrlData)/2,testCt);
ptally = zeros(testCt,2);
for i = 1:testCt
    close all
    testData = ctrlData;
    sigma = testvals(i);
    r = pearsrnd(mu,sigma,skew,kurt,m,n);
%     plotData = r;
%     plot(plotData,rand(numel(plotData),1),'.')
%     hold all
%     plotData = ctrlData;
%     plot(plotData,rand(numel(plotData),1)+1,'.')
%     drawnow
    test25 = prctile(r,25);
    test75 = prctile(r,75);
    ptally(i,:) = [test25/ctrl25 test75/ctrl75];
end
%%
[~,minRef] = min(abs(mean(ptally,2)-1));
sigma = testvals(minRef);
plot(abs(mean(ptally,2)-1))
%%
Y = random('Normal',mu,sigma,1000000,1);
% plot(Y,'.')
prctile(Y,97.7)
sigma*2



