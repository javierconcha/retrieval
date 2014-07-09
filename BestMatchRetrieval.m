function [XResults,IMatrix] = BestMatchRetrieval(waterpixels,LUT,LUTconc)
%% Best Match
tic

matlabpool open 4

parfor index = 1:size(waterpixels,1)
    WaterPxMatrix = ones(size(LUT,1),1)*waterpixels(index,1:5);
    RMSE = sqrt(mean((LUT(:,1:5)-WaterPxMatrix).^2,2));
    [~,I] = min(RMSE);
    XResults(index,:) = LUTconc(I,:);
    IMatrix(index) = I;
    
end

IMatrix = IMatrix';
matlabpool close

disp('Elapsed time is (min):')
disp(toc/60)
