function [bi,minNIR] = glint_removal(ROIbandi,ROIbandNIR)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Method base in the Technical Note: Simple and robust removal of sun glint
% for mapping shallow-water benthos, by J. D. Hedley, A. R. Harborne, and P.
% J. Mumby
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Determining the slope

% Ri = zeros(size(ROIbandi));

        
        a = ROIbandNIR(:);
        b = ROIbandi(:);
        
        %% Deglinting
            p = polyfit(double(a(a~=0)),double(b(a~=0)),1);
            bi = p(1);
            minNIR = min(min(a(a~=0)));
%             Ri(:,column:column+stripsize-1) = bandioriginal(:,column:column+stripsize-1)-bi.*(ROIbandNIR(:,column:column+stripsize-1)-MinNIR);
