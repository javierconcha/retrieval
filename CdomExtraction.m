function a_CDOM = CdomExtraction(CDOM_RawCurve)


ODcel = CDOM_RawCurve(:,2)-CDOM_RawCurve(CDOM_RawCurve(:,1)==800,2); % direct from the filter w/o biass

L =0.1; % [m] (Cell length = 10cm)
a_CDOM = 2.303*ODcel/L;