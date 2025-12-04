function selStride = stride_direction(selStride, MrkL, side, Data)
Heel = MrkL.HEEL;
TOE = MrkL.TOE;
COM = MrkL.COM;
HeelG = Data.Mrks.(strcat(side,'HEEL'));
    for i = 1:size(selStride,2)
        IC_now = selStride(i).ICEvents;
        Heel_now = Heel(IC_now(1,1):IC_now(1,2),:);
        HeelG_now = HeelG(IC_now(1,1):IC_now(1,2),:);
        RangeAP = range(Heel_now(:,1));
        RangeML = range(Heel_now(:,3));
        if RangeAP > RangeML
            selStride(i).ForwardsStride = true;
            hFig = figure;
            plot(Heel_now(:,1))
            hold on
            plot(Heel_now(:,3))
            plot(HeelG_now(:,end))
            legend('AP_L','ML_L','V_G')
            close(hFig)
        else
            selStride(i).ForwardsStride = false;
            hFig = figure;
            plot(Heel_now(:,1))
            hold on
            plot(Heel_now(:,3))
            plot(HeelG_now(:,end))
            legend('AP_L','ML_L','V_G')
            close(hFig)
        end
    end
t = 1;
end