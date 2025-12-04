function plotEvents(TurnM, TurnDur, maxV, minV)
    for i = 1:size(TurnM,1)
%         for s = 1:2         
%             line([TurnDur(i,s),TurnDur(i,s)],[minV, maxV],'Color','b','LineStyle',':','LineWidth',1)
%         end
        X2=([TurnDur(i,1), TurnDur(i,2)]);
        Y2=([minV,minV]);
        Y3=([maxV,maxV]);
        I = patch([X2 fliplr(X2)],[Y2 fliplr(Y3)], 'b');
        alpha(0.01);
        
        pos = maxV-50*(i-1);
        value = num2str(TurnM(i));
        value = value(1:(max(strfind(value, '.')))+1);
        text(TurnDur(i,1),pos,value)
    end
end
