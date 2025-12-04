function AngleFoot = sigToFill(AngleFoot,posToFill)    
    portionToFill = [];
    s = find(diff(posToFill)>1);
    if ~isempty(s)
        for i = 1:size(s,1)
            if i == 1
                portionToFill = [portionToFill; posToFill(1),posToFill(s(i))];
            else
                portionToFill = [portionToFill; posToFill(s(i-1)+1),posToFill(s(i))];
            end
        end
        portionToFill = [portionToFill; posToFill(s(i)+1),posToFill(end)];
    else
        portionToFill = [portionToFill; posToFill(1),posToFill(end)];
    end
    
    for i = 1:size(portionToFill,1)
        if portionToFill(i,1)>1
            AngleFoot(portionToFill(i,1):portionToFill(i,2))=...
                ones(length(portionToFill(i,1):portionToFill(i,2)),1)*AngleFoot(portionToFill(i,1)-1);
        else
            if portionToFill(i,end)<length (AngleFoot)
                AngleFoot(portionToFill(i,1):portionToFill(i,2))=...
                    ones(length(portionToFill(i,1):portionToFill(i,2)),1)*AngleFoot(portionToFill(i,2)+1);
            else
                warning('signal cannot be filled');
            end
        end
    end
end