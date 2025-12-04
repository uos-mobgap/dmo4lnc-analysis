function [Event_IC, Event_FC] = trueAndFalseEvent(IC, FC, deltaEventHS, deltaEventTO)

Event_IC.TP = sum(~isnan(deltaEventHS));
Event_IC.FP = length(IC)-sum(~isnan(deltaEventHS));
Event_IC.FN = sum(isnan(deltaEventHS));

Event_FC.TP = sum(~isnan(deltaEventTO));
Event_FC.FP = length(FC)-sum(~isnan(deltaEventTO));
Event_FC.FN = sum(isnan(deltaEventTO));
end