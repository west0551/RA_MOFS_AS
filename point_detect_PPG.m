function [SP,SPT,DP,DPT,DN,DNT] = point_detect_original(ppgSeg)


SP = zeros(1,length(ppgSeg));
SPT = zeros(1,length(ppgSeg));
DP = zeros(1,length(ppgSeg));
DPT = zeros(1,length(ppgSeg));
DN = zeros(1,length(ppgSeg));
DNT = zeros(1,length(ppgSeg));


for k=1:length(ppgSeg)

    [SP(k), index]=max(ppgSeg(k).s);
    %SPT(k) = ppgSeg(k).t(index);
    SPT(k) = index;

    %h=(m-min(ppgSeg(k).s))*0.01+min(ppgSeg(k).s);
    h=min(ppgSeg(k).s);
    [v,index2]=findpeaks(ppgSeg(k).s(index+1:end),'MinPeakHeight',h);
    [h,index3]=findpeaks(-ppgSeg(k).s(1:round(0.95*end)));

        if ~isempty(v)
           DP(k) = v(1);
           %DPT(k) = ppgSeg(k).t(index2(1)+index);
           DPT(k) = index2(1)+index;
        else 
           DP(k) = nan;
           DPT(k) = nan;
        end

        if ~isempty(h)
           DN(k) = -h(1);
           %DNT(k) = ppgSeg(k).t(index3(1));
           DNT(k) = index3(1);
           else 
           DN(k) = nan;
           DNT(k) = nan;
        end

end


