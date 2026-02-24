function [ ppgSeg ] = function_ppgSeg( ppgs_preprocess,fps,interval,ppgt )

ppg1derive = diff(ppgs_preprocess);

[m,index]=findpeaks(ppg1derive,'MinPeakHeight',0);%find all the ascending slopes

% 找一阶导最大值点
ii=0;
crit = round(interval*2/3);
temp = index(1);
for k=1:(length(index)-1)
    if (index(k+1)-temp)<crit
        if m(k+1) > m(index==temp)
            temp = index(k+1);
        end
    else
        ii=ii+1;
        c(ii)=temp;
        temp = index(k+1);
    end
end


%% search the beginning points of pulses in one interval before the systolic ascending slopes
p = zeros(1,length(c));

%find the foot of first segment
if c(1)>=3
    [~,b]=findpeaks(-ppgs_preprocess(1:c(1)));
    if ~isempty(b)
                p(1)=b(end);
            else
                p(1)=nan;
    end
else
    p(1)=nan;
end

%find the beginning points of every pulse
for k=2:length(c)
    [~,q]=findpeaks(-ppgs_preprocess(c(k-1)+1:c(k)));
    if ~isempty(q)
        p(k)=q(end)+c(k-1);
    else
        p(k)=nan;
    end
end

p(isnan(p))=[];

% figure;
% plot(ppgt,ppgs_preprocess,ppgt(p),ppgs_preprocess(p),'r*');


%% segmentation base on the beginning points
ppgSeg(1).t=ppgt(1:p(1));
ppgSeg(1).s=ppgs_preprocess(1:p(1));
for x=2:length(p)
    ppgSeg(x).t=ppgt(p(x-1)+1:p(x));
    ppgSeg(x).s=ppgs_preprocess(p(x-1)+1:p(x));
%     plot(ppgSeg(x).t,ppgSeg(x).s)
end

ppgSeg(1)=[];%the first segment is always not complete

%% remove the segment that not suitable

mm=zeros(1,length(ppgSeg));
mmi=zeros(1,length(ppgSeg));

for rr=1:length(ppgSeg)
    mm(rr)=max(ppgSeg(rr).s);
    mmi(rr)=min(ppgSeg(rr).s);
end
mm_average=mean(mm);
mmi_average=mean(mmi);

r=1;
while r<=length(ppgSeg)
    ma = findpeaks(ppgSeg(r).s);
    mai = min(ppgSeg(r).s);
    if ma(1)>2.5*mm_average||ma(1)<0.3*mm_average ...
            ||mai<2.5*mmi_average||mai>0.3*mmi_average ...
            ||length(ppgSeg(r).s)>round(interval*4/3)|| ...
            length(ppgSeg(r).s)<round(interval*2/3)|| ...
            ma(1) <= 0.8*max(ppgSeg(r).s)
                       % ma(1) ~= max(ppgSeg(r).s)
%            abs(ppgSeg(r).s(1)-ppgSeg(r).s(end))>0.5*(max(ppgSeg(r).s)-min(ppgSeg(r).s)) ...

        ppgSeg(r)=[];
        r=r-1;
    end
    r=r+1;
end 


end

