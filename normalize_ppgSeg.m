function ppgSeg_norm = normalize_ppgSeg(ppgSeg)
% 对 ppgSeg 中的每段信号 s 进行 Min-Max 归一化
% 输入：ppgSeg（结构体数组，包含字段 .s）
% 输出：ppgSeg_norm（结构体数组，.s 为归一化后信号，其它字段保留）

    ppgSeg_norm = ppgSeg;  % 复制结构体框架

    for k = 1:length(ppgSeg)
        s = ppgSeg(k).s;
        s_min = min(s);
        s_max = max(s);

        if s_max > s_min
            s_norm = (s - s_min) / (s_max - s_min);
        else
            s_norm = zeros(size(s));  % 处理平坦信号
        end

        ppgSeg_norm(k).s = s_norm;  % 替换为归一化结果
    end
end

% %Z-score 标准化
%     for k = 1:length(ppgSeg)
%     s = ppgSeg(k).s;
%     mu = mean(s);
%     sigma = std(s);
%     if sigma > 0
%         ppgSeg(k).s = (s - mu) / sigma;
%     else
%         ppgSeg(k).s = zeros(size(s));  % 防止除以 0
%     end
%     end