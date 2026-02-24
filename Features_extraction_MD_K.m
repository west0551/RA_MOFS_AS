function [MD_K_up, MD_K_down, MD_K_total] = Features_extraction_MD_K(ppgSeg, SPT)
% 提取 PPG 片段的 K 值特征（幅值分布比值）
% 输入：
%   ppgSeg：结构体数组，字段包含 .s（PPG 幅值）
%   SPT：峰值点索引数组（每个片段一个）
% 输出：
%   MD_K_up、MD_K_down、MD_K_total：1×N 数组，分别代表上升期、下降期、整体 K 特征

    % N = length(ppgSeg);
    N=5;
    % 初始化输出变量
    MD_K_up = zeros(1, N);
    MD_K_down = zeros(1, N);
    MD_K_total = zeros(1, N);

    for k = 1:N
        % 上升期：SPT(k) 到末尾
        seg_up = ppgSeg(k).s(1:SPT(k));
        max_up = max(seg_up);
        min_up = min(seg_up);
        mean_up = mean(seg_up);
        if max_up > min_up
            MD_K_up(k) = (mean_up - min_up) / (max_up - min_up);
        else
            MD_K_up(k) = NaN;  % 避免除以 0
        end
        % 下降期：开头到 SPT(k)
        seg_down = ppgSeg(k).s(SPT(k):end);
        max_down = max(seg_down);
        min_down = min(seg_down);
        mean_down = mean(seg_down);
        if max_down > min_down
            MD_K_down(k) = (mean_down - min_down) / (max_down - min_down);
        else
            MD_K_down(k) = NaN;
        end

        % 整体段
        seg_total = ppgSeg(k).s;
        max_total = max(seg_total);
        min_total = min(seg_total);
        mean_total = mean(seg_total);
        if max_total > min_total
            MD_K_total(k) = (mean_total - min_total) / (max_total - min_total);
        else
            MD_K_total(k) = NaN;
        end
    end


end
