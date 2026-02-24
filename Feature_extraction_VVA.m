function [VVA_T_up_med, VVA_T_down_med, VVA_width_med, ...
          VVA_width_ratio_LR_med, VVA_width_ratio_med, ...
          VVA_amplitude_med, VVA_width_triangle_area_med, ...
          VVA_width_triangle_area_ratio_med, ...
          VVA_area_up_real_med, VVA_area_down_real_med, ...
          VVA_area_ratio_LR_med, VVA_area_ratio_up_med, ...
          VVA_area_ratio_down_med, VVA_area_total_med, ...
          VVA_area_ratio_total_med] = Feature_extraction_VVA(ppgSeg_norm, SPT)



% 输入：
%   ppgSeg_norm：结构体数组，每个包含已归一化的 s、t
%   SPT：峰值点索引数组，长度与 ppgSeg_norm 相同
%
% 输出：
%   每个输出是 N×19 的矩阵（N为片段数量，19是幅值等级数）

    num_levels = 20;
    % N = length(ppgSeg_norm);
    N=5;
    num_pts = num_levels - 1;

    % 初始化所有输出矩阵（N 行，19 列）
    VVA_T_up_all = zeros(N, num_pts);
    VVA_T_down_all = zeros(N, num_pts);
    VVA_width_all = zeros(N, num_pts);
    VVA_width_ratio_LR_all = zeros(N, num_pts);
    VVA_width_ratio_all = zeros(N, num_pts);
    VVA_amplitude_all = zeros(N, num_pts);
    VVA_width_triangle_area_all = zeros(N, num_pts);
    VVA_width_triangle_area_ratio_all = zeros(N, num_pts);

    VVA_area_up_real = zeros(N, num_pts);
    VVA_area_down_real = zeros(N, num_pts);
    VVA_area_ratio_LR = zeros(1, num_pts);
    VVA_area_ratio_up = zeros(1, num_pts);
    VVA_area_ratio_down = zeros(1, num_pts);
    VVA_area_total = zeros(1, num_pts);
    VVA_area_ratio_total = zeros(1, num_pts);


    for k = 1:N
        s = ppgSeg_norm(k).s;
        t = ppgSeg_norm(k).t;
        peak_idx = SPT(k);
        FA = s(peak_idx);

        % 等幅值（0.05 ~ 0.95）
        h = linspace(0, FA, num_levels + 1);
        h = h(2:end);

        % 上升期与下降期
        s_up = s(1:peak_idx); t_up = t(1:peak_idx);
        s_down = s(peak_idx:end); t_down = t(peak_idx:end);

        % 插值获得时间点
        t_interp_up = interp1(s_up, t_up, h, 'linear', 'extrap');
        t_interp_down = interp1(s_down, t_down, h, 'linear', 'extrap');

        % ppgSeg_k=k;
        % % 提取第2个PPG片段数据
        % s = ppgSeg_norm(ppgSeg_k).s;
        % t = ppgSeg_norm(ppgSeg_k).t;
        % 
        % % 插值点（注意加 :）
        % % h = h(2);   % 等幅值（高度）
        % t_up = t_interp_up;     % 上升期插值时间
        % t_down = t_interp_down; % 下降期插值时间
        % 
        % % 开始绘图
        % figure; hold on; grid on;
        % 
        % % 绘制PPG原始波形
        % plot(t, s, 'b-', 'LineWidth', 1.5);
        % 
        % % 标出等幅值插值点
        % plot(t_up, h, 'ro', 'MarkerSize', 6, 'DisplayName', '上升期插值点');
        % plot(t_down, h, 'go', 'MarkerSize', 6, 'DisplayName', '下降期插值点');
        % 
        % % 加文字标注（可选）
        % for i = 1:length(h)
        %     text(t_up(i), h(i)+0.02, sprintf('%.2f', h(i)), 'Color', 'r', 'FontSize', 8);
        %     text(t_down(i), h(i)+0.02, sprintf('%.2f', h(i)), 'Color', 'g', 'FontSize', 8);
        % end
        % 
        % % 图形美化
        % xlabel('时间 t (s)');
        % ylabel('PPG归一化幅值');
        % title('PPG信号与上升/下降等幅值点插值标注');
        % legend({'PPG原始波形', '上升期点', '下降期点'}, 'Location', 'best');

        % 时间相关参数
        t_peak = t(peak_idx);
        t_start = t(1);
        t_end = t(end);
        duration_total = t_end - t_start;

        % 整个上升段面积（从起点到峰值）
        t_up_total = linspace(t_start, t_peak, 200);
        s_up_total = interp1(t, s, t_up_total, 'linear');
        area_up_total = trapz(t_up_total, s_up_total);

        % 整个下降段面积（从峰值到末端）
        t_down_total = linspace(t_peak, t_end, 200);
        s_down_total = interp1(t, s, t_down_total, 'linear');
        area_down_total = trapz(t_down_total, s_down_total);

        % 整段 PPG 面积
        area_total_full = trapz(t, s);

        for i = 1:num_pts
            % 表示血管在扩张过程中，从起始点扩张到当前状态需要的时间
            VVA_T_up_all(k, i) = t_peak - t_interp_up(i);
            % 表示血管在收缩过程中，从峰值点收缩到当前状态需要的时间
            VVA_T_down_all(k, i) = t_interp_down(i) - t_peak;
            % 宽度特征
            VVA_width_all(k, i) = VVA_T_up_all(k, i) + VVA_T_down_all(k, i);
            % 上下对称性比值
            % 避免除以0
            if VVA_T_down_all(k, i) ~= 0
                VVA_width_ratio_LR_all(k, i) = VVA_T_up_all(k, i) / VVA_T_down_all(k, i);
            else
                VVA_width_ratio_LR_all(k, i) = NaN;
            end
            % 宽度占整段时长的比例
            VVA_width_ratio_all(k, i) = VVA_width_all(k, i) / duration_total;
            % 血管当前归一化幅值
            VVA_amplitude_all(k, i) = h(i);
            % 简化三角形面积
            VVA_width_triangle_area_all(k, i) = h(i) * VVA_width_all(k, i) / 2;
            % 面积占最大理论面积比例
            % 避免除以0
            if FA ~= 0 && duration_total ~= 0
                VVA_width_triangle_area_ratio_all(k, i) = ...
                    VVA_width_triangle_area_all(k, i) / (FA * duration_total / 2);
            else
                VVA_width_triangle_area_ratio_all(k, i) = NaN;
            end

             % 真实面积（上升期）
            t_range_up = linspace(t_interp_up(i), t_peak, 50);  % 50个点细化
            s_range_up = interp1(t, s, t_range_up, 'linear');
            VVA_area_up_real(k, i) = trapz(t_range_up, s_range_up);

            % 真实面积（下降期）
            t_range_down = linspace(t_peak, t_interp_down(i), 50);
            s_range_down = interp1(t, s, t_range_down, 'linear');
            VVA_area_down_real(k, i) = trapz(t_range_down, s_range_down);

            A_up = VVA_area_up_real(k,i);
            A_down = VVA_area_down_real(k,i);
            A_all = A_up + A_down;

            % 新增特征计算
            VVA_area_ratio_LR(i) = A_up / A_down;
            VVA_area_ratio_up(i) = A_up / area_up_total;
            VVA_area_ratio_down(i) = A_down / area_down_total;
            VVA_area_total(i) = A_all;
            VVA_area_ratio_total(i) = A_all / area_total_full;
        end
    end
        VVA_T_up_med = median(VVA_T_up_all, 1);
        VVA_T_down_med = median(VVA_T_down_all, 1);
        VVA_width_med = median(VVA_width_all, 1);
        VVA_width_ratio_LR_med = median(VVA_width_ratio_LR_all, 1);
        VVA_width_ratio_med = median(VVA_width_ratio_all, 1);
        VVA_amplitude_med = median(VVA_amplitude_all, 1);
        VVA_width_triangle_area_med = median(VVA_width_triangle_area_all, 1);
        VVA_width_triangle_area_ratio_med = median(VVA_width_triangle_area_ratio_all, 1);
        VVA_area_up_real_med = median(VVA_area_up_real, 1);
        VVA_area_down_real_med = median(VVA_area_down_real, 1);
        VVA_area_ratio_LR_med = median(VVA_area_ratio_LR, 1);
        VVA_area_ratio_up_med = median(VVA_area_ratio_up, 1);
        VVA_area_ratio_down_med = median(VVA_area_ratio_down, 1);
        VVA_area_total_med = median(VVA_area_total, 1);
        VVA_area_ratio_total_med = median(VVA_area_ratio_total, 1);
end

