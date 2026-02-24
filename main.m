folder_path = 'C:\Users\SZY\Desktop\signal preprocess and feature extraction\Sample';
folder_path_save = 'C:\Users\SZY\Desktop\signal preprocess and feature extraction\Date_feature.xlsx';
row_num = 2;  % 从第二行开始写数据，第一行保留给header
% 获取所有MAT文件
mat_files = dir(fullfile(folder_path, '*.mat'));

% 检查是否有MAT文件
if isempty(mat_files)
    error('文件夹中没有找到.mat文件！');
end

for files_i = 1:length(mat_files)

     close all;
    % 当前文件名和路径
    filename = fullfile(folder_path, mat_files(files_i).name);
    matData = load(filename);
    % 检查是否存在extracted_data变量
    if ~isfield(matData, 'PPG_data')
        warning('文件 %s 中未找到 PPG_data 变量，跳过处理', mat_files(files_i).name);
        continue;
    end
    
    data = matData.PPG_data;
    data_number=300000;
    % 检查数据列数是否为4
    if size(data, 2) ~= 4
        warning('文件 %s 数据列数不是4列，跳过处理', mat_files(files_i).name);
        continue;
    end
      % fname = mat_files(files_i).name;
      [~, fname, ~] = fileparts(mat_files(files_i).name);

    % 提取时间和PPG信号
    for data_column = 1:1
    time = data(1:data_number, 1); % 时间列
    ppg_original_IR = data(1:data_number, 2*data_column-1); % PPG信号列（近红外）
    ppg_original_Red = data(1:data_number, 2*data_column); % PPG信号列（近红外）

    ppg_right=ppg_original_IR;
    ppg_right = abs(ppg_right);
    % 反向信号（根据需要调整）
    ppg_right = -ppg_right;

    % 设计Butterworth高通滤波器
    order = 4; % 滤波器阶数
    fs = 1000; % 采样频率（单位：Hz）
    [b, a] = butter(order,  0.5/  (fs / 2), 'high');  % 高通滤波器设计
    [c, d] = butter(order, 10 / (fs / 2), 'low');  % 低通滤波器设计

    % 设计50Hz陷波滤波器(Notch Filter)
    f0 = 50; % 要滤除的频率(Hz)
    bw = 2;  % 带宽(Hz)
    [b_notch, a_notch] = iirnotch(f0/(fs/2), bw/(fs/2));

    % 使用 filtfilt 滤波
    ppg_right = medfilt1(ppg_right, 5); 
    ppg_right = filtfilt(b_notch, a_notch, ppg_right); 
    ppg_filtered_IR = filtfilt(b, a, ppg_right);
    ppg_filtered_IR = filtfilt(c, d, ppg_filtered_IR);

    filtered_diff_signal = signal_diff(ppg_filtered_IR,data_number-2);

    %计算滤波后的一阶导的二阶导，局部为0的点即为拐点
    filtered_diff2_signal = signal_diff2(filtered_diff_signal,data_number-2);
    %二阶导信号的一阶导数，用于辅助寻找二阶导的局部极值
    filtered_diff_diff2_ppg_right_IR = signal_diff(filtered_diff2_signal,data_number-2);

    ppgs = ppg_filtered_IR;%定义需要滤波的信号Right_finger
    N = length(ppgs);  % 信号长度
    f = (0:N-1) * (fs / N);      % 频率轴（0 ~ fs Hz）
    % 计算FFT（取绝对值得到幅度谱）
    ppg_fft = abs(fft(ppgs));  
    % 由于FFT对称，只需取前一半频率
    f_half = f(1:floor(N/2));
    ppg_fft_half = ppg_fft(1:floor(N/2));    
    % 找到幅度最大的频率
    [max_magnitude, max_idx] = max(ppg_fft_half);
    dominant_freq = f_half(max_idx);  % 主频率（Hz）
    interval=1/dominant_freq*1000+50;

%% 选择使用静息PPG还是释压PPG
    ppgs_preprocess=ppg_filtered_IR(1000:100000);
    ppgt = 0:(1/fs):(length(ppgs_preprocess)-1)/fs;
%PPG切片    
    ppgSeg = function_ppgSeg(ppgs_preprocess,fs,interval,ppgt);
 
%% 基于PWA方法找PPG信号特征点定位
    [SP,SPT,DP,DPT,DN,DNT] = point_detect_PPG(ppgSeg);
% 基于PWA方法找APG和VPG特征点
    [A,B,C,D,E,F,U,V,W,AT,BT,CT,DT,ET,FT,UT,VT,WT] = point_detect_APG(ppgSeg,fs,interval);
   % 补全 DPT 和 DNT，如果为空则用 WT 或 CT 替代
    for k = 1:length(DPT)
        if isnan(DPT(k)) && ~isnan(WT(k))
            DPT(k) = WT(k)+30;
        end
        if isnan(DNT(k)) && ~isnan(ET(k))
            DNT(k) = ET(k)+30;
        end
    end

% === 绘图 ===
    plot_mode=3;
    ppgSeg_num = 20;
    if plot_mode > 0       
        ppgSegsss = ppgSeg(ppgSeg_num).s;
        NN = length(ppgSegsss);
        
        % 一阶、二阶导数
        diff_ppgSegs = signal_diff(ppgSegsss, NN);
        diff2_ppgSegs = signal_diff2(diff_ppgSegs, NN);
        
        % 特征点
        sp = SPT(ppgSeg_num); dp = DPT(ppgSeg_num); dn = DNT(ppgSeg_num);
        ut = UT(ppgSeg_num); vt = VT(ppgSeg_num); wt = WT(ppgSeg_num);
        at = AT(ppgSeg_num); bt = BT(ppgSeg_num); ct = CT(ppgSeg_num);
        dt = DT(ppgSeg_num); et = ET(ppgSeg_num); ft = FT(ppgSeg_num);
        t = 1:NN;
        switch plot_mode
            case 1  % === 分开绘图 ===
                % PPG
                figure; plot(ppgSegsss, 'b', 'LineWidth', 1.5); hold on;
                plot(1, ppgSegsss(1), 'ro'); 
                if ~isnan(sp), plot(sp, ppgSegsss(sp), 'go'); end
                if ~isnan(dp), plot(dp, ppgSegsss(dp), 'ko'); end
                if ~isnan(dn), plot(dn,ppgSegsss(dn), 'mo'); end
                title(sprintf('PPG | 文件名: %s', fname));
                xlabel('Time'); ylabel('PPG'); legend('PPG', 'O','S', 'N', 'D'); grid on;

                % VPG
                figure; plot(diff_ppgSegs, 'b', 'LineWidth', 1.5); hold on;
                if ~isnan(ut), plot(ut, diff_ppgSegs(ut), 'ro'); end
                if ~isnan(vt), plot(vt, diff_ppgSegs(vt), 'go'); end
                if ~isnan(wt), plot(wt, diff_ppgSegs(wt), 'ko'); end
                title(sprintf('VPG | 文件名: %s', fname));
                xlabel('Time'); ylabel('VPG'); legend('VPG', 'U', 'V', 'W'); grid on;
        
                % APG
                figure;
                plot(diff2_ppgSegs, 'b', 'LineWidth', 1.5); hold on;
                
                legend_entries = {'APG'};  % 初始化图例
                h = [];                    % 图例句柄数组
                
                % A点
                if ~isnan(at)
                    h(end+1) = plot(at, diff2_ppgSegs(at), 'ro', 'MarkerSize', 6);
                    legend_entries{end+1} = 'A';
                end
                
                % B点
                if ~isnan(bt)
                    h(end+1) = plot(bt, diff2_ppgSegs(bt), 'go', 'MarkerSize', 6);
                    legend_entries{end+1} = 'B';
                end
                
                % C点
                if ~isnan(ct)
                    h(end+1) = plot(ct, diff2_ppgSegs(ct), 'ko', 'MarkerSize', 6);
                    legend_entries{end+1} = 'C';
                end
                
                % D点
                if ~isnan(dt)
                    h(end+1) = plot(dt, diff2_ppgSegs(dt), 'mo', 'MarkerSize', 6);
                    legend_entries{end+1} = 'D';
                end
                
                % E点
                if ~isnan(et)
                    h(end+1) = plot(et, diff2_ppgSegs(et), 'co', 'MarkerSize', 6);
                    legend_entries{end+1} = 'E';
                end
                
                % F点
                if ~isnan(ft)
                    h(end+1) = plot(ft, diff2_ppgSegs(ft), 'yo', 'MarkerSize', 6);
                    legend_entries{end+1} = 'F';
                end
                
                title(sprintf('APG | 文件名: %s', fname));
                xlabel('Time');
                ylabel('APG');
                legend(['APG', legend_entries(2:end)], 'Location', 'best');  % 忽略 'APG' 的重复句柄
                grid on;

        
            case 2  % === 两两组合图（双Y轴） + 特征点 ===
        
                % === PPG + VPG ===
                figure;
                yyaxis left
                plot(t, ppgSegsss, 'b'); ylabel('PPG'); hold on;
                if ~isnan(sp), plot(sp, ppgSegsss(sp), 'ro'); end
                if ~isnan(dp), plot(dp, ppgSegsss(dp), 'go'); end
                if ~isnan(dn), plot(dn, ppgSegsss(dn), 'ko'); end
                
                yyaxis right
                plot(t, diff_ppgSegs, 'r'); ylabel('VPG');
                if ~isnan(ut), plot(ut, diff_ppgSegs(ut), 'ro'); end
                if ~isnan(vt), plot(vt, diff_ppgSegs(vt), 'go'); end
                if ~isnan(wt), plot(wt, diff_ppgSegs(wt), 'ko'); end
        
                title(sprintf('PPG + VPG | 文件名: %s', fname)); xlabel('Time'); grid on;
        
                % === PPG + APG ===
                figure;
                yyaxis left
                plot(t, ppgSegsss, 'b'); ylabel('PPG'); hold on;
                if ~isnan(sp), plot(sp, ppgSegsss(sp), 'ro'); end
                if ~isnan(dp), plot(dp, ppgSegsss(dp), 'go'); end
                if ~isnan(dn), plot(dn, ppgSegsss(dn), 'ko'); end
        
                yyaxis right
                plot(t, diff2_ppgSegs, 'g'); ylabel('APG');
                if ~isnan(at), plot(at, diff2_ppgSegs(at), 'ro'); end
                if ~isnan(bt), plot(bt, diff2_ppgSegs(bt), 'go'); end
                if ~isnan(ct), plot(ct, diff2_ppgSegs(ct), 'ko'); end
                if ~isnan(dt), plot(dt, diff2_ppgSegs(dt), 'mo'); end
                if ~isnan(et), plot(et, diff2_ppgSegs(et), 'co'); end
                if ~isnan(ft), plot(ft, diff2_ppgSegs(ft), 'yo'); end
        
                title(sprintf('PPG + APG | 文件名: %s', fname)); xlabel('Time'); grid on;
        
                % === VPG + APG ===
                figure;
                yyaxis left
                plot(t, diff_ppgSegs, 'r'); ylabel('VPG'); hold on;
                if ~isnan(ut), plot(ut, diff_ppgSegs(ut), 'ro'); end
                if ~isnan(vt), plot(vt, diff_ppgSegs(vt), 'go'); end
                if ~isnan(wt), plot(wt, diff_ppgSegs(wt), 'ko'); end
        
                yyaxis right
                plot(t, diff2_ppgSegs, 'g'); ylabel('APG');
                if ~isnan(at), plot(at, diff2_ppgSegs(at), 'ro'); end
                if ~isnan(bt), plot(bt, diff2_ppgSegs(bt), 'go'); end
                if ~isnan(ct), plot(ct, diff2_ppgSegs(ct), 'ko'); end
                if ~isnan(dt), plot(dt, diff2_ppgSegs(dt), 'mo'); end
                if ~isnan(et), plot(et, diff2_ppgSegs(et), 'co'); end
                if ~isnan(ft), plot(ft, diff2_ppgSegs(ft), 'yo'); end
        
                title(sprintf('VPG + APG | 文件名: %s', fname)); xlabel('Time'); grid on;
        
            case 3  % === 三合一图 + 所有特征点 ===
                figure;
                plot(t, ppgSegsss, 'b-', 'LineWidth', 1.2); hold on;
                plot(t, diff_ppgSegs, 'r--', 'LineWidth', 1);
                plot(t, diff2_ppgSegs, 'g:', 'LineWidth', 1);
        
                % 标 PPG 点
                if ~isnan(sp), plot(sp, ppgSegsss(sp), 'ro'); end
                if ~isnan(dp), plot(dp, ppgSegsss(dp), 'go'); end
                if ~isnan(dn), plot(dn, ppgSegsss(dn), 'ko'); end
        
                % 标 VPG 点
                if ~isnan(ut), plot(ut, diff_ppgSegs(ut), 'ro'); end
                if ~isnan(vt), plot(vt, diff_ppgSegs(vt), 'go'); end
                if ~isnan(wt), plot(wt, diff_ppgSegs(wt), 'ko'); end
        
                % 标 APG 点
                if ~isnan(at), plot(at, diff2_ppgSegs(at), 'ro'); end
                if ~isnan(bt), plot(bt, diff2_ppgSegs(bt), 'go'); end
                if ~isnan(ct), plot(ct, diff2_ppgSegs(ct), 'ko'); end
                if ~isnan(dt), plot(dt, diff2_ppgSegs(dt), 'mo'); end
                if ~isnan(et), plot(et, diff2_ppgSegs(et), 'co'); end
                if ~isnan(ft), plot(ft, diff2_ppgSegs(ft), 'yo'); end
        
                title(sprintf('PPG + VPG + APG | 文件名: %s', fname));
                xlabel('Time'); ylabel('Amplitude');
                legend('PPG', 'VPG', 'APG'); grid on;
         end
    end
    %% 基于PWA方法的特征计算
    % 计算每个PPG片段的长度
    % 计算每个PPG片段的特征
    N = 10;  % 或 N = length(ppgSeg);
    PWA_PPG_T1 = zeros(1, N);        % 峰值点索引
    PWA_PPG_T2 = zeros(1, N);        % DNT索引
    PWA_PPG_T3 = zeros(1, N);        % DPT索引
    PWA_PPG_T4 = zeros(1, N);        % 单个PPG周期长度
    PWA_PPG_T3_T1 = zeros(1, N);     % DPT - SPT
    PWA_PPG_AMP_OS = zeros(1, N);    % 收缩峰幅值
    PWA_PPG_AMP_ON = zeros(1, N);    % 切迹点幅值
    PWA_PPG_AMP_OD = zeros(1, N);    % 舒张谷幅值
    PWA_PPG_HR = zeros(1, N);        % 心率估计
    PWA_PPG_CT = zeros(1, N);        % 上升时间
    PWA_PPG_CS = zeros(1, N);        % 上升斜率
    PWA_PPG_ASI = zeros(1, N);       % 动脉硬化指数
    PWA_PPG_RI = zeros(1, N);        % 反射指数
    
    PWA_APG_c_a        = zeros(1, N);  % 比值c/a：动脉僵硬
    PWA_APG_e_a        = zeros(1, N);  % 比值e/a：动脉僵硬
    PWA_APG_b_a        = zeros(1, N);  % 比值b/a：反映硬度增加，随年龄增长增加
    PWA_APG_d_a        = zeros(1, N);  % 比值d/a：与动脉硬度降低有关
    PWA_APG_bcde_a     = zeros(1, N);  % (b-c-d-e)/a：血管衰老与动脉硬化疾病指数
    PWA_APG_be_a       = zeros(1, N);  % (b-e)/a：APG老化指数
    PWA_APG_cdb_a     = zeros(1, N);  % (c+d-b)/a：综合老化指数




for k = 1:N
    s = ppgSeg(k).s;
    len_s = length(s);

% PPG信号特征计算
        PWA_PPG_T1(k) = SPT(k);
        PWA_PPG_AMP_OS(k) = s(SPT(k))-s(1);
        PWA_PPG_CT(k) = SPT(k);  % 上升时间（点数）


    % DNT 检查，若为 NaN 可不设置幅值
    if ~isnan(DNT(k))
        PWA_PPG_T2(k) = DNT(k);
        PWA_PPG_AMP_ON(k) = s(DNT(k))-s(1);
    else
        PWA_PPG_T2(k) = NaN;
        PWA_PPG_AMP_ON(k) = NaN;
    end

    % DPT 检查，若为 NaN 不进行相关计算
    if ~isnan(DPT(k))
        PWA_PPG_T3(k) = DPT(k);
        PWA_PPG_AMP_OD(k) = s(DPT(k))-s(1);
        PWA_PPG_T3_T1(k) = PWA_PPG_T3(k) - PWA_PPG_T1(k);
    else
        PWA_PPG_T3(k) = NaN;
        PWA_PPG_T3_T1(k) = NaN;
        PWA_PPG_AMP_OD(k) = NaN;
    end

    % 总长度
    PWA_PPG_T4(k) = len_s;

    % 心率（粗略估计，60秒 / 每拍点数）
    if len_s > 0
        PWA_PPG_HR(k) = 60 / (len_s / fs);  % 注意单位一致性，len_s是点数，除以采样率得到秒
    end

    % 上升斜率（振幅 / 时间）
       PWA_PPG_CS(k) = PWA_PPG_AMP_OS(k) / PWA_PPG_CT(k);


    % 动脉硬化指数 ASI
    if ~isnan(PWA_PPG_T3_T1(k)) && ~isnan(PWA_PPG_T1(k))
        PWA_PPG_ASI(k) = PWA_PPG_T3_T1(k) / PWA_PPG_T1(k);
    else
        PWA_PPG_ASI(k) = NaN;
    end

    % 反射指数 RI（收缩峰 / 舒张谷）
    if ~isnan(PWA_PPG_AMP_OD(k))
        PWA_PPG_RI(k) = PWA_PPG_AMP_OS(k) / PWA_PPG_AMP_OD(k);
    else
        PWA_PPG_RI(k) = NaN;
    end

    % APG信号特征计算
    % 确保参考值a存在，避免除零
    if ~isnan(A(k))
        % c/a：表示动脉僵硬（血管弹性差）
        if ~isnan(C(k))
            PWA_APG_c_a(k) = C(k) / A(k);
        else
            PWA_APG_c_a(k) = nan;
        end
        % e/a：表示动脉僵硬（晚期反射）
        if ~isnan(E(k))
            PWA_APG_e_a(k) = E(k) / A(k);
        else
            PWA_APG_e_a(k) = nan;
        end
        % b/a：随动脉硬化程度增加而上升，常随年龄增长
        if ~isnan(B(k))
            PWA_APG_b_a(k) = B(k) / A(k);
        else
            PWA_APG_b_a(k) = nan;
        end
        % d/a：与动脉硬度下降相关，反映舒张功能改善
        if ~isnan(D(k))
            PWA_APG_d_a(k) = D(k) / A(k);
        else
            PWA_APG_d_a(k) = nan;
        end
        % (b-c-d-e)/a：血管衰老与动脉硬化相关疾病指标
        if ~any(isnan([B(k), C(k), D(k), E(k)]))
            PWA_APG_bcde_a(k) = (B(k) - C(k) - D(k) - E(k)) / A(k);
        else
            PWA_APG_bcde_a(k) = nan;
        end
        % (b-e)/a：APG老化指数，值越高表示反射波提前出现
        if ~any(isnan([B(k), E(k)]))
            PWA_APG_be_a(k) = (B(k) - E(k)) / A(k);
        else
            PWA_APG_be_a(k) = nan;
        end
        % (c+d-b)/a：综合评价老化程度
        if ~any(isnan([C(k), D(k), B(k)]))
            PWA_APG_cdb_a(k) = (C(k) + D(k) - B(k)) / A(k);
        else
            PWA_APG_cdb_a(k) =nan;
        end
    end
   
end
    % ==== 中位数特征提取 ====
    %PPG原始信号特征
    PWA_PPG_T1_median        = median(PWA_PPG_T1, 'omitnan');
    PWA_PPG_T2_median        = median(PWA_PPG_T2, 'omitnan');
    PWA_PPG_T3_median        = median(PWA_PPG_T3, 'omitnan');
    PWA_PPG_T4_median        = median(PWA_PPG_T4, 'omitnan');
    PWA_PPG_T3_T1_median     = median(PWA_PPG_T3_T1, 'omitnan');
    PWA_PPG_AMP_OS_median    = median(PWA_PPG_AMP_OS, 'omitnan');
    PWA_PPG_AMP_ON_median    = median(PWA_PPG_AMP_ON, 'omitnan');
    PWA_PPG_AMP_OD_median    = median(PWA_PPG_AMP_OD, 'omitnan');  
    PWA_PPG_HR_median        = median(PWA_PPG_HR, 'omitnan');
    PWA_PPG_CT_median        = median(PWA_PPG_CT, 'omitnan');
    PWA_PPG_CS_median        = median(PWA_PPG_CS, 'omitnan');
    PWA_PPG_ASI_median       = median(PWA_PPG_ASI, 'omitnan');
    PWA_PPG_RI_median        = median(PWA_PPG_RI, 'omitnan');

    %VPG信号特征


    %APG信号特征
    PWA_APG_c_a_median    = median(PWA_APG_c_a, 'omitnan');
    PWA_APG_e_a_median    = median(PWA_APG_e_a, 'omitnan');
    PWA_APG_b_a_median    = median(PWA_APG_b_a, 'omitnan');
    PWA_APG_d_a_median    = median(PWA_APG_d_a, 'omitnan');
    PWA_APG_bcde_a_median = median(PWA_APG_bcde_a, 'omitnan');
    PWA_APG_be_a_median   = median(PWA_APG_be_a, 'omitnan');
    PWA_APG_cdb_a_median  = median(PWA_APG_cdb_a, 'omitnan');

    %特征点保存
    % 前提：变量 SP、DP、DN、SPT、DPT、DNT、U、V、W、UT、VT、WT、A、B、C、D、E、F、AT、BT、CT、DT、ET、FT
    
    % === 基于原始PPG信号 ===
    SP_median   = median(SP(1:N),   'omitnan');   % 收缩峰值
    DP_median   = median(DP(1:N),   'omitnan');   % 舒张谷值
    DN_median   = median(DN(1:N),   'omitnan');   % 切迹点值
    
    SPT_median  = median(SPT(1:N),  'omitnan');   % 收缩峰时间
    DPT_median  = median(DPT(1:N),  'omitnan');   % 舒张谷时间
    DNT_median  = median(DNT(1:N),  'omitnan');   % 切迹点时间
    
    % === 基于一阶导（VPG） ===
    U_median    = median(U(1:N),    'omitnan');   % U 点值
    V_median    = median(V(1:N),    'omitnan');   % V 点值
    W_median    = median(W(1:N),    'omitnan');   % W 点值
    
    UT_median   = median(UT(1:N),   'omitnan');   % U 时间
    VT_median   = median(VT(1:N),   'omitnan');   % V 时间
    WT_median   = median(WT(1:N),   'omitnan');   % W 时间
    
    % === 基于二阶导（APG） ===
    A_median    = median(A(1:N),    'omitnan');   % A 点值
    B_median    = median(B(1:N),    'omitnan');   % B 点值
    C_median    = median(C(1:N),    'omitnan');   % C 点值
    D_median    = median(D(1:N),    'omitnan');   % D 点值
    E_median    = median(E(1:N),    'omitnan');   % E 点值
    F_median    = median(F(1:N),    'omitnan');   % F 点值
    
    AT_median   = median(AT(1:N),   'omitnan');   % A 时间
    BT_median   = median(BT(1:N),   'omitnan');   % B 时间
    CT_median   = median(CT(1:N),   'omitnan');   % C 时间
    DT_median   = median(DT(1:N),   'omitnan');   % D 时间
    ET_median   = median(ET(1:N),   'omitnan');   % E 时间
    FT_median   = median(FT(1:N),   'omitnan');   % F 时间


        %% 假设已经有两个数组：valid_OS 和 PWA_PPG_T4_median
        % 步骤 1：去除 NaN
    valid_OS = PWA_PPG_AMP_OS(~isnan(PWA_PPG_AMP_OS));
    % 血管变异性（VV）指标计算
    PPG_FD_VV_mean = mean(valid_OS);                 % 平均值
    PPG_FD_VV_std = std(valid_OS);                   % 标准差
    PPG_FD_VV_cv = PPG_FD_VV_std / PPG_FD_VV_mean;                 % 变异系数
    
    % Poincare 分析
    PPG_FD_VV_diff = diff(valid_OS);
    PPG_FD_VV_SD1 = std(PPG_FD_VV_diff)/sqrt(2);            % 短期波动
    PPG_FD_VV_SD2 = sqrt(2*PPG_FD_VV_std^2 - PPG_FD_VV_SD1^2);    % 长期波动
    
    % fprintf('血管变异性指标：\n');
    % fprintf('均值: %.3f, 标准差: %.3f, 变异系数: %.3f\n', PPG_FD_VV_mean, PPG_FD_VV_std, PPG_FD_VV_cv);
    % fprintf('Poincare SD1: %.3f, SD2: %.3f\n\n', PPG_FD_VV_SD1, PPG_FD_VV_SD2);
    
    %% 心率变异性（HRV）指标计算
    RR = PWA_PPG_T4;  % 假设单位为 ms
    
    PPG_FD_HR_mean = mean(RR);
    PPG_FD_HR_std = std(RR);        % SDNN
    PPG_FD_HR_diff = diff(RR);
    PPG_FD_RMSSD = sqrt(mean(PPG_FD_HR_diff.^2));
    PPG_FD_NN50 = sum(abs(PPG_FD_HR_diff) > 50);  % 以50ms为阈值
    PPG_FD_pNN50 = PPG_FD_NN50 / length(PPG_FD_HR_diff) * 100;
    
    % Poincare 分析
    PPG_FD_HR_SD1 = std(PPG_FD_HR_diff)/sqrt(2);
    PPG_FD_HR_SD2 = sqrt(2*PPG_FD_HR_std^2 - PPG_FD_HR_SD1^2);
    
    % fprintf('心率变异性指标：\n');
    % fprintf('均值RR: %.3f ms, SDNN: %.3f ms\n', PPG_FD_HR_mean, PPG_FD_HR_std);
    % fprintf('RMSSD: %.3f ms, NN50: %d, pNN50: %.2f%%\n', PPG_FD_RMSSD, PPG_FD_NN50, PPG_FD_pNN50);
    % fprintf('Poincare SD1: %.3f ms, SD2: %.3f ms\n', PPG_FD_HR_SD1, PPG_FD_HR_SD2);

    % 步骤 1：去除 NaN
    % valid_OS = PWA_PPG_AMP_OS(~isnan(PWA_PPG_AMP_OS));
    
    % 步骤 2：归一化为概率分布（用直方图近似）
    [counts, edges] = histcounts(valid_OS, 'Normalization', 'probability');
    p = counts(counts > 0);  % 删除为0的项以避免 log2(0)
    
    % 步骤 3：计算信息熵（单位：比特）
    PWA_OS_entropy = -sum(p .* log2(p));
    % fprintf('PWA_AMP_OS 熵值（信息熵）: %.4f bits\n', PWA_OS_entropy);
    
    amp_mean = mean(valid_OS);
    amp_min = min(valid_OS);
    amp_max = max(valid_OS);
    
    if amp_max ~= amp_min
        PWA_VV = (amp_mean - amp_min) / (amp_max - amp_min);
    else
        PWA_VV = NaN;  % 避免除以0
    end
    
    % fprintf('PWA_AMP_OS VV值: %.4f\n', PWA_VV);




    
%% 基于微循环的K值计算
    [MD_K_up, MD_K_down, MD_K_total] = Features_extraction_MD_K(ppgSeg, SPT);
     MD_K_up_median = median(MD_K_up);  % KK的中位数
    MD_K_down_median = median(MD_K_down);  % KK的中位数
    MD_K_total_median = median(MD_K_total);  % KK的中位数
    if(1)
        fprintf('文件名: %s\n', fname);
        fprintf('MD_K_up 中位数: %.4f\n', MD_K_up_median);
        fprintf('MD_K_down 中位数: %.4f\n', MD_K_down_median);
        fprintf('MD_K_total 中位数: %.4f\n', MD_K_total_median);
    end

%% 归一化
    ppgSeg_norm = normalize_ppgSeg(ppgSeg);
%% 基于容积变化的波形宽度和波形面积特征计算
    [VVA_T_up_med, VVA_T_down_med, VVA_width_med, ...
      VVA_width_ratio_LR_med, VVA_width_ratio_med, ...
      VVA_amplitude_med, VVA_width_triangle_area_med, ...
      VVA_width_triangle_area_ratio_med, ...
      VVA_area_up_real_med, VVA_area_down_real_med, ...
      VVA_area_ratio_LR_med, VVA_area_ratio_up_med, ...
      VVA_area_ratio_down_med, VVA_area_total_med, ...
      VVA_area_ratio_total_med] = Feature_extraction_VVA(ppgSeg_norm, SPT);


%% 数据保存
        % 构造保存用数组（写入 Excel）

feature_data = {fname, ...
    SP_median, SPT_median, DP_median, DPT_median, DN_median, DNT_median,U_median, UT_median, V_median, VT_median, W_median,  WT_median, ...
    A_median,  AT_median,B_median,  BT_median,C_median,  CT_median, D_median,  DT_median, E_median,  ET_median, F_median,  FT_median ...
    MD_K_up_median, MD_K_down_median, MD_K_total_median, ...
    PWA_PPG_T1_median, PWA_PPG_T2_median, PWA_PPG_T3_median, PWA_PPG_T4_median, PWA_PPG_T3_T1_median, ...
    PWA_PPG_AMP_OS_median, PWA_PPG_AMP_ON_median, PWA_PPG_AMP_OD_median, ...
    PWA_PPG_HR_median, PWA_PPG_CT_median, PWA_PPG_CS_median, PWA_PPG_ASI_median, PWA_PPG_RI_median, ...
    PWA_APG_c_a_median,PWA_APG_e_a_median,PWA_APG_b_a_median,PWA_APG_d_a_median,PWA_APG_bcde_a_median,PWA_APG_be_a_median,PWA_APG_cdb_a_median, ...
    VVA_T_up_med(1), VVA_T_up_med(2), VVA_T_up_med(3), VVA_T_up_med(4), VVA_T_up_med(5), ...
    VVA_T_up_med(6), VVA_T_up_med(7), VVA_T_up_med(8), VVA_T_up_med(9), VVA_T_up_med(10), ...
    VVA_T_up_med(11), VVA_T_up_med(12), VVA_T_up_med(13), VVA_T_up_med(14), VVA_T_up_med(15), ...
    VVA_T_up_med(16), VVA_T_up_med(17), VVA_T_up_med(18), VVA_T_up_med(19), ...
    VVA_T_down_med(1), VVA_T_down_med(2), VVA_T_down_med(3), VVA_T_down_med(4), VVA_T_down_med(5), ...
    VVA_T_down_med(6), VVA_T_down_med(7), VVA_T_down_med(8), VVA_T_down_med(9), VVA_T_down_med(10), ...
    VVA_T_down_med(11), VVA_T_down_med(12), VVA_T_down_med(13), VVA_T_down_med(14), VVA_T_down_med(15), ...
    VVA_T_down_med(16), VVA_T_down_med(17), VVA_T_down_med(18), VVA_T_down_med(19), ...
    VVA_width_med(1), VVA_width_med(2), VVA_width_med(3), VVA_width_med(4), VVA_width_med(5), ...
    VVA_width_med(6), VVA_width_med(7), VVA_width_med(8), VVA_width_med(9), VVA_width_med(10), ...
    VVA_width_med(11), VVA_width_med(12), VVA_width_med(13), VVA_width_med(14), VVA_width_med(15), ...
    VVA_width_med(16), VVA_width_med(17), VVA_width_med(18), VVA_width_med(19), ...
    VVA_width_ratio_LR_med(1), VVA_width_ratio_LR_med(2), VVA_width_ratio_LR_med(3), VVA_width_ratio_LR_med(4), VVA_width_ratio_LR_med(5), ...
    VVA_width_ratio_LR_med(6), VVA_width_ratio_LR_med(7), VVA_width_ratio_LR_med(8), VVA_width_ratio_LR_med(9), VVA_width_ratio_LR_med(10), ...
    VVA_width_ratio_LR_med(11), VVA_width_ratio_LR_med(12), VVA_width_ratio_LR_med(13), VVA_width_ratio_LR_med(14), VVA_width_ratio_LR_med(15), ...
    VVA_width_ratio_LR_med(16), VVA_width_ratio_LR_med(17), VVA_width_ratio_LR_med(18), VVA_width_ratio_LR_med(19), ...
    VVA_width_ratio_med(1), VVA_width_ratio_med(2), VVA_width_ratio_med(3), VVA_width_ratio_med(4), VVA_width_ratio_med(5), ...
    VVA_width_ratio_med(6), VVA_width_ratio_med(7), VVA_width_ratio_med(8), VVA_width_ratio_med(9), VVA_width_ratio_med(10), ...
    VVA_width_ratio_med(11), VVA_width_ratio_med(12), VVA_width_ratio_med(13), VVA_width_ratio_med(14), VVA_width_ratio_med(15), ...
    VVA_width_ratio_med(16), VVA_width_ratio_med(17), VVA_width_ratio_med(18), VVA_width_ratio_med(19), ...
    VVA_amplitude_med(1), VVA_amplitude_med(2), VVA_amplitude_med(3), VVA_amplitude_med(4), VVA_amplitude_med(5), ...
    VVA_amplitude_med(6), VVA_amplitude_med(7), VVA_amplitude_med(8), VVA_amplitude_med(9), VVA_amplitude_med(10), ...
    VVA_amplitude_med(11), VVA_amplitude_med(12), VVA_amplitude_med(13), VVA_amplitude_med(14), VVA_amplitude_med(15), ...
    VVA_amplitude_med(16), VVA_amplitude_med(17), VVA_amplitude_med(18), VVA_amplitude_med(19), ...
    VVA_width_triangle_area_med(1), VVA_width_triangle_area_med(2), VVA_width_triangle_area_med(3), VVA_width_triangle_area_med(4), VVA_width_triangle_area_med(5), ...
    VVA_width_triangle_area_med(6), VVA_width_triangle_area_med(7), VVA_width_triangle_area_med(8), VVA_width_triangle_area_med(9), VVA_width_triangle_area_med(10), ...
    VVA_width_triangle_area_med(11), VVA_width_triangle_area_med(12), VVA_width_triangle_area_med(13), VVA_width_triangle_area_med(14), VVA_width_triangle_area_med(15), ...
    VVA_width_triangle_area_med(16), VVA_width_triangle_area_med(17), VVA_width_triangle_area_med(18), VVA_width_triangle_area_med(19), ...
    VVA_width_triangle_area_ratio_med(1), VVA_width_triangle_area_ratio_med(2), VVA_width_triangle_area_ratio_med(3), VVA_width_triangle_area_ratio_med(4), VVA_width_triangle_area_ratio_med(5), ...
    VVA_width_triangle_area_ratio_med(6), VVA_width_triangle_area_ratio_med(7), VVA_width_triangle_area_ratio_med(8), VVA_width_triangle_area_ratio_med(9), VVA_width_triangle_area_ratio_med(10), ...
    VVA_width_triangle_area_ratio_med(11), VVA_width_triangle_area_ratio_med(12), VVA_width_triangle_area_ratio_med(13), VVA_width_triangle_area_ratio_med(14), VVA_width_triangle_area_ratio_med(15), ...
    VVA_width_triangle_area_ratio_med(16), VVA_width_triangle_area_ratio_med(17), VVA_width_triangle_area_ratio_med(18), VVA_width_triangle_area_ratio_med(19), ...
    VVA_area_up_real_med(1), VVA_area_up_real_med(2), VVA_area_up_real_med(3), VVA_area_up_real_med(4), VVA_area_up_real_med(5), ...
    VVA_area_up_real_med(6), VVA_area_up_real_med(7), VVA_area_up_real_med(8), VVA_area_up_real_med(9), VVA_area_up_real_med(10), ...
    VVA_area_up_real_med(11), VVA_area_up_real_med(12), VVA_area_up_real_med(13), VVA_area_up_real_med(14), VVA_area_up_real_med(15), ...
    VVA_area_up_real_med(16), VVA_area_up_real_med(17), VVA_area_up_real_med(18), VVA_area_up_real_med(19), ...
    VVA_area_down_real_med(1), VVA_area_down_real_med(2), VVA_area_down_real_med(3), VVA_area_down_real_med(4), VVA_area_down_real_med(5), ...
    VVA_area_down_real_med(6), VVA_area_down_real_med(7), VVA_area_down_real_med(8), VVA_area_down_real_med(9), VVA_area_down_real_med(10), ...
    VVA_area_down_real_med(11), VVA_area_down_real_med(12), VVA_area_down_real_med(13), VVA_area_down_real_med(14), VVA_area_down_real_med(15), ...
    VVA_area_down_real_med(16), VVA_area_down_real_med(17), VVA_area_down_real_med(18), VVA_area_down_real_med(19), ...
    VVA_area_ratio_LR_med(1), VVA_area_ratio_LR_med(2), VVA_area_ratio_LR_med(3), VVA_area_ratio_LR_med(4), VVA_area_ratio_LR_med(5), ...
    VVA_area_ratio_LR_med(6), VVA_area_ratio_LR_med(7), VVA_area_ratio_LR_med(8), VVA_area_ratio_LR_med(9), VVA_area_ratio_LR_med(10), ...
    VVA_area_ratio_LR_med(11), VVA_area_ratio_LR_med(12), VVA_area_ratio_LR_med(13), VVA_area_ratio_LR_med(14), VVA_area_ratio_LR_med(15), ...
    VVA_area_ratio_LR_med(16), VVA_area_ratio_LR_med(17), VVA_area_ratio_LR_med(18), VVA_area_ratio_LR_med(19), ...
    VVA_area_ratio_up_med(1), VVA_area_ratio_up_med(2), VVA_area_ratio_up_med(3), VVA_area_ratio_up_med(4), VVA_area_ratio_up_med(5), ...
    VVA_area_ratio_up_med(6), VVA_area_ratio_up_med(7), VVA_area_ratio_up_med(8), VVA_area_ratio_up_med(9), VVA_area_ratio_up_med(10), ...
    VVA_area_ratio_up_med(11), VVA_area_ratio_up_med(12), VVA_area_ratio_up_med(13), VVA_area_ratio_up_med(14), VVA_area_ratio_up_med(15), ...
    VVA_area_ratio_up_med(16), VVA_area_ratio_up_med(17), VVA_area_ratio_up_med(18), VVA_area_ratio_up_med(19), ...
    VVA_area_ratio_down_med(1), VVA_area_ratio_down_med(2), VVA_area_ratio_down_med(3), VVA_area_ratio_down_med(4), VVA_area_ratio_down_med(5), ...
    VVA_area_ratio_down_med(6), VVA_area_ratio_down_med(7), VVA_area_ratio_down_med(8), VVA_area_ratio_down_med(9), VVA_area_ratio_down_med(10), ...
    VVA_area_ratio_down_med(11), VVA_area_ratio_down_med(12), VVA_area_ratio_down_med(13), VVA_area_ratio_down_med(14),VVA_area_ratio_down_med(15), ...
    VVA_area_ratio_down_med(16), VVA_area_ratio_down_med(17), VVA_area_ratio_down_med(18), VVA_area_ratio_down_med(19), ...
    VVA_area_total_med(1), VVA_area_total_med(2), VVA_area_total_med(3), VVA_area_total_med(4), VVA_area_total_med(5), ...
    VVA_area_total_med(6), VVA_area_total_med(7), VVA_area_total_med(8), VVA_area_total_med(9), VVA_area_total_med(10), ...
    VVA_area_total_med(11), VVA_area_total_med(12), VVA_area_total_med(13), VVA_area_total_med(14),VVA_area_total_med(15), ...
    VVA_area_total_med(16), VVA_area_total_med(17), VVA_area_total_med(18), VVA_area_total_med(19), ...
    VVA_area_ratio_total_med(1), VVA_area_ratio_total_med(2), VVA_area_ratio_total_med(3), VVA_area_ratio_total_med(4), VVA_area_ratio_total_med(5), ...
    VVA_area_ratio_total_med(6), VVA_area_ratio_total_med(7), VVA_area_ratio_total_med(8), VVA_area_ratio_total_med(9), VVA_area_ratio_total_med(10), ...
    VVA_area_ratio_total_med(11), VVA_area_ratio_total_med(12), VVA_area_ratio_total_med(13), VVA_area_ratio_total_med(14),VVA_area_ratio_total_med(15), ...
    VVA_area_ratio_total_med(16), VVA_area_ratio_total_med(17), VVA_area_ratio_total_med(18), VVA_area_ratio_total_med(19), ...
    PWA_OS_entropy, PWA_VV,PPG_FD_VV_mean, PPG_FD_VV_std, PPG_FD_VV_cv,PPG_FD_VV_SD1, PPG_FD_VV_SD2, PPG_FD_HR_mean, PPG_FD_HR_std, PPG_FD_RMSSD, PPG_FD_NN50, PPG_FD_pNN50, PPG_FD_HR_SD1, PPG_FD_HR_SD2};    
       

header = {'Filename', ...
    'SP_median', 'SPT_median','DP_median','DPT_median','DN_median','DNT_median','U_median', 'UT_median','V_median',  'VT_median', 'W_median', 'WT_median', ...
    'A_median',  'AT_median', 'B_median',  'BT_median', 'C_median',  'CT_median','D_median',  'DT_median','E_median',  'ET_median', 'F_median',  'FT_median' ...
    'MD_K_up_median', 'MD_K_down_median', 'MD_K_total_median', ...
    'PWA_T1_median', 'PWA_T2_median', 'PWA_T3_median', 'PWA_T4_median', 'PWA_T3_T1_median', ...
    'PWA_AMP_OS_median', 'PWA_AMP_ON_median', 'PWA_AMP_OD_median', ...
    'PWA_HR_median', 'PWA_CT_median', 'PWA_CS_median', 'PWA_ASI_median', 'PWA_RI_median', ...
    'PWA_APG_c_a_median','PWA_APG_e_a_median','PWA_APG_b_a_median','PWA_APG_d_a_median','PWA_APG_bcde_a_median','PWA_APG_be_a_median','PWA_APG_cdb_a_median', ...
    'VVA_T_up_med_5', 'VVA_T_up_med_10', 'VVA_T_up_med_15', 'VVA_T_up_med_20', 'VVA_T_up_med_25', ...
    'VVA_T_up_med_30', 'VVA_T_up_med_35', 'VVA_T_up_med_40', 'VVA_T_up_med_45', 'VVA_T_up_med_50', ...
    'VVA_T_up_med_55', 'VVA_T_up_med_60', 'VVA_T_up_med_65', 'VVA_T_up_med_70', 'VVA_T_up_med_75', ...
    'VVA_T_up_med_80', 'VVA_T_up_med_85', 'VVA_T_up_med_90', 'VVA_T_up_med_95', ...
    'VVA_T_down_med_5', 'VVA_T_down_med_10', 'VVA_T_down_med_15', 'VVA_T_down_med_20', 'VVA_T_down_med_25', ...
    'VVA_T_down_med_30', 'VVA_T_down_med_35', 'VVA_T_down_med_40', 'VVA_T_down_med_45', 'VVA_T_down_med_50', ...
    'VVA_T_down_med_55', 'VVA_T_down_med_60', 'VVA_T_down_med_65', 'VVA_T_down_med_70', 'VVA_T_down_med_75', ...
    'VVA_T_down_med_80', 'VVA_T_down_med_85', 'VVA_T_down_med_90', 'VVA_T_down_med_95', ...
    'VVA_width_med_5', 'VVA_width_med_10', 'VVA_width_med_15', 'VVA_width_med_20', 'VVA_width_med_25', ...
    'VVA_width_med_30', 'VVA_width_med_35', 'VVA_width_med_40', 'VVA_width_med_45', 'VVA_width_med_50', ...
    'VVA_width_med_55', 'VVA_width_med_60', 'VVA_width_med_65', 'VVA_width_med_70', 'VVA_width_med_75', ...
    'VVA_width_med_80', 'VVA_width_med_85', 'VVA_width_med_90', 'VVA_width_med_95', ...
    'VVA_width_ratio_LR_med_5', 'VVA_width_ratio_LR_med_10', 'VVA_width_ratio_LR_med_15', 'VVA_width_ratio_LR_med_20', 'VVA_width_ratio_LR_med_25', ...
    'VVA_width_ratio_LR_med_30', 'VVA_width_ratio_LR_med_35', 'VVA_width_ratio_LR_med_40', 'VVA_width_ratio_LR_med_45', 'VVA_width_ratio_LR_med_50', ...
    'VVA_width_ratio_LR_med_55', 'VVA_width_ratio_LR_med_60', 'VVA_width_ratio_LR_med_65', 'VVA_width_ratio_LR_med_70', 'VVA_width_ratio_LR_med_75', ...
    'VVA_width_ratio_LR_med_80', 'VVA_width_ratio_LR_med_85', 'VVA_width_ratio_LR_med_90', 'VVA_width_ratio_LR_med_95', ...
    'VVA_width_ratio_med_5', 'VVA_width_ratio_med_10', 'VVA_width_ratio_med_15', 'VVA_width_ratio_med_20', 'VVA_width_ratio_med_25', ...
    'VVA_width_ratio_med_30', 'VVA_width_ratio_med_35', 'VVA_width_ratio_med_40', 'VVA_width_ratio_med_45', 'VVA_width_ratio_med_50', ...
    'VVA_width_ratio_med_55', 'VVA_width_ratio_med_60', 'VVA_width_ratio_med_65', 'VVA_width_ratio_med_70', 'VVA_width_ratio_med_75', ...
    'VVA_width_ratio_med_80', 'VVA_width_ratio_med_85', 'VVA_width_ratio_med_90', 'VVA_width_ratio_med_95', ...
    'VVA_amplitude_med_5', 'VVA_amplitude_med_10', 'VVA_amplitude_med_15', 'VVA_amplitude_med_20', 'VVA_amplitude_med_25', ...
    'VVA_amplitude_med_30', 'VVA_amplitude_med_35', 'VVA_amplitude_med_40', 'VVA_amplitude_med_45', 'VVA_amplitude_med_50', ...
    'VVA_amplitude_med_55', 'VVA_amplitude_med_60', 'VVA_amplitude_med_65', 'VVA_amplitude_med_70', 'VVA_amplitude_med_75', ...
    'VVA_amplitude_med_80', 'VVA_amplitude_med_85', 'VVA_amplitude_med_90', 'VVA_amplitude_med_95', ...
    'VVA_width_triangle_area_med_5', 'VVA_width_triangle_area_med_10', 'VVA_width_triangle_area_med_15', 'VVA_width_triangle_area_med_20', 'VVA_width_triangle_area_med_25', ...
    'VVA_width_triangle_area_med_30', 'VVA_width_triangle_area_med_35', 'VVA_width_triangle_area_med_40', 'VVA_width_triangle_area_med_45', 'VVA_width_triangle_area_med_50', ...
    'VVA_width_triangle_area_med_55', 'VVA_width_triangle_area_med_60', 'VVA_width_triangle_area_med_65', 'VVA_width_triangle_area_med_70', 'VVA_width_triangle_area_med_75', ...
    'VVA_width_triangle_area_med_80', 'VVA_width_triangle_area_med_85', 'VVA_width_triangle_area_med_90', 'VVA_width_triangle_area_med_95', ...
    'VVA_width_triangle_area_ratio_med_5', 'VVA_width_triangle_area_ratio_med_10', 'VVA_width_triangle_area_ratio_med_15', 'VVA_width_triangle_area_ratio_med_20', 'VVA_width_triangle_area_ratio_med_25', ...
    'VVA_width_triangle_area_ratio_med_30', 'VVA_width_triangle_area_ratio_med_35', 'VVA_width_triangle_area_ratio_med_40', 'VVA_width_triangle_area_ratio_med_45', 'VVA_width_triangle_area_ratio_med_50', ...
    'VVA_width_triangle_area_ratio_med_55', 'VVA_width_triangle_area_ratio_med_60', 'VVA_width_triangle_area_ratio_med_65', 'VVA_width_triangle_area_ratio_med_70', 'VVA_width_triangle_area_ratio_med_75', ...
    'VVA_width_triangle_area_ratio_med_80', 'VVA_width_triangle_area_ratio_med_85', 'VVA_width_triangle_area_ratio_med_90', 'VVA_width_triangle_area_ratio_med_95', ...
    'VVA_area_up_real_med_5', 'VVA_area_up_real_med_10', 'VVA_area_up_real_med_15', 'VVA_area_up_real_med_20', 'VVA_area_up_real_med_25', ...
    'VVA_area_up_real_med_30', 'VVA_area_up_real_med_35', 'VVA_area_up_real_med_40', 'VVA_area_up_real_med_45', 'VVA_area_up_real_med_50', ...
    'VVA_area_up_real_med_55', 'VVA_area_up_real_med_60', 'VVA_area_up_real_med_65', 'VVA_area_up_real_med_70', 'VVA_area_up_real_med_75', ...
    'VVA_area_up_real_med_80', 'VVA_area_up_real_med_85', 'VVA_area_up_real_med_90', 'VVA_area_up_real_med_95', ...
    'VVA_area_down_real_med_5', 'VVA_area_down_real_med_10', 'VVA_area_down_real_med_15', 'VVA_area_down_real_med_20', 'VVA_area_down_real_med_25', ...
    'VVA_area_down_real_med_30', 'VVA_area_down_real_med_35', 'VVA_area_down_real_med_40', 'VVA_area_down_real_med_45', 'VVA_area_down_real_med_50', ...
    'VVA_area_down_real_med_55', 'VVA_area_down_real_med_60', 'VVA_area_down_real_med_65', 'VVA_area_down_real_med_70', 'VVA_area_down_real_med_75', ...
    'VVA_area_down_real_med_80', 'VVA_area_down_real_med_85', 'VVA_area_down_real_med_90', 'VVA_area_down_real_med_95', ...
    'VVA_area_ratio_LR_med_5', 'VVA_area_ratio_LR_med_10', 'VVA_area_ratio_LR_med_15', 'VVA_area_ratio_LR_med_20', 'VVA_area_ratio_LR_med_25', ...
    'VVA_area_ratio_LR_med_30', 'VVA_area_ratio_LR_med_35', 'VVA_area_ratio_LR_med_40', 'VVA_area_ratio_LR_med_45', 'VVA_area_ratio_LR_med_50', ...
    'VVA_area_ratio_LR_med_55', 'VVA_area_ratio_LR_med_60', 'VVA_area_ratio_LR_med_65', 'VVA_area_ratio_LR_med_70', 'VVA_area_ratio_LR_med_75', ...
    'VVA_area_ratio_LR_med_80', 'VVA_area_ratio_LR_med_85', 'VVA_area_ratio_LR_med_90', 'VVA_area_ratio_LR_med_95', ...
    'VVA_area_ratio_up_med_5', 'VVA_area_ratio_up_med_10', 'VVA_area_ratio_up_med_15', 'VVA_area_ratio_up_med_20', 'VVA_area_ratio_up_med_25', ...
    'VVA_area_ratio_up_med_30', 'VVA_area_ratio_up_med_35', 'VVA_area_ratio_up_med_40', 'VVA_area_ratio_up_med_45', 'VVA_area_ratio_up_med_50', ...
    'VVA_area_ratio_up_med_55', 'VVA_area_ratio_up_med_60', 'VVA_area_ratio_up_med_65', 'VVA_area_ratio_up_med_70', 'VVA_area_ratio_up_med_75', ...
    'VVA_area_ratio_up_med_80', 'VVA_area_ratio_up_med_85', 'VVA_area_ratio_up_med_90', 'VVA_area_ratio_up_med_95', ...
    'VVA_area_ratio_down_med_5', 'VVA_area_ratio_down_med_10', 'VVA_area_ratio_down_med_15', 'VVA_area_ratio_down_med_20', 'VVA_area_ratio_down_med_25', ...
    'VVA_area_ratio_down_med_30', 'VVA_area_ratio_down_med_35', 'VVA_area_ratio_down_med_40', 'VVA_area_ratio_down_med_45', 'VVA_area_ratio_down_med_50', ...
    'VVA_area_ratio_down_med_55', 'VVA_area_ratio_down_med_60', 'VVA_area_ratio_down_med_65', 'VVA_area_ratio_down_med_70', 'VVA_area_ratio_down_med_75', ...
    'VVA_area_ratio_down_med_80', 'VVA_area_ratio_down_med_85', 'VVA_area_ratio_down_med_90', 'VVA_area_ratio_down_med_95', ...
    'VVA_area_total_med_5', 'VVA_area_total_med_10', 'VVA_area_total_med_15', 'VVA_area_total_med_20', 'VVA_area_total_med_25', ...
    'VVA_area_total_med_30', 'VVA_area_total_med_35', 'VVA_area_total_med_40', 'VVA_area_total_med_45', 'VVA_area_total_med_50', ...
    'VVA_area_total_med_55', 'VVA_area_total_med_60', 'VVA_area_total_med_65', 'VVA_area_total_med_70', 'VVA_area_total_med_75', ...
    'VVA_area_total_med_80', 'VVA_area_total_med_85', 'VVA_area_total_med_90', 'VVA_area_total_med_95', ...
    'VVA_area_ratio_total_med_5', 'VVA_area_ratio_total_med_10', 'VVA_area_ratio_total_med_15', 'VVA_area_ratio_total_med_20', 'VVA_area_ratio_total_med_25', ...
    'VVA_area_ratio_total_med_30', 'VVA_area_ratio_total_med_35', 'VVA_area_ratio_total_med_40', 'VVA_area_ratio_total_med_45', 'VVA_area_ratio_total_med_50', ...
    'VVA_area_ratio_total_med_55', 'VVA_area_ratio_total_med_60', 'VVA_area_ratio_total_med_65', 'VVA_area_ratio_total_med_70', 'VVA_area_ratio_total_med_75', ...
    'VVA_area_ratio_total_med_80', 'VVA_area_ratio_total_med_85', 'VVA_area_ratio_total_med_90', 'VVA_area_ratio_total_med_95', ...
    'PWA_OS_entropy', 'PWA_VV',"PPG_FD_VV_mean", "PPG_FD_VV_std", "PPG_FD_VV_cv", "PPG_FD_VV_SD1", "PPG_FD_VV_SD2", "PPG_FD_HR_mean", "PPG_FD_HR_std", "PPG_FD_RMSSD", "PPG_FD_NN50", "PPG_FD_pNN50", "PPG_FD_HR_SD1", "PPG_FD_HR_SD2"};




        % 如果文件不存在，创建并写入表头
        if ~exist(folder_path_save, 'file')
            writecell(header, folder_path_save, 'Sheet', 1, 'Range', 'A1');
        end
        if(data_column==1) % 写入右侧数据
        writecell(feature_data, folder_path_save, 'Sheet', 1, 'Range', ['A' num2str(row_num)]);                   
        % 写入左侧数据
        else
        % 写入左侧数据
         writecell(feature_data, folder_path_save, 'Sheet', 1, 'Range', ['A' num2str(row_num+1)]);
        end
        
        fprintf('已处理文件: %s \n', fname);
    end
    row_num = row_num + 1 ;  % 每次增加2行
end

