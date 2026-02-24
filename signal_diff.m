function diff_signal = signal_diff(signal_hat,PPG_data_len)
%SIGNAL_ 此处显示有关此函数的摘要
%   此处显示详细说明
    diff_signal = zeros(PPG_data_len,1);
     for j = 1+2 : PPG_data_len-2
        diff_signal(j) =  signal_hat(j+1)-signal_hat(j-1) + 2 *( signal_hat(j+2)-signal_hat(j-2)) ;
        diff_signal(j) = diff_signal(j)*6;
     end
     diff_signal(1) = diff_signal(3); 
     diff_signal(2) = diff_signal(3);
     diff_signal(PPG_data_len-1)=diff_signal(PPG_data_len-2);
     diff_signal(PPG_data_len)=diff_signal(PPG_data_len-2);
end

