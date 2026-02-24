function diff2_signal = signal_diff2(diff_signal,signal_len)
     diff2_signal = zeros(signal_len,1);
     for j = 1+2 : signal_len-2
        diff2_signal(j)=  diff_signal(j+1)-diff_signal(j-1) + 2 *( diff_signal(j+2)-diff_signal(j-2)) ;
        diff2_signal(j) =diff2_signal(j)*3;
     end
     diff2_signal(1) = diff2_signal(3); 
     diff2_signal(2) = diff2_signal(3);
     diff2_signal(signal_len-1)=diff2_signal(signal_len-2);
     diff2_signal(signal_len)=diff2_signal(signal_len-2);
end

