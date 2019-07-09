%% Initialization
N_test = 1:5;                                       % 1:10 orginal
%N_test = 12:16;
SVM_error = N_test*0;
mmse_error = SVM_error;
SNR_dB = -5;
ave_num = 10000;                                   % 100000 original
T = 300;


%% Simulation
for loop1 = 1:length(N_test)
    N = N_test(loop1);                              % Number of antennas
    b = 0:2^N-1;                                    % Transmitted signals
    B = (de2bi(b,'left-msb'))';                     % Convert to binary and transpose
    B = 2*B - 1;                                    % Convert 0 to 1
    B = B(:,1:2^(N));
    var_noise = 10^(-SNR_dB/10)*N;                  % Consider signal power to be N
    for loop2 = 1:ave_num
        H = randn(1,N);
        Y = sign(H*B+randn(1,length(B(1,:)))*sqrt(var_noise));
        B_hat = B.*repmat(Y,N,1);
        theta = sum(B_hat,2)/norm(B_hat);
        theta = theta';
        H_est = theta*0;
        for loop3 = 1:T
            H_est = H_est + theta/loop3;
            index = find(H_est*B_hat<=0);
            if(isempty(index))
                break;
            else
                i_select = ceil(rand(1)*length(index));
                theta = theta+(B_hat(:,i_select)).';
            end
            %     i_select = ceil(rand(1)*length(B(1,:)));
            %     if(H_est*B_hat(:,i_select)<0)
            %         theta = theta + (B_hat(:,i_select)).';
            %     end
        end
        H_est = H_est/T;
        if(norm(H_est)>0)
            H_est = H_est/norm(H_est)*sqrt(N);
            % loop2
        end
        SVM_alg_err = norm(H_est-H);
        SVM_error(loop1) = SVM_error(loop1)+SVM_alg_err^2;
        
        %% MMSE
        
        N_SIM = 1000;                          %origin 100000
        H_test = randn(N_SIM,N);
        y_test = (H_test*B);
        y_hat = y_test.*repmat(Y,N_SIM,1);
        prob_y = qfunc(-y_hat);
        prob_y1 = prod(prob_y,2);
        H_mmse = prob_y1'*H_test/sum(prob_y1);
        mmse = norm(H_mmse-H);
        mmse_error(loop1) = mmse_error(loop1)+mmse^2;
        
        if(mod(loop2,1000)==0)
            loop2
        end
    end
    SVM_error(loop1) = SVM_error(loop1)/ave_num/N;
    mmse_error(loop1) = mmse_error(loop1)/ave_num/N;
    N
end

%% Plot
figure
plot(N_test, mmse_error,'r');
hold on
plot(N_test, SVM_error,'b');
legend('mmse','SVM')
title(strcat('SNR=',num2str(SNR_dB),'dB'))
