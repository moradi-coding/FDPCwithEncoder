close all; clc; clear; format short;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% You may change the parfor in MonteCarlo to for loop and comment this part
CoreNum = 4; % depends on your computer specifics.
poolobj = gcp('nocreate');
if isempty(poolobj)
    parpool(CoreNum);
else
    disp('matlab pool already started');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% output folder
dirname = 'output';
if ~exist(dirname, 'dir')
    mkdir(dirname)
end


% parameters
% H is constructed based on the paper's algorithm (look at the paper)
load('H_n256_k192_t16_rankH64_systematic_H_base1_FULLMatrix_NoCycle.mat')
H = H_base1.H;

[M,N] = size(H);  
K = N - rank(H);
Rc = K/N; % code rate

maxBP_iter = 50;  % maximum number of the decoding iterations

EbNo_vec = 0:.5:6;
maxRun = 1E12;
maxFE = [300*ones(1,length(EbNo_vec)-4) 300 100 100 100]; % maximum number of error frames
ClustNum = [1e3 1e3 1e3 1e3 2e3 1e4*ones(1,length(EbNo_vec)-5)]*10;

FER = zeros(1,length(EbNo_vec));
FE = zeros(1,length(EbNo_vec));  % #of frame errors
BER = zeros(1,length(EbNo_vec));
Nruns = zeros(1,length(EbNo_vec)); % #of actual runs

H_sparse = sparse(H);
H_sparse = logical(H_sparse);
decodercfg = ldpcDecoderConfig(logical(H_sparse),'norm-min-sum');  % 'bp', 'layered-bp', 'norm-min-sum', 'offset-min-sum'
disp(decodercfg.Algorithm)
fprintf('-------------------------------------\n');
for EbNo_count = 1:length(EbNo_vec)
    tic;
    EbNo_dB = EbNo_vec(EbNo_count);
    EbNo = 10^(EbNo_dB/10);
    sigma = 1/sqrt(2*Rc*EbNo); % for variance of the noise

    Nblkerrs = 0; % initialization for counting block errors
    Nbiterrs = 0; % initialization for counting bit errors
    fprintf('[%02d:%02d] Starting! SNR = %.1f\n',0,0, EbNo_dB);
    for i = 1:maxRun/ClustNum(EbNo_count)
        parfor j = 1:ClustNum(EbNo_count)

            msg = randi([0, 1], 1, K);
            c = encode_message_2(msg, H); % look at the paper for encoding
            % sum(mod(c*H',2))
            % % BPSK modulation  
            modulated = 1 - 2*c';
  
            % AWGN
            r = modulated + randn(N, 1)*sigma;
            % r = ones(N,1) + randn(N, 1)*sigma; 
            % Decoding
            Lch = 2*r./sigma.^2;  % calc LLRs
            
            [LLRxhat,~,~] = ldpcDecode(Lch,decodercfg,maxBP_iter,...
                DecisionType='soft', OutputFormat= 'whole', Termination = 'early');
%             numiter(i) = iter;    % save iteration count

            xhat = (-sign(LLRxhat))/2 + 0.5;
            msg_hat = xhat(M+1:end)';
            
            % NON systematic (obtaining FER based on only encoding the Zero vector
            % S = sum(msg_hat); 
            % if S
            %     Nblkerrs = Nblkerrs + 1;
            %     Nbiterrs = Nbiterrs + S;
            % end

            % systematic
            Nblkerrs = Nblkerrs + any(msg_hat~=msg);
            Nbiterrs = Nbiterrs + sum(msg_hat~=msg);
        end
        if Nblkerrs >= maxFE(EbNo_count)
            break;
        end
        t = toc;
        elapsed_m = t/60;
        elapsed_h = floor(elapsed_m/60);
        elapsed_m = floor(elapsed_m - elapsed_h*60);
        fprintf('[%02d:%02d] SNR = %.1f, Frame = %d, FE = %d\n',elapsed_h,elapsed_m,EbNo_dB,(i*ClustNum(EbNo_count)),Nblkerrs);
    end
    t = toc;
    elapsed_m = t/60;
    elapsed_h = floor(elapsed_m/60);
    elapsed_m = floor(elapsed_m - elapsed_h*60);
    temp = (i*ClustNum(EbNo_count));
    fprintf(2,'[%02d:%02d] SNR = %.1f, Frame = %d, FE = %d\n',elapsed_h,elapsed_m,EbNo_dB,temp,Nblkerrs);
    
    temp = (i*ClustNum(EbNo_count));
    FER(EbNo_count) = Nblkerrs/temp;
    FE(EbNo_count) = Nblkerrs;
    BER(EbNo_count) = Nbiterrs./ temp/K;
    Nruns(EbNo_count) = temp;
    fprintf('-------------------------------------\n');
end



% Creating sim_state, saving and plotting the results
res.Nruns = Nruns;
res.FE = FE;
res.FER = FER;
res.BER = BER;
res.maxBP_iter = maxBP_iter;

res.SNR = EbNo_vec;
res.K = K;
res.N = N;
res.maxRun = maxRun;
res.maxFE = maxFE;



semilogy(EbNo_vec(1:end),FER(1:end),'-*r','LineWidth',2)
hold on
semilogy(EbNo_vec(1:end),BER(1:end),'-ob','LineWidth',2)
grid on

xlabel('Eb/No (dB)','FontSize', 14)
ylabel('FER/BER','FontSize', 14);
legend('wEnc - norm-min-sum, (128, 82), Iter = 5','Interpreter', 'LaTeX','FontSize', 14);

filename = [dirname sprintf('/H_n%d_k%d_t32_systematic_base11_maxBP_iter%d_norm-min-sum_NoCycle.mat',N,K,maxBP_iter)];
save(filename,'res');
