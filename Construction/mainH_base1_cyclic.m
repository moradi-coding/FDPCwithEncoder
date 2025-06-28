clear; clc; close all
% in this code, we remove 2t-3, and the the last column 1 so that we can replace 2t-2 zero
% columns after the first 2t-1 columns. So our columns will be 2t-1, 2t-2,
% 2t-5, 2t-7, 2t-9, ... 9, 7, 5, 3.
% In this case, t should be at least 3.

% output folder
dirname = 'output';
if ~exist(dirname, 'dir')
    mkdir(dirname)
end

rng(101)

t = 16; % row_sums = t; t must be >= 3s
num_per = 1;
n_target = 256; % n should be bigger than n_target


m_base = 2*t;
n = t^2;
m_size = m_base + num_per*m_base;

H = generate_H_matrix_2(t); % base matrix
C = H(:,m_size + 1:end);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% permutations
for i = 1:num_per
    col_permutation = randperm(size(C, 2));
    C_permuted = C(:, col_permutation); % Apply the column permutation to the matrix
    H_permuted = zeros(m_base, n);
    H_permuted(:,m_size + 1:end) = C_permuted;
    H = [H; H_permuted];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% double diagonal for the encoder
H(1:m_size,1:m_size) = 0;
for i = 1: m_size-1
    H(i,i) = 1;
    H(i+1, i) = 1;
end
H(m_size,m_size) = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% removing extra columns to get n_target columns
if n_target ~= n
    H(:,m_size+1:m_size+1 + n-n_target  - 1) = [];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



m = rank(H);
K = n_target - m;
disp(K/n_target)
msg = randi([0, 1], 1, K);
c = encode_message_2(msg, H);
sum(mod(c*H',2))

H_base1.H = H;
H_base1.t = t;
H_base1.m = m;


filename = [dirname sprintf('/H_n%d_k%d_t%d_rankH%d_systematic_H_base1_FULLMatrix_NoCycle.mat',n_target,K,t,m)];
save(filename,'H_base1');







