% 3-D SVD ~ Mastorakis (1996)
% Calculate the singular values of a 3-D matrix based on Mastorakis (1996) algorithm
% Input variable is a 3-D matrix
% When n ≥ m, the matrix Σ has at most m non-zero elements on the diagonal, and may be written as Σ=[ˆΣ  0]'.
% Therefore, it is possible to exactly represent X using the economy SVD [1].

% [1] Brunton SL, Kutz JN. «Data-Driven Science and Engineering: Machine Learning, Dynamical Systems, and Control». Cambridge University Press; 2019.

function [sigma,u,v,w] = SVD3D(matrix3D)

    % Get matrix dimensions
    [m1, m2, m3] = size(matrix3D);
    A = matrix3D;
    B = zeros(m1,m2*m3);
    
    % a) From 3-D matrix A [m1 x m2 x m3], obtain equivalent 2-D matrix B [m1 x m2m3]
    for i2=1:m2
        B(:,1+(i2-1)*m3:m3+(i2-1)*m3) = A(:,i2,:);
    end
    
    % b) Obtain SVD of [m1 x m2m3] matrix B
    % U_B [m1 x m1] unitary matrix
    % V_B [m2m3 x m2m3] unitary matrix
    % S_B [m1 x m2m3] rectangular diagonal matrix with non-negative real numbers sig_B on its diagonal
    % sig_B: sig_1B,sig_2B,...,sig_rnB are the non-negative singular values of S_B
    % r_B: rank of S_B
    % u_sB is the sth column of U_B for s=1,2,...,r_B
    % v_sB is the sth column of V_B for s=1,2,...,r_B
    [U_B,S_B,V_B] = svd(B,'econ'); % use economy svd
    sig_B = diag(S_B); % sig_B: sig_1B,sig_2B,...,sig_r_BB
    r_B = length(sig_B); % Rank of S_B
    
    % c) inverse transform  sth column of V_B (V_sB [m2m3 x 1]) to get B_sB [m2 x m3]
    % B_sB [m2 x m3]
    % v_s [m2m3 x 1] sth column of V_B
    % B_B: [r_B x m2 x m3] B_1B,B_2B,...,B_sB
    B_B = zeros(r_B,m2,m3); % [r_B x m2 x m3]
    for s=1:r_B
        B_sB = zeros(m2,m3);
        v_sB = V_B(:,s); % v_sB is the sth columb of V_B
        for i2=1:m2
            B_sB(i2,1:m3) = v_sB(1+(i2-1)*m3:i2*m3); % [m2 x m3]
        end
        B_B(s,:,:) = B_sB; % B_B: [r_B x m2 x m3] B_1B,B_2B,...,B_sB
    end
    
    % d) For s=1,..,r_B obtain SVD of B_sB
    % S_sB [m2 x m3] is rectangular diagonal matrix with non-negative real numbers sig_sB on its diagonal
    % sig_sB: sig_s1,sig_s2,...,sig_stao_s are the non-negative singular values of S_sB
    % tao_s rank of S_B
    % v_sj is the jth column of V_sB
    % w_sj is the jth column of W_sB
    r = 0; % r = tao_1+tao_2+...+tao_r_B
    tao = zeros(1,r_B+1); % tao = [tao_0=0, tao_1, tao_2, ..., tao_r_B], tao_0=0
    % A = sig_1*u_1ov_1ow_1 + sig_2*u_2ov_2ow_2 + ... + sig_r*u_rov_row_r
    
    sig = cell(r_B,1); % e) sig_ttao_s
    for s = 1:r_B
        % Obtain SVD of B_sB
        B_sB = squeeze(B_B(s,:,:));
        [V_sB,S_sB,W_sB] = svd(B_sB);
        sig_sB = diag(S_sB); % sig_sB: sig_s1,sig_s2,...,sig_stao_s
        tao_s = length(sig_sB); % tao_s rank of S_B
        sig{s} = sig_sB; % e) sig_ttao_s
        tao(s+1) = tao_s; % tao_0=0
        r = r + tao_s; % r = tao_1+tao_2+...+tao_r_B
    end
    
    sigma = zeros(r,1);
    for t=1:r
        for s=1:r_B
            sig_tB = sig_B(s); % sig_sB or sig_B(s) is one element from sig_B : [sig_1B,sig_2B,...,sig_r_BB]
            sig_s = sig{s};  % sig{s} is sig_sB from sig{s} = [sig_1tao_1,sig_2tao_2,...,sig_stao_s,...,sig_rBtaor_B]
            sigma(1+sum(tao(1:s)):sum(tao(1:s+1))) = sig_tB*sig_s; % sigma_t = sigma_tB * sigma_tT_s
        end
    end
    sigma = sort(sigma/max(sigma), 'descend');
    
    % Initialize vectors. 2nd dimension is the maximum index that spans the array
    u = zeros(m1,sum(tao));
    v = zeros(m2,sum(tao));
    w = zeros(m3,sum(tao)); 
    for s=1:r_B
        sum_tao_k = sum(tao(1:s)); % sum_k_to_s-1(tao_k), tao_0=0
        u_sB = U_B(:,s); % [m1 x 1] u_sB is the sth columb of U_B [m1 x m1] : [u_1B, u_2B, ..., u_sB, ..., u_r_BB | ... u_m2B]
        tao_s = tao(s+1); % tao = [tao_0=0, tao_1, tao_2, ..., tao_r_B], tao_0=0+
        for j=1:tao_s
            v_sj = V_sB(:,j); % [m2 x 1] v_sj is the jth column of V_sB [m2 x m2] : [v_s1, v_s2, ..., v_sj, ..., v_stao_s | ... v_m2]
            w_sj = W_sB(:,j); % [m3 x 1] w_sj is the jth column of W_sB [m3 x m3] : [w_s1, w_s2, ..., w_sj, ..., w_stao_s | ... w_m2]
    
            % e) SVD of A is achieved through the following. (Substituting B_sB for sth column of V_B)
            % sig: sig_1,sig_2,...,sig_r are the non-negative singular values of 3-D matrix A
    	    jsum_idx = j+sum_tao_k;
            u(:,jsum_idx) = u_sB; % [m1 x sum(tao)]
            v(:,jsum_idx) = v_sj; % [m2 x sum(tao)]
            w(:,jsum_idx) = w_sj; % [m3 x sum(tao)]
       end
    end
end