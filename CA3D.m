% 3-D Coherence Analysis Algorithm
% Generate the mode-by-mode correlation matrix of a 3-D radiation pattern (theta,phi,frequency).
% Modes can be interpreted as radiation patterns at specific frequencies in specific masks.
% input variable patterns is a 3-D matrix

function R = CA3D(patterns)

    % Initialize: Define Frequency Sweep and Spatial Parameters   
    H = patterns;
    [y_dim, x_dim, Nk] = size(H);

    % Generate the correlation matrix (R) for current mask
    R = zeros(Nk,Nk); % Create container for pattern correlation (maskNo)
    % Iterate frequencies/modes in current mask
    for f_i=1:Nk
        H_i = H(:,:,f_i); % Get 2D-pattern (theta,phi) for freq_i, H_i: (THETA,PHI;f_i)
        H_i_norm = H_i./abs(H_i);
        
        for f_j=1:Nk
            H_j = H(:,:,f_j); % Get 2D-pattern (theta,phi) for freq_j, H_j: (THETA,PHI;f_j)
            H_j_norm = H_j./abs(H_j);
            
            % Perform dot-product operation of the two patterns
            R_ij= sum(dot(H_i_norm,H_j_norm)/y_dim)/x_dim; % divide by theta_dim and phi_dim to rescale the resulting large value and to prevent storage issues. No effect otherwise.

            % Fill in the correlation matrix
            R(f_i,f_j) = R_ij;
            if f_j ~= f_i
                R(f_j,f_i) = R_ij;
            end
        end
    end
end
