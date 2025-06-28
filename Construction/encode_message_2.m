function c = encode_message_2(msg, H)
    % Encode a binary message 'm' using the parity-check matrix 'H'
    
    [num_rows, ~] = size(H);
    
    parity_len = num_rows;   
    p = false(1, parity_len);  % Initialize binary vector for parity bits
    
    % Compute the parity bits by solving the linear system from H
    for i = 1:parity_len
        % Find which bits of 'msg' correspond to this row in H
        parity_eq = H(i, parity_len+1:end); % The part of H involving message bits
        relevant_m = msg(logical(parity_eq)); % Extract relevant bits of 'm'
        
        % Parity equation: sum of relevant message bits plus current parity
        if i == 1 
            p(i) = mod(sum(relevant_m), 2); % Ensure binary format (mod 2)
        else
            p(i) = mod(sum(relevant_m) + p(i-1), 2); % last row also gives p(end)
        end
    end
    
    c = [p msg]; % Concatenate parity bits and message bits
end
