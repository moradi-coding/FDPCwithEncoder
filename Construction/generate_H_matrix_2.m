function H = generate_H_matrix_2(t)
    H = zeros(2*t, t^2);        
    col = 1; % Index to start filling columns in the matrix        
    gap = 0; % Start the gap between ones, initially no gap (0 zeros between ones)
    
    while col <= t^2
        % Fill columns according to the current gap between ones
        for row = 1:(2*t-gap-1)
            % If column exceeds t^2, break the loop
            if col > t^2
                break;
            end            
            % Place the ones in the column based on the current row and gap
            H(row, col) = 1;
            H(row+gap+1, col) = 1;                        
            col = col + 1; % Move to the next column
        end                
        gap = gap + 2; % Increase the gap by 2 for the next block of columns
    end
end

