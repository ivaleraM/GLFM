function x = ford(y, theta)
    % Mapping function for ordinal data
    % Inputs:
    %       y: [1*R] Pseudo-observations
    %   theta: [1*(R-1)] Thresholds that divide the real line into R regions
    for r = 1:length(theta)
        val = theta();
        if (y < val)
            x = r;
            break;
        end
        if r == size(theta,3);
            x = r+1;
        end
    end 
end