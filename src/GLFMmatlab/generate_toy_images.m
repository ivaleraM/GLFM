function [data,gT] = generate_toy_images(N,s2x)
    % Inputs:
    %   s2x: observation noise
    %     N: number of samples
    %
    % Outputs:
    %  data: data structure with X: N*D obs matrix and C: 1*D datatype
    %  string vector
    %
    Btrue = 2* [0,1.0,0,0,0,0,  1,1,1,0,0,0, 0,1,0,0,0,0, ...
        0,0,0,0,0,0, 0,0,0,0,0,0, 0,0,0,0,0,0; ...
        0,0.0,0,1,1,1,  0,0,0,1,0,1, 0,0,0,1,1,1, ...
        0,0,0,0,0,0, 0,0,0,0,0,0, 0,0,0,0,0,0; ...
        0,0,0,0,0,0, 0,0,0,0,0,0, 0,0,0,0,0,0, ...
        1,0,0,0,0,0, 1,1,0,0,0,0, 1,1,1,0,0,0; ...
        0,0,0,0,0,0, 0,0,0,0,0,0, 0,0,0,0,0,0, ...
        0,0,0,1,1,1, 0,0,0,0,1,0, 0,0,0,0,1,0];
    
    K = size(Btrue,1); % number of latent features
    D = size(Btrue,2); % number of dimensions
    
    Ztrue = rand(N,K) < 0.5;
    
    data.X = sqrt(s2x) * randn(N,D) + ( Ztrue * Btrue );
    data.C = repmat('g',1,D);
    
    gT.B = Btrue;
    gT.Z = Ztrue;
    
end