function X_pred = complete_matrix(params,postDist)
    % Returns N x D matrix with filled out positions
    % params.beta: [Dall×K double]
    % postDist.mean: [K×N double]
    [K,N] = size(postDist.mean);
    y_mu_post = (params.beta  * postDist.mean)';
    y_pred = mvnrnd(y_mu_post,postDist.noiseCovMat)'; % D x N
    %y_pred = postDist.y;
    if isfield(params,'betaCont')
        Dcon = size(params.betaCont,1);
    else
        Dcon = 0;
    end
    if isfield(params,'psi')
        Dcat = size(params.psi,1);
        Ddis = length(params.A);
    else
        Dcat = 0;
        Ddis = 0;
    end
    X_pred.continuous = nan * ones(Dcon,N);
    X_pred.categorical = nan * ones(Dcat,N);
    X_pred.discrete = nan * ones(Ddis,N);
    % fill out X_pred
    X_pred.continuous = y_pred(1:Dcon,:);
    offset = 0;
    for i=1:length(params.nClass)
        nc = params.nClass(i)-1;
        Xcat_tmp = params.A{i} * y_pred((Dcon+offset+1):(Dcon+offset+nc),:) ...
            - postDist.b_psi{i}; % L x N
        X_pred.categorical(offset+1:(offset+nc),:) = Xcat_tmp;
        % fill out X_pred.discrete
        M_candidates = [eye(nc), zeros(nc,1)]; % M+1 x M+1
        for n=1:N
            distance_vec = compute_dist(Xcat_tmp(:,n),M_candidates); % 1 x (M+1)
            [~,idx] = min(distance_vec);
            X_pred.discrete(i,n) = idx;
        end
        offset = offset + nc;
    end
    
    function distance_vec = compute_dist(X_vec,Y)
        M = size(Y,2);
        distance_vec = sqrt(sum((repmat(X_vec,1,M) - Y).^2,1));