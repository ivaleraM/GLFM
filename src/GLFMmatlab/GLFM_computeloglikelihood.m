function loglik = GLFM_computeloglikelihood(Xtrue_ND, C, hidden, params)
    %% This function calculates the log-lik of
    %     the held-out data in the GLFM. The function uses the trained
    %     parameters,
    %     to calculate the PDF at the points of missing values, using the true values.
    %     Inputs:
    %         - Xtrue_ND: N*D array of true data. The function computes the log-lik
    %                    values for all no nan values in Xtrue
    %         - C: D      string with data types
    %         - hidden:   inferred latent variables from the GLFM
    %         - params:   hyperparameters from the GLFM
    %     Output:
    %         - loglik: N*D matrix with the evaluation of the log-likelihood for every value in Xtrue_ND.

    D = size(Xtrue_ND,2);
    N = size(Xtrue_ND,1);
    
    loglik = zeros(N,D);
    
    Zp = hidden.Z;
    
    K = size(hidden.B,2);
    if (size(Zp,2) ~= K)
        error('Incongruent sizes between Zp and hidden.B: number of latent variables should not be different');
    end
    
    for d=1:D
        xd = Xtrue_ND(:,d);
        if isfield(params,'t') % if external transformations have been defined
            if ~isempty(params.t{d}) % there is an external transform for data type d
                C(d) = params.ext_dataType{d}; % set new type of data
            end
        end
        switch C(d)
            case 'g', 
                pdf = pdf_g(Xtrue_ND(:,d),Zp, squeeze(hidden.B(d,:,1))', hidden.mu(d), hidden.w(d), hidden.s2Y(d), params);
                if isfield(params,'t')
                    if ~isempty(params.t{d}) % we have used a special transform beforehand
                        xd = params.t{d}(xd); % if there was an additional external transformation, transform pdf
                        pdf = pdf .* abs( params.dt_1{d}(xd) );
                    end
                end
                loglik(:,d) = log(pdf);
            case 'p'
                pdf = pdf_p(xd,Zp, squeeze(hidden.B(d,:,1))', hidden.mu(d), hidden.w(d), hidden.s2Y(d), params);
                if isfield(params,'t')
                    if ~isempty(params.t{d}) % we have used a special transform beforehand
                        xd = params.t{d}(xd); % if there was an additional external transformation, transform pdf
                        pdf = pdf .* abs( params.dt_1{d}(xd) );
                    end
                end
                loglik(:,d) = log(pdf);
            case 'n'
                pdf = pdf_n(xd,Zp, squeeze(hidden.B(d,:,1))', hidden.mu(d), hidden.w(d), hidden.s2Y(d), params);
                loglik(:,d) = log(pdf);
            case 'c' 
                pdf = pdf_c(Zp, squeeze(hidden.B(d,:,1:hidden.R(d))), hidden.s2Y(d));
                idx= R.*[0:N-1]'+xd;
                idx(isnan(idx))=1;
                loglik(:,d) = log(pdf(idx)'); 
                loglik(isnan(Xtrue_ND(:,d)),d) = nan;
            case 'o' 
                pdf = pdf_o(Zp, squeeze(hidden.B(d,:,1))', hidden.theta(d,1:(hidden.R(d)-1)), hidden.s2Y(d));
                idx= R.*[0:N-1]'+ xd;
                idx(isnan(idx))=1;
                loglik(:,d) = log(pdf(idx)'); 
                loglik(isnan(xd)) = nan;
            otherwise
                error('Unknown data type');
        end
        
    end