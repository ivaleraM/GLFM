function [xd, pdf] = IBPsampler_computePDF(data, Zp, hidden, params, d)
    % Function to generate the PDF solutions corresponding to patterns in
    % Zp, and dimension d
    % Inputs:
    %   data.X: N*D data matrix (necessary to compute domain of x
    %   data.C: 1*D string with data types, D = number of dimensions
    %   Zp: P * K matrix of patterns (P is the number of patterns)
    %   hidden.B: latent feature matrix (D * K * maxR)   
    %   hidden.mu: 1*D shift parameter for internal transformation
    %   hidden.w:  1*D scale parameter for internal transformation
    %   hidden.Theta: D*maxR matrix of aux. variables for ordinal data
    %                 where R is the max number of categories
    %
    % Outputs:
    %   xd: 1*numS where numS is the number of points to be considered
    %  pdf: P*numS where P is the number of patterns to consider

    data.X(isnan(data.X(:,d)),d) = params.missing;

    % compute x-domain [mm MM] to compute pdf
    mm = min(data.X(data.X(:,d) ~= params.missing, d)); % min value
    MM = max(data.X(data.X(:,d) ~= params.missing, d)); % max value

    if ~isempty(params.t{d}) % if there is an external transformation
        data.C(d) = params.ext_dataType{d};
        mm = params.t_1{d}(mm);
        MM = params.t_1{d}(MM);
    end

    P = size(Zp,1);
    K = size(hidden.B,2);
    if (size(Zp,2) ~= K)
        error('Incongruent sizes between Zp and hidden.B');
    end
    if (data.C(d) == 'g') || (data.C(d) == 'p')
        if ~isfield(params,'numS')
            params.numS = 100;
        end
        xd = linspace(mm, MM, params.numS);
    elseif (data.C(d) == 'n')
        xd = mm:MM;
        params.numS = length(xd);
    else
        xd = unique(data.X(data.X(:,d) ~= params.missing, d));
        params.numS = length(xd); % number of labels for categories or ordinal data
    end
    pdf = zeros(P,params.numS);
    for p=1:P
        switch data.C(d)
            case 'g', pdf(p,:) = pdf_g(xd,Zp(p,:), squeeze(hidden.B(d,:,1))', hidden.mu(d), hidden.w(d), hidden.s2Y(d), params);
            case 'p', pdf(p,:) = pdf_p(xd,Zp(p,:), squeeze(hidden.B(d,:,1))', hidden.mu(d), hidden.w(d), hidden.s2Y(d), params);
            case 'n', pdf(p,:) = pdf_n(xd,Zp(p,:), squeeze(hidden.B(d,:,1))', hidden.mu(d), hidden.w(d), hidden.s2Y(d), params);
            case 'c', pdf(p,:) = pdf_c(Zp(p,:), squeeze(hidden.B(d,:,1:hidden.R(d))), hidden.s2Y(d));
            case 'o', pdf(p,:) = pdf_o(Zp(p,:), squeeze(hidden.B(d,:,1))', hidden.theta(d,1:(hidden.R(d)-1)), hidden.s2Y(d));
            otherwise
                error('Unknown data type');
        end
        if (sum(isnan(pdf)) > 0)
            error('Some values are nan!');
        end
    end
%     if ( data.C(d) == 'g' || data.C(d) == 'p' )
%         int = integrate(xd, pdf)
%     elseif (data.C(d) == 'n')
%         int = sum(pdf,2)
%     end

    if isfield(params,'t')
        if ~isempty(params.t{d}) % we have used a special transform beforehand
            xd = params.t{d}(xd); % if there was an additional external transformation, transform pdf
            pdf = pdf .* abs( params.dt_1{d}(xd) );
        end
    end
end
