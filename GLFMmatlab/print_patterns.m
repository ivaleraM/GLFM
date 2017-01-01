function print_patterns(fid,X,patterns,C)

str = sprintf('meanTPFD=%.2f, medianTFPD=%.2f', ...
    mean_not_isnan(X(:,end),0), median_not_isnan(X(:,end),0) );
fprintf(fid,['\t\t\t\t\t',str, '\n']);
for r=1:size(patterns,1)
    str = sprintf('%2d. Pattern: %s numPat=%2d, meanTPFD=%.2f, medianTFPD=%.2f\n', r, num2str(patterns(r,:)), sum(C == r), ...
        mean_not_isnan(X(C==r,end),0), median_not_isnan(X(C==r,end),0) );
    fprintf(fid,str);
end