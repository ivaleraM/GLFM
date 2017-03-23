function print_patterns(fid,Y,patterns,C,label)

str = sprintf('%35s mean=%.2f, median=%.2f\n', strtrim(label), ...
    mean_not_isnan(Y,0), median_not_isnan(Y,0) );
fprintf(fid,str);
for r=1:size(patterns,1)
    str = sprintf('%2d. Pattern: %s numObs=%2d, mean=%.2f, median=%.2f\n', r, num2str(patterns(r,:)), sum(C(~isnan(Y)) == r), ...
        mean_not_isnan(Y(C(~isnan(Y))==r),0), median_not_isnan(Y(C(~isnan(Y))==r),0) );
    %str = sprintf('%2d. Pattern: %s numPat=%2d, meanTPFD=%.2f, medianTFPD=%.2f\n', r, num2str(patterns(r,:)), sum(C == r), ...
    %    mean_not_isnan(X(C==r,end),0), median_not_isnan(X(C==r,end),0) );
    fprintf(fid,str);
end