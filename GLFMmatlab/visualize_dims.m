%f_1=@(x,w) log(exp(w*x)-1);

figure();
for d=1:D
    %if (data.C(d) == 'p')
        mask = Xmiss(:,d)~=missing;
        
        hist(Xmiss(mask,d), 100); hold off;

                
        M = max(Xmiss(mask,d));
        m = min(Xmiss(mask,d));
        fprintf('Feature %2d: %s, max=%f, min=%.2f\n',d,data.C(d),M,m);
%         w = 2/max(Xmiss(:,d));
%         max( f_1(Xmiss(mask,d),w) )
%         min( f_1(Xmiss(mask,d),w) )
%         if sum(isinf(f_1(Xmiss(mask,d),w))) > 0
%             fprintf('d=inf for dimension %d',d);
%         elseif sum(isnan(f_1(Xmiss(mask,d),w))) > 0
%             fprintf('d=nan for dimension %d',d); 
%         end
        pause;
    %elseif (data.C(d) == 'n')
    %    mask = Xmiss(:,d)~=missing;
    %    
    %    hist(Xmiss(mask,d), 100); hold off;
    %end
end
        