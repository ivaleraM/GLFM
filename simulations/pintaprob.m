%load (['results/20QmasKini3add13y14_it50'])
[N D]=size(X);
Zest=Zest';
Kest = Kest-1;

Zest'*Zest

ii=0;

prob=zeros(2^Kest+1,maxR,D);
for i=0:2^Kest-1
    ii=ii+1;
    z_dib=[1 de2bi(i,Kest)];
    for d=1:D
        Br=squeeze(B(d,:,:));
        for r=1:R(d)-1
           for r2=1:R(d)
                if r2~=r
                    prob(ii,r,d) =prob(ii,r,d)+ log(normcdf(z_dib*(Br(:,r)-Br(:,r2)),0,1));
                end
           end
        end
        prob(ii,1:R(d)-1,d)=exp(prob(ii,1:R(d)-1,d));
        prob(ii,R(d),d)=1-sum(prob(ii,1:R(d)-1,d),2);
    end
end
baseline=zeros(D,maxR);
for d=1:D
    for r=1:R(d)
        baseline(d,r)=sum(X(:,d)==r)/N;
    end
end

% figure,
% for k=1:Kest+1
%     subplot(1,Kest+1,k)
%     imshow(1-squeeze(prob(k,:,:)).',[0 1])
% end   
% for i=1:2^Kest+1
%     figure,plot(1:D,squeeze(sum(prob([i],[1],:),2))','x-',1:D,baseline(:,1),'o--')
% end

