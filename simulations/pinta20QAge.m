%s2y=0.5;
minX=0;
vecADD=vec;
a=0.5;
%%a*x^2;
% f_1=@(x) sqrt(x/a);
% df_1=@(x) 1./(2*sqrt(a*x)); 

%%(a*x)^2;
f_1=@(x) (1/a)*sqrt(x); %a*x^2;
df_1=@(x) 1./(2*a*sqrt(x));

[N D]=size(X);
Zest=Zest';
Kest = Kest-1;

Zest'*Zest

ii=0;
maxR=max(R(1:D));

baseline=zeros(D,maxR);

prob=zeros(Kest+1,maxR,D);
prob2=[];

vec=[0 1 2 4 8 16 32 64];
for i=vec(1:Kest+1)%8 16 32 64
    ii=ii+1;
    z_dib=[1 de2bi(i,Kest)];
    for d=1:20
        Br=squeeze(B(d,:,:));
        prob(ii,1,d) =normcdf(z_dib*(Br(:,1)-Br(:,2))/sqrt(2),0,1);
        %prob(ii,1,d)=exp(prob(ii,1,d));
        prob(ii,2,d)=1-prob(ii,1,d);
        baseline(d,1)=sum(X(:,d)==1)/N;
        baseline(d,2)=sum(X(:,d)==2)/N;
    end
    
end

d=21; %edad
h_1=@(h,mu) normcdf((f_1(h+1)-mu)/s2y,0,1)-normcdf((f_1(h)-mu)/s2y,0,1);


Br=squeeze(B(d,:,:));
vec=[0 1 2 4];
i=vec(1);
z_dib=[1 de2bi(i,Kest)]                  
mu=z_dib*Br(:,1);

j=0;
for h=18:98
   j=j+1;
   y(1,j)=h_1(h-minX,mu); 
end

i=vec(2);
z_dib=[1 de2bi(i,Kest)]                  
mu=z_dib*Br(:,1);
j=0;
for h=18:98
   j=j+1;
   y(2,j)=h_1(h-minX,mu); 
end

i=vec(3);
z_dib=[1 de2bi(i,Kest)]                  
mu=z_dib*Br(:,1);
j=0;
for h=18:98
   j=j+1;
   y(3,j)=h_1(h-minX,mu); 
end

i=vec(4);
z_dib=[1 de2bi(i,Kest)]                  
mu=z_dib*Br(:,1);
j=0;
for h=18:98
   j=j+1;
   y(4,j)=h_1(h-minX,mu); 
end

j=0;
for h=18:98
   j=j+1;
   yb(1,j)=sum(X(:,d)==h-minX)/N; 
end



figure,
h=semilogy(1:5, squeeze(prob(1,1,(1:5))), 'rs-' ,1:5,squeeze(prob(2,1,(1:5))), 'cx-' ,1:5, squeeze(prob(3,1,(1:5))), 'mo-',1:5, squeeze(prob(4,1,(1:5))), 'b^-',1:5, baseline(1:5,1),'k+:' ,6:13,squeeze(prob(1,1,(6:13))), 'rs-' ,6:13, squeeze(prob(2,1,(6:13))), 'cx-'  ,6:13, squeeze(prob(3,1,(6:13))), 'mo-' ,6:13, squeeze(prob(4,1,(6:13))), 'b^-',6:13, baseline(6:13,1),'k+:',  14:20, squeeze(prob(1,1,(14:20))), 'rs-' , 14:20, squeeze(prob(2,1,(14:20))), 'cx-' , 14:20, squeeze(prob(3,1,(14:20))), 'mo-', 14:20, squeeze(prob(4,1,(14:20))), 'b^-',14:20, baseline(14:20,1) ,'k+:', 'linewidth',1.5, 'MarkerSize',10);
set(h([1 6 11]),'color',[0 .5 0]);
set(h([2 7 12]),'color',[1 0.1 0])
set(h([3 8 13]),'color','m')
set(h([4 9 14]),'color','b')
set(gca,'fontsize',20)
legend('[000]','[100]','[010]','[001]','Baseline','location','southwest')
set(h,'MarkerFaceColor','auto')
grid on

label={'1. Alcohol abuse', '2. Alcohol depend.', '3. Drug abuse', '4. Drug depend.', '5. Nicotine depend.', '6. MDD','7. Bipolar disorder','8. Dysthimia','9. Panic disorder','10. SAD','11. Specific phobia','12. GAD','13. PG', '14. Avoidant PD','15. Dependent PD','16. OCPD', '17. Paranoid PD', '18. Schizoid PD', '19. Histrionic PD','20. Antisocial PD'}';
hText = xticklabel_rotate(1:20,65,label,'Fontsize',20);
figurapdf(20, 8)
% print -dpdf 20Q+_4.pdf


figure,h1=stem(18:98,y(1,:),'rs')
set(h1([1]),'color',[0 .5 0]);
hold on,  stem(18:98,y(2,:), 'rx')
%set(h1([1]),'color',[1 0.1 0])
hold on, stem(18:98,y(3,:), 'mo')
%set(h1([1]),'color','m')
hold on, stem(18:98,y(4,:), 'b^')
hold on, stem(18:98,yb(1,:), 'k+')
set(gca,'fontsize',20)
legend('[000]','[100]','[010]','[001]','Baseline','location','northeast')

figurapdf(15, 5)
% print -dpdf 21Q_4.pdf

% for i=1:size(y,2)
% cdfy(:,i)=sum(y(:,1:i),2);
% cdfyb(1,i)=sum(yb(:,1:i),2);
% end
% 
% figure,stem(18:98,cdfy')
% hold on,stem(18:98,cdfyb','k')

