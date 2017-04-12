function int = integrate(xd, pdf)
    int = zeros(size(pdf,1),1);
    
    for i=1:(size(pdf,2)-1)
        int = int + pdf(:,i)*(xd(i+1)-xd(i));
    end