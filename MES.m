
function MES(n)
    a = 0;
    b = 2;
    h = (b-a)/n;

    M = sparse(n, n);

    M(1, 1) = B(0,0) ;
    M(1, 2) = B(0,1) ;
    for k = 2 : (n-1)
        M(k, k-1) = B(k-1,k-2) ;
        M(k, k) = B(k-1,k-1) ;
        M(k, k+1) = B(k-1,k) ;
    end
    M (n, n-1) = B(n-1,n-2) ;
    M (n, n) = B(n-1,n-1) ;


    L = zeros(n, 1) ;
    for k = 1 : (n)
        L(k , 1) = -30*ei(k-1,0);
    end

    wynik = M\L;
    disp(full(M));
    disp(L);
    disp(wynik);

    x = [a : h : b];

    y = zeros(1 , n+1) ;
    for k = 1 : n
        y(k) = wynik(k);
    end

    plot(x, y) ;
    grid on

    function y=xi(i)
        y=(2*i)/n;
    end
    
    function y=ei(i,x)
        if x>=xi(i-1) && x<=xi(i) && i-1>=0
            y=(x-xi(i-1))/(xi(i)-xi(i-1));
        elseif x>=xi(i) && x<=xi(i+1)
            y=(xi(i+1)-x)/(xi(i+1)-xi(i));
        else
            y=0;
        end
    end

    function y=eiDiv(i,x)
        if x>=xi(i-1) && x<=xi(i) && i-1>=0
            y=1/(xi(i)-xi(i-1));
        elseif x>=xi(i) && x<=xi(i+1)
            y=-1/(xi(i+1)-xi(i));
        else
            y=0;
        end
    end

    function y=B(i1, i2)
        y=integrall(i1, i2) - 3*ei(i1,0)*ei(i2,0);
    end

    function y=low(i,mod)
        if i-1>=0 && mod==0
            y=xi(i-1);
        elseif i<=n
            y=xi(i);
        else
            y=-inf;
        end
    end

    function y=high(i,mod)
        if mod==0
            y=xi(i);
        elseif i+1<=n
            y=xi(i+1);
        else
            y=-inf;
        end
    end

    function y=integrall(i1,i2)
        res=    integralBounded(max(low(i1,0),low(i2,0)),min(high(i1,0),high(i2,0)),i1,i2);
        res=res+integralBounded(max(low(i1,0),low(i2,1)),min(high(i1,0),high(i2,1)),i1,i2);
        res=res+integralBounded(max(low(i1,1),low(i2,0)),min(high(i1,1),high(i2,0)),i1,i2);
        res=res+integralBounded(max(low(i1,1),low(i2,1)),min(high(i1,1),high(i2,1)),i1,i2);
        y=res;
    end

    function y=integralBounded(a,b,i1,i2)
        if a<b && a~=-inf
            shft = (a+b)/2;
            shft2 = (b-a)/2;
            xy = eiDiv(i1,shft-sqrt(3/5)*shft2)*eiDiv(i2,shft-sqrt(3/5)*shft2)*(5/9) + ...
                eiDiv(i1,shft)*eiDiv(i2,shft)*(8/9) + ...
                eiDiv(i1,shft+sqrt(3/5)*shft2)*eiDiv(i2,shft+sqrt(3/5)*shft2)*(5/9);
            if a>=1
                y=shft2*xy*5;
            elseif b<=1
                y=shft2*xy*3;
            else
                shft1 = (a+1)/2;
                shft21 = (1-a)/2;
                xy1 = shft21*(...
                    eiDiv(i1,shft1-sqrt(3/5)*shft21)*eiDiv(i2,shft1-sqrt(3/5)*shft21)*(5/9) + ...
                    eiDiv(i1,shft1)*eiDiv(i2,shft1)*(8/9) + ...
                    eiDiv(i1,shft1+sqrt(3/5)*shft21)*eiDiv(i2,shft1+sqrt(3/5)*shft21)*(5/9));
                shft21 = (1+b)/2;
                shft22 = (b-1)/2;
                xy2 = shft22*(...
                    eiDiv(i1,shft21-sqrt(3/5)*shft22)*eiDiv(i2,shft21-sqrt(3/5)*shft22)*(5/9) + ...
                    eiDiv(i1,shft21)*eiDiv(i2,shft21)*(8/9) + ...
                    eiDiv(i1,shft21+sqrt(3/5)*shft22)*eiDiv(i2,shft21+sqrt(3/5)*shft22)*(5/9));
                y=3*xy1+5*xy2;
            end
        else
            y=0;
        end
    end
end
