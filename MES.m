
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
%     disp(B(0,0));
% disp("CDDCD");
%     disp(integrall(0,0));
% %     disp(integralBounded(max(low(0,0),low(0,0)),min(high(0,0),high(0,0)),0,0));
% %     disp(integralBounded(max(low(0,0),low(0,1)),min(high(0,0),high(0,1)),0,0));
% %     disp(integralBounded(max(low(0,1),low(0,0)),min(high(0,1),high(0,0)),0,0));
%     disp(integralBounded(max(low(0,1),low(0,1)),min(high(0,1),high(0,1)),0,0));
%     disp("cdcdcd");
%     disp(eiDiv(0,1/2-sqrt(3/5)*1/2)*eiDiv(0,1/2-sqrt(3/5)*1/2)*(5/9) + ...
%                 eiDiv(0,1/2)*eiDiv(0,1/2)*(8/9) + ...
%                 eiDiv(0,1/2+sqrt(3/5)*1/2)*eiDiv(0,1/2+sqrt(3/5)*1/2)*(5/9));
%     disp(3*ei(i1,0)*ei(i2,0));
% disp(eiDiv(0,0.5));
% disp((3/2)*(8/9)*eiDiv(0,0.5)*eiDiv(0,0.5));
    disp(full(M));
    disp(L);
    disp(wynik);
    % tworzymy wektor punktow (x)
    x = [a : h : b];

    % oraz wartosci (y)
    y = zeros(1 , n+1) ;
    for k = 1 : n
        y(k) = wynik(k);
    end

    plot(x, y) ;
    grid on

    function y=xi(i)
        y=2*i/n;
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
%         if max(xi(i1),xi(i2))<1
%             y=3*integralBounded(0,1,i1,i2) + 5*integralBounded(1,2,i1,i2);
%         elseif min(xi(i1),xi(i2))>1
%             y=3*integralBounded(0,1,i1,i2) + 5*integralBounded(1,2,i1,i2);
%         end
        res=0;
        res=res+integralBounded(max(low(i1,0),low(i2,0)),min(high(i1,0),high(i2,0)),i1,i2);
        res=res+integralBounded(max(low(i1,0),low(i2,1)),min(high(i1,0),high(i2,1)),i1,i2);
        res=res+integralBounded(max(low(i1,1),low(i2,0)),min(high(i1,1),high(i2,0)),i1,i2);
        res=res+integralBounded(max(low(i1,1),low(i2,1)),min(high(i1,1),high(i2,1)),i1,i2);
%         res=res+integralBounded(max(xi(i1-1),xi(i2)),min(xi(i1),xi(i2+1)),i1,i2);
%         res=res+integralBounded(max(xi(i1),xi(i2-1)),min(xi(i1+1),xi(i2)),i1,i2);
%         res=res+integralBounded(max(xi(i1),xi(i2)),min(xi(i1+1),xi(i2+1)),i1,i2);
%         y=integralBounded(0,1,i1,i2) + integralBounded(1,2,i1,i2);
        y=res;
    end

    function y=integralBounded(a,b,i1,i2)
%         disp(a);
%         disp(b);
        if a<b && a~=-inf
            shft = (a+b)/2;
            shft2 = (b-a)/2;
%             disp("begin");
%                             disp(eiDiv(i1,shft-sqrt(3/5)*shft2)*eiDiv(i2,shft-sqrt(3/5)*shft2)*(5/9));
%                 disp(eiDiv(i1,shft)*eiDiv(i2,shft)*(8/9));
%                 disp(eiDiv(i1,shft+sqrt(3/5)*shft2)*eiDiv(i2,shft+sqrt(3/5)*shft2)*(5/9));
            xy = eiDiv(i1,shft-sqrt(3/5)*shft2)*eiDiv(i2,shft-sqrt(3/5)*shft2)*(5/9) + ...
                eiDiv(i1,shft)*eiDiv(i2,shft)*(8/9) + ...
                eiDiv(i1,shft+sqrt(3/5)*shft2)*eiDiv(i2,shft+sqrt(3/5)*shft2)*(5/9);
%             disp("end");
            if a>=1
                y=shft2*xy*5;
%                 disp("d");
%                 disp(y);
            elseif b<=1
                y=shft2*xy*3;
%                 disp(y);
%                 disp(xy);
%                 disp(shft2);
            else
%                 disp("xx");
                shft1 = (a+1)/2;
                shft21 = (1-a)/2;
                xy1 = shft2*(...
                    eiDiv(i1,shft1-sqrt(3/5)*shft21)*eiDiv(i2,shft1-sqrt(3/5)*shft21)*(5/9) + ...
                    eiDiv(i1,shft1)*eiDiv(i2,shft1)*(8/9) + ...
                    eiDiv(i1,shft1+sqrt(3/5)*shft21)*eiDiv(i2,shft1+sqrt(3/5)*shft21)*(5/9));
                shft21 = (1+b/2);
                shft22 = (b-1/2);
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
