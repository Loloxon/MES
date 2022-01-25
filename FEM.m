% rozwiazuje równanie:  
% -u'' + 3u = f(x)    - f(x) zdefiniwane na dole
% u(1) = 0, u(3) = 0 - jednorodne warunki Dirichleta

function FEM(a,b,c,n)
    % o obliczamy delte miedzy punktami
    h = (b-a)/n;

    % obliczami wartosci elementów macierzy przekątniowej
    P1 = (c*2*h)/3+2/h ;
    P2 = c*h/6-1/h ;

    % tworzymy rzadką macierz o rozmiarze n-1 na n-1
    macierz_M = sparse(n-1, n-1);

    % tworzymy wektor B o rozmiarze n-1
    wektor_B = zeros(n-1, 1) ;

    % wypełniamy macierz
    macierz_M(1, 1) = P1 ;
    macierz_M(1, 2) = P2 ;
    for k = 2 : (n-2)
        macierz_M(k, k-1) = P2 ;
        macierz_M(k, k) = P1 ;
        macierz_M(k, k+1) = P2 ;
    end
    macierz_M (n-1, n-2) = P2 ;
    macierz_M (n-1, n-1) = P1 ;

    % Wypełniamy wektor prawej strony (B)
    % TODO: te współcznniki trzeba policzyć z całek
    for k = 1 : (n-1)
        wektor_B(k , 1) = f(a+(k-1)*h)*h/6 + f(a+k*h)*2*h/3 + f(a+(k+1)*h)*h/6;
    end

    % eliminacja gaussa dla macierzy n-1 na n-1
    % i wektora prawej strony
    % matlab automatycznie używa optymalnego solvera dla macierzy rzadkich
    % TODO trzeba to rozwiązać przy pomocy własnego algorytmu raczej
    wynik = macierz_M\wektor_B ;
    
    % funkcja przesuniecie (shift function)


    % tworzymy wektor punktow (x)
    x = [a : h : b];

    % oraz wartosci (y)
    y = zeros(1 , n+1) ;
    for k = 2 : n
        y(k) = wynik(k-1);
    end
    

    % plotowanie wyniku
    x_u = a : 0.01 : b;
    u = roz_dokladne(x_u);
    plot(x_u , u) ;

    % zachowanie wykresu
    hold on
    plot(x, y) ;
    % koniec zachowania, kolejen rozwiazanie nadpisze wykres
    hold off
end

function y=f(x)
    y = (-x^2 + 3 + 3*(x^2 - 4*x + 3))*exp(x);
end

function y=roz_dokladne(x)
    y = (x.^2 - 4*x + 3).*exp(x);
end

% function y=f(x)
%     y = (-2 + 3 * (x .^2 + x - 2));
% end
% 
% function y=roz_dokladne(x)
%     y = (x .^2 + x - 2);
% end
% 
% function y=shift(y1, y2, a, b, x)
%     y = (y2 - y1)/(b - 2) * (x - a) + y1;
% end