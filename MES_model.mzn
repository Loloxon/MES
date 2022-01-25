int: n;
int: a = 0;
int: b = 2;
float: h = (b-a)/n;

array[1..n,1..n] of var float: M;
array[1..n] of var float: L;
array[1..n] of var float: result;

function var float: xi(var int: i) = (2*i)/n;

function var float: ei(var int: i, var float: x) = 
    if x>=xi(i-1) /\ x<=xi(i) /\ i-1>=0
        then (x-xi(i-1))/(xi(i)-xi(i-1))
        else	
           if x>=xi(i) /\ x<=xi(i+1)
            		then (xi(i+1)-x)/(xi(i+1)-xi(i))
		            else 0
		         endif
	   endif;

function var float: eiDir(var int: i, var float: x) = 
    if x>=xi(i-1) /\ x<=xi(i) /\ i-1>=0
        then 1/(xi(i)-xi(i-1))
        else	
           if x>=xi(i) /\ x<=xi(i+1)
            		then (-1)/(xi(i+1)-xi(i))
		            else 0
		         endif
	   endif;

function var float: B(var int: i1, var int: i2) = integral(i1,i2) - 3*ei(i1,0)*ei(i2,0);

function var float: low(var int: i, var int: mode) = 
    if i-1>=0 /\ mode=0
        then xi(i-1)
        else
            if i<=n
                then xi(i)
                else -infinity
            endif
    endif;

function var float: high(var int: i, var int: mode) = 
    if mode=0
        then xi(i)
        else
            if i+1<=n
                then xi(i+1)
                else -infinity
            endif
    endif;

function var float: integral(var int: i1, var int: i2) = 
    integralB(max(low(i1,0),low(i2,0)), min(high(i1,0), high(i2,0)), i1, i2)+
    integralB(max(low(i1,0),low(i2,1)), min(high(i1,0), high(i2,1)), i1, i2)+
    integralB(max(low(i1,1),low(i2,0)), min(high(i1,1), high(i2,0)), i1, i2)+
    integralB(max(low(i1,1),low(i2,1)), min(high(i1,1), high(i2,1)), i1, i2);

function var float: integralB(var float: a, var float: b, var int: i1, var int: i2) = 
    if a<b /\ a!=-infinity then
        if a<1 /\ b>1 then
        let{
        var float: mid1 = (a+1)/2;
        var float: move1 = (1-a)/2;
        var float: res1 = eiDir(i1,mid1-sqrt(3/5)*move1)*eiDir(i2,mid1-sqrt(3/5)*move1)*(5/9) + 
                          eiDir(i1,mid1)*eiDir(i2,mid1)*(8/9) + 
                          eiDir(i1,mid1+sqrt(3/5)*move1)*eiDir(i2,mid1+sqrt(3/5)*move1)*(5/9);
        var float: mid2 = (1+b)/2;
        var float: move2 = (b-1)/2;
        var float: res2 = eiDir(i1,mid2-sqrt(3/5)*move2)*eiDir(i2,mid2-sqrt(3/5)*move2)*(5/9) + 
                          eiDir(i1,mid2)*eiDir(i2,mid2)*(8/9) + 
                          eiDir(i1,mid2+sqrt(3/5)*move2)*eiDir(i2,mid2+sqrt(3/5)*move2)*(5/9);} in
                          move1*res1*3+move2*res2*5
                      else
        let{
        var float: mid = (a+b)/2;
        var float: move = (b-a)/2;
        var float: res = eiDir(i1,mid-sqrt(3/5)*move)*eiDir(i2,mid-sqrt(3/5)*move)*(5/9) + 
                         eiDir(i1,mid)*eiDir(i2,mid)*(8/9) + 
                         eiDir(i1,mid+sqrt(3/5)*move)*eiDir(i2,mid+sqrt(3/5)*move)*(5/9);} in 
                         if a>=1 then move*res*5 else move*res*3 endif
        endif        
    else 0 endif;

constraint M[1,1] = B(0,0);
constraint M[1,2] = B(0,1);
constraint forall(i,j in 2..n-1)(if i=j then M[i,i-1] = B(i-1,i-2) /\ M[i,i] = B(i-1,i-1) /\ M[i,i+1] = B(i-1,i) else if j<i-1 \/ j>i+1 then M[i,j]=0 endif endif);
constraint M[n,n-1] = B(n-1,n-2);
constraint M[n,n] = B(n-1,n-1);

constraint forall(i in 1..n)(L[i]=-30*ei(i-1,0));

constraint forall(i in 1..n)(sum(j in 1..n)(M[i,j]*result[j])=L[i]);

solve satisfy;

output ["M = \n"] ++ [show(M[i,j]) ++ if j=n then "\n" else " || " endif | i, j in 1..n] ++ ["\n"] ++
       ["L = " ++ show([L[i] | i in 1..n]) ++ ";\n"] ++ ["\n"] ++
       ["result = " ++ show([result[i] | i in 1..n]) ++ ";\n"]
