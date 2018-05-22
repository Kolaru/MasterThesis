t0 = cputime ;

n = 100000 ;
c = 3.4 ;

A = ER(n,c);
B = largestcomponent(A);
display(size(B, 2)) ;
e = cputime - t0 ;
e



