function ks = multinomial_many_draw( ps, N )
cs = cumsum( ps );
rs = rand(N,1) * cs(end );
ks = zeros(1,N);
for n = 1:N
    ks(n) = find( cs > rs(n), 1);
end