function [ f ] = forcingFunction( s,p )

if isfield(p,'scale')
    scale = p.scale; 
else
    scale = 1;
end

f1 = 0;
f2 = 0;
for i = 1:p.n
    f1 = f1 + gaussmf(s, [0.5/p.n (i-1)/(p.n-1)]);
    f2 = f2 + gaussmf(s, [0.5/p.n (i-1)/(p.n-1)])*p.w(i);
end
f = s.*(p.end - p.start)*scale.*f2./(f1);
end

