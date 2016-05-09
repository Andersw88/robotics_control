function [ w ] = learnDMP(d, p)

% if ~isfield(p,'scale') p.scale = p.end - p.start; end
p.scale2 = p.end - p.start;

pS = p.a/p.tau^2*(p.end-d(:,1)) - p.b/p.tau*d(:,2)- d(:,3);
f = zeros(length(p.s),p.n);
for i = 1:p.n
    f(:,i) = gaussmf(p.s, [0.5/p.n (i-1)/(p.n-1)]);
end
f2 = sum(f,2);
G = p.s.*(p.scale2)./(p.tau^2*f2);
G = bsxfun(@times,G,f);

w = G\-pS;
end

