function [ ds ] = canonicalSystem(s,p)
    ds = -s*p.g/p.tau;
end

