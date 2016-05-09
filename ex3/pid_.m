function [ u ] = pid_( e,dt,G)
    persistent I pe;
    if isempty(I)
       I = 0;
    end
    if isempty(pe)
       pe = e;
       u = 0;
       return
    end
    I = I+e;
    de = (e-pe)/dt
    u = G(1)*e + I*G(2) + de*G(3);
    pe = e;
end

