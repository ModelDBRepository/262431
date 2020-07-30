function Out = FunCm(Z,Cm0,a,delta)
if Z~=0
Out = ((Cm0*delta)/a^2)*(Z+((a^2-Z.^2-Z*delta)./(2*Z))...
    *log((2*Z+delta)/delta));
else
Out = Cm0;
end
end