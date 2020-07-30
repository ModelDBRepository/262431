function Out = SimplNICEThRT_nanoMC(ESi,DISPLAY,tNICE,t,Q1,h1,r1,Q2,h2,r2,CmR,Gna,Vna,Gk,Vk,Gl,Vl,...
    GT,VT,hinf,tauh,rinf,taur,minf,pinf,Cm0,aBLS,fBLS,RSI,proteinMode,gateMultip)
if DISPLAY == 2
global reverseStr; %#ok<TLEV>
Progress = 100*(t-tNICE(1))/(tNICE(2)-tNICE(1));  %#ok<*NASGU>
msg = sprintf('Progress: %3.1f', Progress); 
fprintf([reverseStr, msg]);
reverseStr = repmat(sprintf('\b'), 1, length(msg));  
end
switch proteinMode
    case 0, xP = 1; xl = 1;             % (ratio of protein coverage (electrolytes, leak) in the BLS compartment)
            MP = 1; Ml = 1;             % Multipliers of the gate and leakage currents
    case 1, xP = 0; xl = 1;
            MP = gateMultip; Ml = 1;
    case 2, xP = 0; xl = 0;
            MP = gateMultip; Ml = gateMultip;
end

h1 = h1*(h1<=1)+(h1>1); h2 = h2*(h2<=1)+(h2>1);
r1 = r1*(r1<=1)+(r1>1); r2 = r2*(r2<=1)+(r2>1);

V1 = 10^(3)*Q1/CmR(t); V2 = 10^(3)*Q2/Cm0;

Out1 = [ESi(t)+10^(-3)*(1/(pi*aBLS^2*RSI))*(V2-V1)-10^(-3)*(xl*Gl*(V1-Vl)+xP*Gna*minf(V1).^3.*h1.*(V1-Vna)+...
    xP*Gk*(0.75*(1-h1)).^4.*(V1-Vk)+xP*GT*pinf(V1).^2.*r1.*(V1-VT));
(hinf(V1)-h1)/tauh(V1);
(rinf(V1)-r1)/taur(V1)];
Out2 = [ESi(t)+10^(-3)*(1/(pi*aBLS^2*(1/fBLS-1)*RSI))*(V1-V2)-10^(-3)*(Ml*Gl*(V2-Vl)+MP*Gna*minf(V2).^3.*h2.*(V2-Vna)+...
    MP*Gk*(0.75*(1-h2)).^4.*(V2-Vk)+MP*GT*pinf(V2).^2.*r2.*(V2-VT));
(hinf(V2)-h2)/tauh(V2);
(rinf(V2)-r2)/taur(V2)];

Out = [Out1;Out2];
end