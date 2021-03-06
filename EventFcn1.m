function [value,isterminal,direction] = EventFcn1(t,Z,dZ,na,freq,dt,CORRTres) %#ok<INUSL>
isterminal = 1; direction = 1; value = -1;
global temp;
temp = vertcat(temp,[t Z na]);
if (t-temp(1,1))>=2/freq
    % 2 periods of data -> check for periodicity (if periodic -> value=1)
firstIn = find((t-temp(:,1))>=(2/freq),1,'last');
tempor = temp(firstIn:end,:);     % Extract last two periods for analysis
TimeLine = tempor(:,1)-tempor(1,1); % Shift zero in timeline to first value of (Z,na)
% Resample the Z and na values, to allow calculation of autocorrelation
% Resampling is necessary, because the VSVO-ode113 solver has variable
% stepsizes, rendering the autocorrelation meaningless
% + the VSVO solver has the habit to regress locally and slightly in time,
% this behaviour vanishes after resampling
SampleT = (0:dt:2/freq)'; 
Ztemp = interp1(TimeLine,tempor(:,2),SampleT,'spline','extrap');
natemp = interp1(TimeLine,tempor(:,3),SampleT,'spline','extrap');
% Calculate the unbiased normalized autocorrelation
[acorrZ,~] = xcorr(Ztemp-mean(Ztemp),'unbiased');
[acorrna,lags] = xcorr(natemp-mean(natemp),'unbiased');
lags = lags*dt;
acorrZ = acorrZ/acorrZ(ceil(length(acorrZ)/2));  % Normalize autocorrelation to 1
acorrna = acorrna/acorrna(ceil(length(acorrna)/2));
if acorrZ(find(lags>=1/freq,1))>=CORRTres&&acorrna(find(lags>=1/freq,1))>=CORRTres
value = 1;
end
end
end