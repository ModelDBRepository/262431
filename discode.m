function varargout = discode(odetype,ode,tspan,y0,options,disconts,varargin)
%DISCODE encapsulates the MATLAB ode-suite, allowing for discontinuities (disconts 1-D array)
% by breaking up the tspan around the disconts values 

% disconts : 1 D array of time points where the ode-problem is discontinuous
% odetype: string array with available ode solver

% Events (EventFunc) in [disconts-eps,disconts+eps] or [tspan(1),tspan(1)+eps] are not registered.

if nargin < 6
  disconts = [];
  if nargin < 5
    options = [];
    if nargin < 4
        y0 = [];
        if nargin < 3
            tspan = [];
            if nargin < 2
                error('Discode: not enough input arguments');
            end
        end
    end
  end
end

if ~isa(ode,'function_handle')
    error('discode: ode has to be a function handle');      % Possibility of ode a .m, .mex, ... file not supported 
end
if isempty(y0) || isempty(tspan)
    error(message('MATLAB:odearguments:TspanOrY0NotSupplied', odetype));
end

liveDisconts = disconts(tspan(1)<disconts&tspan(2)>=disconts);  % disconts(tspan(1) == disconts(1)) is skipped (see below)

if nargin < 3           % Not supported 
    sol = ode;
else
    extdata = struct('odefun',ode,'options',options); extdata.varargin = varargin;
    sol = struct('solver',odetype,'extdata',extdata,'x',tspan(1),'y',y0);
    sol.stats = struct('nsteps',0,'nfailed',0,'nfevals',0,'npds',0,'ndecomps',0,'nsolves',0);
   
% switch odetype
%  case 'ode113'
%   sol.idata = struct('klastvec',0,'psi2d',zeros(2,1),'phi3d',zeros(3,1)); 
%  case 'ode15i'
%   sol.idata = struct('kvec',0); 
%  case 'ode15s'
%   sol.idata = struct('kvec',0,'dif3d',zeros(3,1));
%  case 'ode23'          
%   sol.idata = struct('f3d',zeros(3,1)); 
%  case 'ode23s' 
%   sol.idata = struct('k1',0,'k2',0); 
%  case 'ode23t'
%   sol.idata = struct('z',0,'znew',0);   
%  case 'ode23tb'
%   sol.idata = struct('t2',0,'y2',0);
%  case 'ode45' 
%   sol.idata = struct('f3d',zeros(3,1));
% end
end

warning('off','MATLAB:odextend:SolutionAlreadyAvailable'); terminalEvent = 0; endDisc = 0;
sol.x = horzcat(sol.x,tspan(1)+eps(tspan(1)));
sol.y = horzcat(sol.y,sol.y(:,end));      % solextends needs at least two points. This also skips the special case tspan(1) = liveDisconts(1)
if ~isempty(liveDisconts)
if liveDisconts(end) == tspan(end)
endDisc = 1;
liveDisconts = liveDisconts(1:end-1);
end
for i = 1:length(liveDisconts)
sol = myOdextend(sol,[],liveDisconts(i)-eps(liveDisconts(i)));    % Special case: two disconts machine precision from eachother - odeextend will do nothing and second discont is skipped
if (sol.x(end) < liveDisconts(i)-eps(liveDisconts(i))) % Equal (==): integration till endpoint. Greater (>): two consecutive discontinuities or discont. following tspan(1) -> Not terminal.
terminalEvent = 1; break;
end
sol.x = horzcat(sol.x,liveDisconts(i)+eps(liveDisconts(i)));
sol.y = horzcat(sol.y,sol.y(:,end));
end
if (endDisc) && ~terminalEvent
sol = myOdextend(sol,[],tspan(end)-eps(tspan(end)));
if (sol.x(end) < tspan(end)-eps(tspan(end)))
terminalEvent = 1;  
else
sol.x = horzcat(sol.x,tspan(end));
sol.y = horzcat(sol.y,sol.y(:,end));
end
end
end
if (~endDisc&&~terminalEvent)
sol = myOdextend(sol,[],tspan(end));
end

switch nargout
    case 1, varargout = sol;
    case 2, varargout = {sol.x',sol.y'}; 
    case 3, varargout = {sol.x',sol.y',sol.xe',sol.ye',sol.ie'};
end

end