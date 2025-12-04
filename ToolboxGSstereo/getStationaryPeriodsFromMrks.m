function [Act, S, V] = getStationaryPeriodsFromMrks(A, fs, varargin)
% 
%   Input: 
%       A   [Nx3]        Acceleration values X, Y, Z [g]
%
%   Optional arguments ('name',value):
%       actThresh   double  Threshold for activity (default: 0.005 [g])
%       wlen        double  Duration of epoch (detault: 1 [s])
%
%   Output:
%       A   [Nx1]       Activity vector: 1 --> activity detected
%                                        0 --> inactivity detected
%       S   [Nx4]       Samples time, acceletation values XYZ from
%                                     stationary periods (mean over epoch).
%       V   [Nx1]       Variance vector
%

p = inputParser;
    
% define optional arguments with default values
addOptional(p,'actThresh', 0.002, @isnumeric); % activity thresholds on variance    
addOptional(p,'wlen',          1, @isnumeric); % epoch duration

% parse varargin
parse(p,varargin{:});
p = p.Results;

% parse again 
actThresh =  p.actThresh;
wlen = p.wlen;
resampleRate = fs;

% initialise return variable
S = zeros(size(A,1),4);
Act = ones(size(A,1),1);
V = ones(size(A,1),1);

% go through data and get periods of stationarity
st = 1/fs;         
en = size(A,1)/fs-wlen;  
    
sCnt = 1;

% interpolate for speed later on
T = st:(1/resampleRate):en;

% turn from fraction of days to sample numbers
wlen = wlen * resampleRate;

%% AMVD
for t=1:length(T)-wlen        
    d = A(t:t+wlen-1,:);    
    mean_a = nanmean(d,1);    
    V_now = sum((d-mean_a).^2)/wlen;
    V(t) = norm(V_now);
    if V(t) <= actThresh %sum(V_now <= actThresh) >= 3
        % stationary!
        S(sCnt,:) = [t/fs mean(d(:,1:3)) ];
        Act(t,1) = zeros(1,1);
        sCnt = sCnt + 1;
    end
    if t == length(T)-wlen   
        if Act(t-1) == 0
            Act(t:end)=zeros(length(Act(t:end)),1);
        end
        V(t:end) = ones(length(V(t:end)),1)*V(t-1);
    end
end

S = S(1:sCnt-1,:);
end