function [c1 c2 c3] = speedinwater(varargin)
if nargin == 1
    T = varargin{1};
	S = 0;
else
    T = varargin{1};
    S = varargin{2};  % parts per thousand
end

% http://resource.npl.co.uk/acoustics/techguides/soundpurewater/content.html#LUBBERS
% range of validity: 10-40 C at atmospheric pressure
c1 = 1405.03 + 4.624*T - 3.83e-2*T.^2;

% Below is only for sea water...
% Coppens (1981) for sea water as a function of T, S, Depth(=0 here)
% Range of validity: T 0 to 35 癈, S 0 to 45 parts per thousand, depth 0 to 4000 m
% http://resource.npl.co.uk/acoustics/techguides/soundseawater/content.html#UNESCO
T1 = T/10;
c2 = 1449.05 + 45.7*T1 - 5.21*T1.^2 + 0.23*T1.^3 + (1.333 - 0.126*T1 + 0.009*T1.^2)*(S - 35);

% NaCl solution:
% <Solar Energy> DEPENDENCE OF SPEED OF SOUND ON SALINITY AND TEMPERATURE IN CONCENTRATED NaCl SOLUTIONS
% speed in NaCl solutions, T = 7-88 C, S in percentage; S<21%
alpha1 = 1403.09 + 4.68391*T - 0.0405388*T.^2 + 1.2955e-5*T.^3 + 6.914851e-7*T.^4;
alpha2 = 14.019 - 0.114996*T + 2.23748e-5*T.^2 + 1.48238e-5*T.^3 - 9.46165e-8*T.^4;
c3 = alpha1 + alpha2*S/10;

end
