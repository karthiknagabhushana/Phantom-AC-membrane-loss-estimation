function [arr_out,t_out] = snip_array(varargin)
arr_in = varargin{1};
t_in = varargin{2};
t1 = varargin{3};
t2 = varargin{4};

t_out = t_in((t_in > t1) & (t_in < t2));
arr_out = arr_in((t_in > t1) & (t_in < t2));
end

