function [hy,y] = hilbertization();

%  
%
%  y = cos(t)
%
%  generate time series
t = (0 : pi/32 : 4 * pi);
y = cos(t);

hy = hilbert(y);     % y + i H(y)  ->  cos(t) + i sin(t) -> exp(it)

plot(t,real(hy),t,imag(hy));


end
