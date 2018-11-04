function out = gauss(x,s)
out = exp(-x.*x/(2*s*s))/sqrt(2*pi*s*s);
