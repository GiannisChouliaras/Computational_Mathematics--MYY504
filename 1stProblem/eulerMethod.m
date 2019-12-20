###########################
#       Euler's Method
###########################

## Initialization
h = 0.1;
t = 0:h:600;

nz = 0;
fy = 0;
fx = -AM;
AM = (2631+3088+2977)/3;
Dx = 11835 + AM;
Dy = 11835 + AM;
Dz = 6000;
m = 425000;
mz = 357000000;
ma = 113000;

x = zeros(size(t));
z = zeros(size(t));
y = zeros(size(t));
x(1) = 0;
y(1) = -AM/1000;
z(1) = 0;
n = numel(z);

u = zeros(size(t));
u(1) = 0;
w = zeros(size(t));
w(1) = 0;
k = zeros(size(t));
k(1) = 0;

## Z
for i = 1:(n-1)
  #Z
  df = (nz - 2 * Dz * abs(u(i)) * u(i))/mz;
  u(i+1) = u(i) + h*df;
  z(i+1) = z(i) + h*u(i);
  #X
  dg = (((fx - Dx * abs(w(i)) * w(i))/(m+3*ma)) + sin(z(i))*u(i)*w(i) - cos(z(i))*u(i)*w(i))/cos(z(i));
  w(i+1) = w(i) + h*dg;
  x(i+1) = x(i) + h*w(i);
  #Y
  dp = (((fy - Dy * abs(k(i)) * k(i))/(m+2*ma)) + cos(z(i))*u(i)*w(i) + sin(z(i))*u(i)*k(i))/cos(z(i));
  k(i+1) = k(i) + h*dp;
  y(i+1) = y(i) + h*k(i);
endfor

plot(t,x);
hold on
plot(t,y);
hold on
plot(t,z);
hold on


