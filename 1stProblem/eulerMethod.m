###########################
#       Euler's Method
###########################

## Initialization
h = 0.1;
t = 0:h:600;

AM = (2631+3088+2977)/3;
Dx = 11835 + AM;
Dy = 11835 + AM;
Dz = 6000;
m = 425000;
mz = 357000000;
ma = 113000;
x(1) = 0;
y(1) = -AM/1000;
z(1) = 0;

###############   [fx,fy,nz] = [-AM,0,0]   #############
nz = 0;
fy = 0;
fx = -AM;


x = zeros(size(t));
z = zeros(size(t));
y = zeros(size(t));

n = numel(z);

u = zeros(size(t));
u(1) = 0;
w = zeros(size(t));
w(1) = 0;
k = zeros(size(t));
k(1) = 0;

##
for i = 1:(n-1)
  #Z (Ø)
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

figure(1);
subplot(3,1,1);
plot(t,x);
title("x(t)");
grid on;

subplot(3,1,2);
plot(t,y);
title("y(t)");
grid on;

subplot(3,1,3);
plot(t,z);
title("ø(t)");
grid on;


###############   [fx,fy,nz] = [0,-AM,0]   #############
nz = 0;
fy = -AM;
fx = 0;


x = zeros(size(t));
z = zeros(size(t));
y = zeros(size(t));

n = numel(z);

u = zeros(size(t));
u(1) = 0;
w = zeros(size(t));
w(1) = 0;
k = zeros(size(t));
k(1) = 0;

##
for i = 1:(n-1)
  #Z (Ø)
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

figure(2);
subplot(3,1,1);
plot(t,x);
title("x(t)");
grid on;

subplot(3,1,2);
plot(t,y);
title("y(t)");
grid on;

subplot(3,1,3);
plot(t,z);
title("ø(t)");
grid on;



###############   [fx,fy,nz] = [0,0,AM]   #############
nz = AM;
fy = 0;
fx = 0;


x = zeros(size(t));
z = zeros(size(t));
y = zeros(size(t));

n = numel(z);

u = zeros(size(t));
u(1) = 0;
w = zeros(size(t));
w(1) = 0;
k = zeros(size(t));
k(1) = 0;

##
for i = 1:(n-1)
  #Z (Ø)
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

figure(3);
subplot(3,1,1);
plot(t,x);
title("x(t)");
grid on;

subplot(3,1,2);
plot(t,y);
title("y(t)");
grid on;

subplot(3,1,3);
plot(t,z);
title("ø(t)");
grid on;
