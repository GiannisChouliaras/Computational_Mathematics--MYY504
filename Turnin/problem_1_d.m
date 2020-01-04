###########################
# Akd/ko etos: 2019-2020
# MYY504 Ypo\ka Mathimatika

# TEAM NUMBER 59 (B'PROJECT)
# Evaggelos Tzortzis   3088
# Chouliaras Ioannis   2631
# Zoganas Dimitrios    2977
###########################

## Initialization
h = 0.1;
t = 0:h:600;

AM = (2631+3088+2977)/3;
Dx = 11835;
Dy = 11835 + AM;
Dz = 6000;
m = 425000;
mz = 357000000;
ma = 113000;
y(1) = 0;
z(1) = 0;
Kpx = 60000;
Kdx = 5000000;
Kpy = 60000 - 5*AM;
Kdy = 5000000;
Kpz = 50000;
Kdz = 7000000 - 100*AM;
xdes = AM/100;
ydes = AM/100;
zdes = -AM/10000;


############# EULER'S METHOD ##############

x = zeros(size(t));
z = zeros(size(t));
y = zeros(size(t));
x(1) = -AM/1000;

n = numel(z);

u = zeros(size(t));
u(1) = 0;
w = zeros(size(t));
w(1) = 0;
k = zeros(size(t));
k(1) = 0;

for i = 1:(n-1)
  nz = Kpz*(zdes-z(i))-Kdz*u(i);
  fy = Kpy*(ydes-y(i))-Kdy*k(i);
  fx = Kpx*(xdes-x(i))-Kdx*w(i);
  #Z (psi)
  df = (nz - 2 * Dz * abs(u(i)) * u(i))/mz;
  u(i+1) = u(i) + h*df;
  z(i+1) = z(i) + h*u(i);
  #X
  dg = (((fx - Dx * abs(w(i)) * w(i))/(m+3*ma)) + sin(z(i))*u(i)*w(i) - cos(z(i))*u(i)*k(i))/cos(z(i));
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
title("\\psi(t)");
grid on;

############### IMPROVED EULER'S METHOD  #############

x = zeros(size(t));
y = zeros(size(t));
z = zeros(size(t));
x(1) = -AM/1000;

n = numel(z);

u = zeros(size(t));
u(1) = 0;
w = zeros(size(t));
w(1) = 0;
v = zeros(size(t));
v(1) = 0;

## may the force, be with you:
for i = 1:(n-1)
  nz = Kpz*(zdes-z(i))-Kdz*u(i);
  fy = Kpy*(ydes-y(i))-Kdy*k(i);
  fx = Kpx*(xdes-x(i))-Kdx*w(i);

  #every k1(z,x,y)
  k1z = u(i);
  k1x = w(i);
  k1y = v(i);

  #every k1 (u,w,v)
  k1u = (nz - (2 * Dz * abs(u(i)) * u(i))) / mz;
  k1w = (((fx - Dx * abs(w(i)) * w(i))/(m+3*ma)) + sin(z(i))*u(i)*w(i) - cos(z(i))*u(i)*v(i))/cos(z(i));
  k1v = (((fy - Dy * abs(v(i)) * v(i))/(m+2*ma)) + cos(z(i))*u(i)*w(i) + sin(z(i))*u(i)*v(i))/cos(z(i));

  #every k2(z,x,y)
  k2z = u(i) + (h * k1u);
  k2x = w(i) + (h * k1w);
  k2y = v(i) + (h * k1v);

  nz1 = Kpz*(zdes-(z(i)+h*k1z))-Kdz*k2z;
  fy1 = Kpy*(ydes-(y(i)+h*k1y))-Kdy*k2y;
  fx1 = Kpx*(xdes-(x(i)+h*k1x))-Kdx*k2x;

  #every k2 (u,w,v)
  k2u = (nz1 - (2*Dz*abs(k2z)*k2z))/mz;
  k2w = ((fx1 - Dx*abs(k2x)*k2x)/(m+3*ma) + sin(z(i)+h*k1z)*k2z*k2x - cos(z(i)+h*k1z)*(u(i)+h*k1u)*k2y)/cos(z(i)+h*k1z);
  k2v = ((fy1 - Dy*abs(k2y)*k2y)/(m+2*ma) + cos(z(i)+h*k1z)*k2z*k2x + sin(z(i)+h*k1z)*(u(i)+h*k1u)*k2y)/cos(z(i)+h*k1z);
  #Z
  u(i+1) = u(i) + (h/2) * (k1u + k2u);
  z(i+1) = z(i) + (h/2) * (k1z + k2z);

  #X
  w(i+1) = w(i) + (h/2) * (k1w + k2w);
  x(i+1) = x(i) + (h/2) * (k1x + k2x);

  #Y
  v(i+1) = v(i) + (h/2) * (k1v + k2v);
  y(i+1) = y(i) + (h/2) * (k1y + k2y);
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
title("\\psi(t)");
grid on;
