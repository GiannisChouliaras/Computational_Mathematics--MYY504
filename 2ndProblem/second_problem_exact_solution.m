

#Initialization
h = 0.1;
t = 0:h:600;

Kpz = 50000;
AM = (2631+3088+2977)/3;
Kdz = 7000000 - 100*AM;
Dz = 6000;
mz = 357000000;
zdes = -AM/10000;

#########################
#    Exact Solution     #
#########################

discr = (Kdz+2*Dz).^2-4*mz*Kpz;
vx = -(Kdz+2*Dz)/(2*mz);
vy = sqrt(abs(discr))/(2*mz);
c2 = zdes*(Kdz+2*Dz)/sqrt(abs(discr));

zt = -zdes.*exp(vx*t).*cos(vy*t) + c2.*exp(vx*t).*sin(vy*t) + zdes;

figure(1);
plot(t,zt);
title("Exact Solution for \\psi");

#########################
#    Euler's Method     #
#########################

z = zeros(size(t));
z(1) = 0;
n = numel(z);

u = zeros(size(t));
u(1) = 0;

for i = 1:(n-1)
  nz = Kpz*(zdes-z(i))-Kdz*u(i);
  #Z (psi)
  df = (nz-2*Dz*u(i))/mz;
  u(i+1) = u(i) + h*df;
  z(i+1) = z(i) + h*u(i);
endfor

figure(2);
plot(t,z);
title("Euler's Method for \\psi");

##################################
#    Improved Euler's Method     #
##################################

z1 = zeros(size(t));
z1(1) = 0;
n = numel(z1);

u1 = zeros(size(t));
u1(1) = 0;

for i = 1:(n-1)
  nz = Kpz*(zdes-z1(i))-Kdz*u1(i);
  
  k1z = u1(i);
  k1u = (nz-(2*Dz*u1(i)))/mz;
  k2z = u1(i)+(h*k1u);
  
  nz1 = Kpz*(zdes-(z1(i)+h*k1z))-Kdz*k2z;
  
  k2u = (nz1 - (2*Dz*k2z))/mz;
  
  u1(i+1) = u1(i) + (h/2) * (k1u + k2u);
  z1(i+1) = z1(i) + (h/2) * (k1z + k2z);
endfor

figure(3);
plot(t,z1);
title("Improved Euler's Method for \\psi");

#Difference between Exact Solution and  Euler's Method
difz = zeros(size(t));
for i = 1:n
  difz(i) = zt(i)-z1(i);
endfor
figure(4);   
plot(t,difz);
title("Difference between Exact Solution and Euler's Method for \\psi");

#Difference between Exact Solution and Improved Euler's Method
difz = zeros(size(t));
for i = 1:n
  difz(i) = zt(i)-z(i);
endfor
figure(5);   
plot(t,difz);
title("Difference between Exact Solution and Improved Euler's Method for \\psi");

### Extra ###
#Difference between Euler's Method and Improved Euler's Method
#{
difz = zeros(size(t));
for i = 1:n
  difz(i) = z(i)-z1(i);
endfor
figure(6);        
plot(t,difz);  
title("Difference between Euler's Method and Improved Euler's Method for \\psi");
#}
