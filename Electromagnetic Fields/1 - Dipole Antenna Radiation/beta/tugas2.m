%%TUGAS 2
clear all;
clc;
%%Masukan Variabel
eta = 377;
f = 3e9;
c = 3e8;
lambda = c/f;
I0 = 1;
r = 10*lambda;
theta = 0:0.01:2*pi;
k = 2*pi/lambda;
l = 1.5*lambda;
%%Menghitung Medan Radiasi1
var1 = i*eta*I0*exp(-i*k*r);
var2 = 2*(pi)*r;
var3 = cos((k*l/2)*cos(theta))-cos(k*l/2);
var4 = sin(theta);
Etheta = abs((var1/var2)*(var3./var4));
EthetadB = 10*log10(Etheta);
n = 256;
d = lambda/4;
BETA = -3*pi/2;
%%Menghitung Medan Radiasi2
psi = (k*d*cos(theta))+BETA;
var21 = 1/n;
var22 = sin((n/2)*psi);
var23 = sin((1/2)*psi);
af = abs((var21)*(var22./var23));
afdb = 10*log10(af);
Etheta2 = EthetadB.*afdb;

%%Membuat Grafik 2D
figure(3);
polarplot(theta,Etheta2);
title('2D Array Pola Radiasi Antenna Dipole: n=256');

%%Mapping Koordinat
azimuth = 0:0.01:2*pi;
nEtheta2 = Etheta2 - min(Etheta2);
nEtheta2(1,1) = 0;
matrixaz2 = [];
matrixth2 = [];
matrixet2 = [];
for i = 1:629
    matrixaz2(i,:) = azimuth(1,:);
    matrixth2(:,i) = theta(1,:);
    matrixet2(:,i) = nEtheta2(1,:);
end

x2 = matrixet2.*cos(matrixth2).*cos(matrixaz2);
y2 = matrixet2.*cos(matrixth2).*sin(matrixaz2);
z2 = matrixet2.*sin(matrixth2);
r_theta = matrixet2;
%%Membuat Grafik 3D
figure(4);
h2 = surf(x2,y2,z2,r_theta,'EdgeColor','none','FaceColor','interp');
cm=colormap; cm=flipud(cm); colormap(cm);
colorbar;
title('3D Array Pola Radiasi Antenna Dipole: n=256');