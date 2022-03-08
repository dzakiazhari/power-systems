%%TUGAS 1
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
l = 5*lambda;

%%Menghitung Medan Radiasi
var1 = i*eta*I0*exp(-i*k*r);
var2 = 2*(pi)*r;
var3 = cos((k*l/2)*cos(theta))-cos(k*l/2);
var4 = sin(theta);
Etheta = abs((var1/var2)*(var3./var4));
EthetadB = 10*log10(Etheta);

%%Membuat Grafik 2D
figure(1);
polarplot(theta,EthetadB);
title('2D Pola Radiasi Antenna Dipole: l=5 lambda');

%%Mapping Koordinat
azimuth = 0:0.01:2*pi;
nEthetadB = EthetadB - min(EthetadB);
nEthetadB(1,1) = 0;
matrixaz1 = [];
matrixth1 = [];
matrixet1 = [];
for i = 1:629
    matrixaz1(i,:) = azimuth(1,:);
    matrixth1(:,i) = theta(1,:);
    matrixet1(:,i) = nEthetadB(1,:);
end

x1 = matrixet1.*cos(matrixth1).*cos(matrixaz1);
y1 = matrixet1.*cos(matrixth1).*sin(matrixaz1);
z1 = matrixet1.*sin(matrixth1);

%%Membuat Grafik 3D
figure(2);
r_theta = matrixet1;
h1 = surf(x1,y1,z1,r_theta,'EdgeColor','none','FaceColor','interp');
cm=colormap; cm=flipud(cm); colormap(cm);
colorbar;
title('3D Pola Radiasi Antenna Dipole: l=5 lambda');