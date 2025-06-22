%% Pruebas varias

close all;
clearvars;
clc;

% Esta linea ya esta preparada para que le de usted y pueda seleccionar las
% carpetas de los experimentos (Si le pide formato de imagen,este  es tif ), 
% posteriormente le pedira las carpetas de los
% datos electricos

id = 0;
num_cons = 1;
ancho_cons = 7.5;
fs = 100;
Q = 0.05;
fondo = true;
girado = false;

device = Dispositivo(id,num_cons,ancho_cons,fs,Q,fondo,girado); 

% La siguiente linea le permite sacar todos los eventos de la prueba
device.find_events;

% La siguiente linea marca los eventos que si son eventos y los analiza
device.analyse_events;

% Con esta linea puede guardar todos los eventos en un txt
device.writeEventsTxt;

% Con esta linea puede guardar solo los eventos que son validos en un txt
device.writeRealEventsTxt()


%% Pruebas IA

evento_prueba = device.Canales.Eventos(4);
corte_imagen=evento_prueba.Canal.position_izquierda_superior;

frame = evento_prueba.frame_inicial;

n = frame + 1;
m = n - evento_prueba.frame_inicial;

load('resunet_propia_v7.mat')
img = imread(device.direc_imagenes(n));

figure
imshow(img)
hold on
linea = evento_prueba.Canal.linea;
xline(linea)
xline(evento_prueba.vector_Al(n-evento_prueba.frame_inicial,1))
hold off

img_resized = img(corte_imagen(1,2):corte_imagen(1,2)+511,corte_imagen(1,1):corte_imagen(1,1)+511);
imshow(img_resized)
C = semanticseg(img_resized, net_n);

B = labeloverlay(img_resized, C);
figure
imshow(B);
celula=C=="C3";

figure
imshow(celula)

[regiones,j]= bwlabel(celula,4);

propiedades=regionprops(regiones==1,"FilledArea");
areas=zeros(j,1);
areas(1)=propiedades.FilledArea;

for i=2:1:j
    region_n=regiones==i;
    propiedades=regionprops(region_n,"FilledArea");
    areas(i)=propiedades.FilledArea;
end

[valor,pos]=max(areas);

imshow(regiones==pos)
resaltada = img_resized;
resaltada(regiones == pos)=255;

linea_verdadera=linea+corte_imagen(1,1);

imshow(evento_prueba.Canal.Dispositivo.imagen_base_ori)
hold on
xline(linea_verdadera)

Al = evento_prueba.vector_Al(m,1);


imshowpair(device.imagen_base,img_resized)
imshowpair(img_resized,regiones==pos)

imshow(resaltada)
hold on
% xline(linea)
frente = xline(floor(evento_prueba.vector_Al(m,1)+linea),"g");
scatter(floor(evento_prueba.vector_Al(m,1))+linea,floor(evento_prueba.vector_Al(m,2)),"r.")
back = xline(floor(evento_prueba.vector_back(m,1)+linea),"b");
scatter(floor(evento_prueba.vector_back(m,1))+linea,floor(evento_prueba.vector_back(m,2)),"r.")
legend([frente, back],["Frente de la célula", "Parte trasera"])
