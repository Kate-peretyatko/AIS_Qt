clear;
clc;
%***************** ОТКРЫТИЕ ФАЙЛА *****************************************
fid = fopen('ais138.dat', 'rb');  % открытие файла на чтение
if fid == -1
    error('File is not opened');
end

%f=fopen('my.txt','wt');

frame=0;             % инициализация переменной
cnt=1;              % инициализация счетчика
while ~feof(fid)    % цикл, пока не достигнут конец файла
    [V,N] = fread(fid, 1, 'float');  %считывание одного
% значения float (V содержит значение
% элемента, N – число считанных элементов)
    if N > 0        % если элемент был прочитан успешно, то
        frame(cnt)=V;% формируем вектор-строку из значений V     
        %fprintf(f, '%f', V);
        %fprintf(f, ' ');
        cnt=cnt+1;  % увеличиваем счетчик на 1
    end
end

%fclose(f);
% disp(data);        % отображение результата на экран
% plot(frame);
% length(frame)
fclose(fid);        % закрытие файла 

%***************** ПРЕОБРАЗОВАНИЕК ДИФФ. ВИДУ *****************************
%Передискретизация в 4 раза
re1 = frame(1: 8: end);
im1 = frame(2: 8: end);
re2 = frame(3: 8: end);
im2 = frame(4: 8: end);
re3 = frame(5: 8: end);
im3 = frame(6: 8: end);
re4 = frame(7: 8: end);
im4 = frame(8: 8: end);
signal1 = complex(re1, im1);
signal2 = complex(re2, im2);
signal3 = complex(re3, im3);
signal4 = complex(re4, im4);

length_signal = 4*length(signal1);
signal = zeros(1, length_signal);

k = 1 : 4 : length(signal);
i = 1 : length(signal1);

signal(k) = signal1(i);
signal(k+1) = signal2(i);
signal(k+2) = signal3(i);
signal(k+3) = signal4(i);

% scatterplot(signal)
l_s = length(signal);
diff_signal = zeros(1,l_s);

%***1ый ВАРИАНТ***%
row = 1 : length(signal) - 1;
% diff_signal(row) = signal(row) .* conj(signal(1));

%***2ой ВАРИАНТ***%
diff_signal(row) = signal(row) .* conj(signal(row+1));

% l_diff = length(diff_signal)
% scatterplot(diff_signal)

%***************** ПОДГОТОВКА ПРЕАМБУЛЫ ***********************************
preamble = [0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 1 1 1 1 1 0];
%preamble = [1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 0 1 1 1 1 1 1 0];
length_data = length(signal) - 1;
length_preamble = length(preamble);
  
twiddles = zeros(1, 17);
k = 1 : 17;
twiddles(k) = exp(1i * (k - 1) * pi/8);

synchro = [twiddles(4:-1:1), twiddles(2:5), twiddles(4:-1:1), twiddles(2:5), twiddles(4:-1:1), twiddles(2:5), twiddles(4:-1:1), twiddles(2:5), twiddles(4:-1:1), twiddles(2:5),twiddles(4:-1:1), twiddles(2:5), twiddles(4:-1:1), twiddles(2:5), twiddles(4:-1:1), twiddles(2:5), twiddles(4:-1:1), twiddles(2:5), twiddles(4:-1:1), twiddles(2:5), twiddles(4:-1:1), twiddles(2:5), twiddles(4:-1:1), twiddles(2:5), twiddles(4:-1:1), twiddles(2:17), twiddles(2:9), twiddles(8:-1:5)];
%synchro = [twiddles(2:5), twiddles(4:-1:1), twiddles(2:5), twiddles(4:-1:1), twiddles(2:5), twiddles(4:-1:1), twiddles(2:5), twiddles(4:-1:1), twiddles(2:5),twiddles(4:-1:1), twiddles(2:5), twiddles(4:-1:1), twiddles(2:5), twiddles(4:-1:1), twiddles(2:5), twiddles(4:-1:1), twiddles(2:5), twiddles(4:-1:1), twiddles(2:5), twiddles(4:-1:1), twiddles(2:5), twiddles(4:-1:1), twiddles(2:5), twiddles(4:-1:1), twiddles(17:-1:14), twiddles(15:16), twiddles(1:2), twiddles(3:14), twiddles(15:17), twiddles(2:6), twiddles(5:-1:2)];

l_syn = length(synchro);

%**************************** БЫСТРАЯ КОРРЕЛЯЦИЯ **************************
% Длина фрейма
length_frame = 1024;
% Длина пакета данных
length_data = 1024;
% Вспомгательная переменная
i = 1;
% Заполняем до конца преамбулу нулями
synchro1 = [synchro, zeros(1,length_frame-length(synchro))];

p = 1 : length_frame : l_s; 
kol_frame = length(p);

for k = 1 : length_frame : l_s;    
    k_corr(i,:) = ifft(fft(signal(k:(k+length_frame-1))).*conj(fft(synchro1)));
    [max_corr(i),id(i)] = max(k_corr(i,:));
    i = i +1;
end
k_corr2 = k_corr(:);
% plot(abs(k_corr2))
i = 1 : kol_frame;
norm_max_corr(i) = max_corr(i)./max(max(max_corr(i)));

[m,i] = max(norm_max_corr);
start_package = (i-1)*length_frame + id(i);
% Выделение пакета данных
data = signal(start_package : (start_package + length_data - 1));
%row = 1 : 1023;
%data(row) = data(row) ./ conj(data(row+1));
%************************** ДЕМОДУЛЯТОР ***********************************
sig1= [1 exp(1i*pi/8) exp(1i*pi/4) exp(1i*3*pi/8) 1i]; %[ 1 exp(1i*pi/4) 1i exp(1i*3*pi/4) 1]; 
sig2=conj(sig1);
k=1;

bits = uint8(zeros(1, length(data')/4));
for index=3:4:length(data)-5
x=[data(index),data(index+1),data(index+2),data(index+3),data(index+4)];
%Products-- solver outputs
pr1=0;
pr2=0;

%const=conj(data(index))/abs(data(index));
for step=1:5
    pr1=pr1+x(step)*sig1(step);
    pr2=pr2+x(step)*sig2(step);
end;
if real(pr1) > real(pr2)
        bits(k) = 0;
    else
        bits(k) = 1;
end;
k=k+1;
end;

% Для первых 32 бит
for y = 1 : 32
if bits(y) == bits(y+1)
    bits(y) = 1;
else 
    bits(y) = 0;
end;
end;

%***************************** ПАРСИНГ ************************************
coordinate_data = bits(33:200);
l_c = length(coordinate_data);

FCR_data = bits(201:216);
l_FCR = length(FCR_data);

%************************ ПОДСЧЁТ CRC *************************************
summ = bitxor(coordinate_data(153:168), FCR_data);

order = 0;
polynomial = 0;
m = l_FCR;
syms x;
for t = l_c : -1 : (l_c-l_FCR+1)
    polynomial = polynomial + (x^order)*summ(m);
    order = order + 1;
    m = m - 1;
end
for k = (l_c-l_FCR) : -1 : 1
    polynomial = polynomial + (x^order)*coordinate_data(k);
    order = order + 1;
end;

% order_FCR = 0;
% poly_FCR = 0;
% syms x;
% for e = l_FCR : -1 : 1
%     poly_FCR = poly_FCR + (x^order_FCR)*FCR_data(e);
%     order_FCR = order_FCR + 1;
% end;

CRC = 0;
order_CRC= 0;
for h = 16 : -1 : 1
    CRC = CRC + x^order_CRC;
    order_CRC = order_CRC + 1;
end;

vector_data = [1];
CRC_data = [1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1];
coordinate_data(153:168) = summ;
for q = 1 : 168
    vector_data(q) = coordinate_data(q);
end

[q,r] = deconv(vector_data, CRC_data);

%************************** ПЕРЕВОД ИЗ 2ой в 10ую *************************
%*** ID
id_data = coordinate_data(1:32);
ID = 0;
for u = 32 : -1 : 1
    ID = ID + 2^(id_data(u));
end;
ID

%*** latitude
latitude_data = coordinate_data(33:64);
latitude = 0;
for u = 32 : -1 : 1
    latitude = latitude + 2^(latitude_data(u));
end;
latitude

%*** longitude
longitude_data = coordinate_data(65:96);
longitude = 0;
for u = 32 : -1 : 1
    longitude = longitude + 2^(longitude_data(u));
end;
longitude

%*** speed
speed_data = coordinate_data(97:112);
speed = 0;
for u = 16 : -1 : 1
    speed = speed + 2^(speed_data(u));
end;
speed

%*** course
course_data = coordinate_data(113:128);
course = 0;
for u = 16 : -1 : 1
    course = course + 2^(course_data(u));
end;
course