%Exercise7
%Christos Trimas 2016030054
%Piskopos Dionisis 2015030115


clear all;
close all;
clc;

e = 0.05;                                       %stathera gia to fragma

prompt = 'Enter n:';                            %dinoume xeirokinhta times gia to plithos n kai thn pithanothta
n = input(prompt);
prompt_1 = 'Enter p:';                  
p = input(prompt_1);

    if (n == 15) && (p == 0.45)                     %elegxoi periptwshs
        q = 1-p;
        x_axis = [-0.25*10^4 3.5*10^4];             %aksones se kathe periptwsh apo thn ekfwnhsh
    elseif (n == 15) && (p == 0.1)
        q = 1 - p;
        x_axis = [-0.25*10^4 3.5*10^4];
    elseif (n == 50) && (p == 0.45)
        q = 1-p;
        x_axis = [-10^14 12*10^14];
    elseif (n == 50) && (p == 0.1)
        q = 1 - p;
        x_axis = [-10^14 12*10^14];
    end;

H = -p*log2(p) -q*log2(q);                              %diadikh entropia
upper_bound = (2^(-n*(H-e)))*ones(1,2);                 %anw fragma apo thewria
lower_bound = (2^(-n*(H+e)))*ones(1,2);                 %katw fragma apo thewria
numOfTypical = (2^(-n*H))*ones(1,2);                    %plithos tupikwn 

y(1) = q^n;                                             %(p^0)*((1-p)^n)
x(1) = 1;

for k = 1:n
    y(2*k) = (p^(k))*(q^(n-(k)));
    y(2*k+1) = (p^(k))*(q^(n-(k)));
    x(2*k) = x(2*k-1)+1;
    x(2*k+1) = x(2*k) + nchoosek(n,k);
end;

figure;
semilogy(x,y);
hold on;
semilogy(x_axis,upper_bound,'g');
hold on;
semilogy(x_axis,lower_bound,'c');
hold on;
semilogy(x_axis,numOfTypical,'r');
xlabel('Num of Sequences');
ylabel('Probabilities');
title('Exercise 1-4');
