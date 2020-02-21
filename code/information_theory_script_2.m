%Exercise 7 5-7
%Christos Trimas 2016030054
%Piskopos Dionisis 2015030115


clear all;
close all;
clc

e = 0.05;                                       %stathera gia to fragma
n = 10:10:1000;                                 %dianisma ekfwnhshs
number_of_seq = 2.^n;                           %arithmos akolouthiwn mhkous n gia to erwthma 5

%n_1 = 15;                                      %den exei kapoia praktikh aksia,dothike gia na ta ksexorisoume
p_1 = 0.45;                                     %pithanotita 0.45
q_1 = 1-p_1;                                    %pithanotita 0.55
H_1 = -p_1*log2(p_1) -q_1*log2(q_1);            %entropia

%n_2 = 50;                                      %den exei kapoia praktikh aksia,dothike gia na ta ksexorisoume
p_2 = 0.1;                                      %pithanotita 0.1
q_2 = 1 - p_2;                                  %pithanotita 0.9
H_2 = -p_2*log2(p_2) -q_2*log2(q_2);            %entropia

probability_e_1 = (1-e)*ones(1,100);             %erwthma 6 pithanotita 1-e
entropy_1 = H_1*ones(1,100);                     %erwthma 7 entropia
bound_entropy_1 = (H_1+e)*ones(1,100);           %fragma entropia + e
bound_entropy_2 = (H_1-e)*ones(1,100);           %fragma entropia - e
entropy_2 = H_2*ones(1,100);                     %erwthma 7 entropia 
bound_entropy_3 = (H_2+e)*ones(1,100);           %fragma entropia + e
bound_entropy_4 = (H_2-e)*ones(1,100);           %fragma entropia - e


temp_1=0;                                        %voithitikes metavlhtes
sum_1=zeros(100);
prob_1=zeros(1,100);
prob_1_6 = zeros(1,100);
num_OF_bits_1 = zeros(1,100);

temp_2 = 0;
sum_2 = zeros(100);
prob_2 = zeros(100);
prob_2_6 = zeros(1,100);
num_Of_bits_2 = zeros(1,100);

for n_1 = n
    upper_bound_A_1(n_1/10) = 2^(n_1*(H_1+e));              %anw fragma gia to plithos twn tupikwn
    lower_bound_A_2(n_1/10) = (1-e)*(2^(n_1*(H_1-e)));      %katw fragma gia twn tupikwn
    prob_1(n_1/10) = 1/(2^(n_1*H_1));                       %pithanotita tupikhs akolouthias
    
    for k = 0:n_1           
       if(2^(-n_1*(H_1+e)) <= (p_1^k)*(q_1^(n_1-k)) && (p_1^k)*(q_1^(n_1-k)) <= 2^(-n_1*(H_1-e)))   %elexos an einai tupikh
            temp_1 = nchoosek(n_1,k);                       %anathesi sth metavlhth gia epanaxrhsimopoihsh parakatw
            sum_1(n_1/10) = sum_1(n_1/10) + temp_1;         %athroisma tupikwn akolouthiwn, prohsoumeno + n ana k
            prob_1_6(n_1/10) = prob_1_6(n_1/10) + temp_1*(p_1^k)*(q_1^(n_1-k));     %pithanotita tupikhs akolouthias,prohgoumenh + n ana k epi p^k * q^n-k
       end
    end
     num_OF_bits_1(n_1/10) = log((sum_1(n_1/10))/n_1);   %kati paei strava kai den leitourgei sthn grafikh
end

for n_2 = n
    upper_A_2(n_2/10) = 2^(n_2*(H_2+e));
    lower_A_2(n_2/10) = (1-e)*(2^(n_2*(H_2-e))); 
    prob_2(n_2/10) = 1/(2^(n_2*H_2));
    
    for k = 0:n_2
       if(2^(-n_2*(H_2+e)) <= (p_2^k)*(q_2^(n_2-k)) && (p_2^k)*(q_2^(n_2-k)) <= 2^(-n_2*(H_2-e)))
            temp_2 = nchoosek(n_2,k);
            sum_2(n_2/10) = sum_2(n_2/10) + temp_2;
            prob_2_6(n_2/10) = prob_2_6(n_2/10) + temp_2*(p_2^k)*(q_2^(n_2-k));
       end
    end
    num_Of_bits_2(n_2/10) = (sum_2(n_2/10))/n_2;
end

figure;                            %erwthma 5
semilogy(n, sum_1,'b--')           %plithos tupikwn gia p = 0.45
hold on;
semilogy(n,upper_bound_A_1,'r');   %anw fragma
hold on;
semilogy(n,lower_bound_A_2,'g');   %katw fragma
hold on;
semilogy(n,number_of_seq,'b');     %plithos olwn twn diadikwn akolouthiwn
hold on;
semilogy(n,sum_2,'b');             %plithos tupikwn gia p = 0.1
hold on;
semilogy(n,upper_A_2,'r');         %anw fragma
hold on;
semilogy(n,lower_A_2,'g');         %katw fragma
hold on;
semilogy(n,number_of_seq,'c--');   %plithos olwn twn diadikwn akolouthiwn
xlabel('n');
ylabel('Number of sequences');
title('Erwthma 5');
 
figure;                            %axonas x ws pros n, erwthma 6
plot(n,prob_1_6);                  %pithanotita mia akolouthia einai tupikh gia p = 0.45
hold on;
plot(n,prob_2_6,'k');              %pithanotita mia akolouthia einai tupikh gia p=0.1
hold on;
plot(n,probability_e_1,'r');       %1-e
xlabel('n');
ylabel('Probability');
title('Erwthma 6');
 

figure;
semilogy(n,num_OF_bits_1,'b');              %plithos bits
hold on;
semilogy(n,entropy_1,'b');                  %entropia
hold on;
semilogy(n,bound_entropy_1,'b--');          %fragma_entropias
hold on;
semilogy(n,bound_entropy_2,'b-.');          %fragma_entropias
hold on;
semilogy(n,num_Of_bits_2,'k');              %plithos bits
hold on;
semilogy(n,entropy_2,'k');                  %entropia
hold on;
semilogy(n,bound_entropy_3,'k--');          %fragma entropias
hold on;
semilogy(n,bound_entropy_4,'k-.');          %fragma entropias
xlabel('n');
ylabel('bits/symbol');
title('Erwthma 7');
