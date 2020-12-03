clear all;
close all;
tic
M = 16;
b = log2(M);
n = 7;k = 4;
d = 1;
%N = 4.48*10^6; %total bits 
N = 1.12*10^6;
input = randi([0,1],1,N);
% reshape the input signal with 4 a group, so that a column represents
% a symbol(4 bits)
pre_hamming = reshape(input, [4, N/4]); %raw matrix containing only signal columns
uncoded_symbols = QAM_Map_Symbol(pre_hamming, d, N);
% linear block code
identity = diag([1 1 1 1]);
Parity = [1,1,0;0,1,1;1,1,1;1,0,1];
Generator = [Parity,identity];

message = mod((Generator')*pre_hamming, 2); %m--[7*N]

l = length(message(:));
pre_map = reshape(message, [1, l]);

input_qua = reshape(pre_map, [4, l/4]);
E_s4=(4^2-1)*(d.^2)/12; 
Es=2*E_s4;

EsN0 = -2:1:25;
N0 = Es./(10.^(EsN0./10));  
EbN0_db = EsN0 - 10*log10(4);
EbN0 = 10.^(EbN0_db./10); 
std_dev =reshape(sqrt(N0./2),1,1,length(EsN0)); %standard deviation

s1 = [-3*d/2; 3*d/2];s2 = [-d/2; 3*d/2];s3 = [d/2; 3*d/2];s4 = [3*d/2; 3*d/2];s5 = [-3*d/2; d/2];s6 = [-d/2; d/2];s7 = [d/2; d/2];
s8 = [3*d/2; d/2];s9 = [-3*d/2; -d/2];s10 = [-d/2; -d/2];s11 = [d/2; -d/2];s12 = [3*d/2; -d/2];s13 = [-3*d/2; -3*d/2];s14 = [-d/2; -3*d/2];
s15 = [d/2; -3*d/2];s16 = [3*d/2; -3*d/2];

bit1 = [0,0,1,0].'; bit2 = [0,1,1,0].'; bit3 = [1,1,1,0].'; bit4 = [1,0,1,0].'; bit5 = [0,0,1,1].'; bit6 = [0,1,1,1].';
bit7 = [1,1,1,1].'; bit8 = [1,0,1,1].'; bit9 = [0,0,0,1].'; bit10 = [0,1,0,1].'; bit11 = [1,1,0,1].'; bit12 = [1,0,0,1].'; 
bit13 = [0,0,0,0].'; bit14 = [0,1,0,0].'; bit15 = [1,1,0,0].'; bit16 = [1,0,0,0].'; 

signals = QAM_Map_Symbol(input_qua, d, l);

% ------------------------------------------------------------------------
% ------------ Add gaussian white noise part
%-------------------------------------------------------------------------
% build a 3D array, where:
% row is different index and for different noise addition, namely 
% different N0
noise =std_dev.*randn(2,length(signals),length(EsN0));

% AWGN
received = signals + noise;

% 16 signals ------ 16 decision condition
s1_shot = received(1,:,:)<-d & received(2,:,:)>d;
s2_shot = received(1,:,:)>-d & received(1,:,:)<0 & received(2,:,:)>d;
s3_shot = received(1,:,:)>0 & received(1,:,:)<d & received(2,:,:)>d;
s4_shot = received(1,:,:)>d & received(2,:,:)>d;
s5_shot = received(1,:,:)<-d & received(2,:,:)<d & received(2,:,:)>0;
s6_shot = received(1,:,:)<0 & received(1,:,:)>-d & received(2,:,:)<d & received(2,:,:)>0;
s7_shot = received(1,:,:)>0 & received(1,:,:)<d & received(2,:,:)<d & received(2,:,:)>0;
s8_shot = received(1,:,:)>d & received(2,:,:)<d & received(2,:,:)>0;
s9_shot = received(1,:,:)<-d & received(2,:,:)<0 & received(2,:,:)>-d;
s10_shot = received(1,:,:)<0 & received(1,:,:)>-d & received(2,:,:)<0 & received(2,:,:)>-d;
s11_shot = received(1,:,:)>0 & received(1,:,:)<d & received(2,:,:)<0 & received(2,:,:)>-d;
s12_shot = received(1,:,:)>d & received(2,:,:)<0 & received(2,:,:)>-d;
s13_shot = received(1,:,:)<-d & received(2,:,:)<-d;
s14_shot = received(1,:,:)>-d & received(1,:,:)<0 & received(2,:,:)<-d;
s15_shot = received(1,:,:)>0 & received(1,:,:)<d & received(2,:,:)<-d;
s16_shot = received(1,:,:)>d & received(2,:,:)<-d;

dm_s1 = s1_shot.*s1; dm_s2 = s2_shot.*s2; dm_s3 = s3_shot.*s3;
dm_s4 = s4_shot.*s4; dm_s5 = s5_shot.*s5; dm_s6 = s6_shot.*s6;
dm_s7 = s7_shot.*s7; dm_s8 = s8_shot.*s8; dm_s9 = s9_shot.*s9;
dm_s10 = s10_shot.*s10; dm_s11 = s11_shot.*s11; dm_s12 = s12_shot.*s12;
dm_s13 = s13_shot.*s13; dm_s14 = s14_shot.*s14; dm_s15 = s15_shot.*s15;
dm_s16 = s16_shot.*s16;

demodulated = dm_s1 + dm_s2 + dm_s3 + dm_s4 + dm_s5 + dm_s6 + dm_s7 + dm_s8 + dm_s9 + dm_s10 + dm_s11 + dm_s12 + dm_s13 + dm_s14 + dm_s15 + dm_s16;
demodulated_bits = s1_shot.*bit1 + s2_shot.*bit2 + s3_shot.*bit3 + s4_shot.*bit4 + s5_shot.*bit5 + s6_shot.*bit6 + s7_shot.*bit7 + s8_shot.*bit8 + s9_shot.*bit9 + s10_shot.*bit10 + s11_shot.*bit11 + s12_shot.*bit12 + s13_shot.*bit13 + s14_shot.*bit14 + s15_shot.*bit15 + s16_shot.*bit16;
%operation to convert the demodulated bits back to n(seven bits) a group;
%the first four columns are data bits, the last three columns are parity
%check bits
dm_7bits = zeros([l/7, 7, length(EsN0)]);
demodulated_bits_rows = zeros(1,l,length(EsN0));
Syndromes = zeros([l/7, n-k, length(EsN0)]);
decimal_syndromes = zeros([l/7, 1, length(EsN0)]);
for k=1:length(EsN0)
    demodulated_bits_rows(:,:,k) = reshape(demodulated_bits(:,:,k), [1, l]); % [1*l*33]
end
% size(demodulated_bits) -------4*490000*33   490000 is l/4
for i=1:length(EsN0)
    dm_7bits(:,:,i) = reshape(demodulated_bits_rows(:,:,i), [7, l/7])'; % size: [7*(l/7)*33]
end
H = [diag([1 1 1]), Parity'];
H = repmat(H', [1 1 length(EsN0)]);
for j = 1:length(EsN0)
    Syndromes(:,:,j) = mod(dm_7bits(:,:,j)*H(:,:,j),2); 
end

for i=1:length(EsN0)
    decimal_syndromes(:,:,i) = 4*Syndromes(:,1,i) + 2*Syndromes(:,2,i) + Syndromes (:,3,i);
end

syns_2D = [decimal_syndromes(:,:)]; %colums conrrespond to different signals with different noises(total 33)

%use the number in syns_2d as index in this matrix, to map the corresponding syndromes, 
%and then xor it with signals(7 bit)to correct
syndrome_pattern_matrix = [
    0,0,1,0,0,0,0;
    0,1,0,0,0,0,0;
    0,0,0,0,1,0,0;
    1,0,0,0,0,0,0;
    0,0,0,0,0,0,1;
    0,0,0,1,0,0,0;
    0,0,0,0,0,1,0
];
% error correction
for types = 1:7
    for pages = 1:length(EsN0)
        error_shots = find((syns_2D(:,pages)==types)); % errors' shots that correspond to the syndromes
        temp = mod(dm_7bits(error_shots, :, pages)+syndrome_pattern_matrix(types,:), 2);
        dm_7bits(error_shots, :, pages) = temp;
    end
end
toc

% ------------------------------------------------------------------------
% ------------ Evaluation part
%-------------------------------------------------------------------------
P4 = (2*(4-1)/4).*Q(d./sqrt(2.*N0));        
noECC_Theoretical_SER = P4+P4-P4.*P4;

noECC_Theoretical_BER = (4/4)*(1-1/sqrt(M)).*Q(sqrt(3*4*EbN0/(M-1)));

is_equal = (demodulated == signals);
symbol_error =  reshape(is_equal(1,:,:) & is_equal(2,:,:),length(received),1,length(EsN0));
symbol_error_num = squeeze(length(received) - sum(symbol_error));
figure(1);
semilogy(EsN0,symbol_error_num./(length(signals)),'mx--');

hold on 
semilogy(EsN0,noECC_Theoretical_SER,'b.-');

%bit error rate without ECC code
is_equal = reshape(demodulated_bits == input_qua,length(signals),4,length(EsN0));
% size(is_equal)
bit_error_num = squeeze(l - (sum(is_equal(:,1,:))+sum(is_equal(:,2,:))+sum(is_equal(:,3,:))+sum(is_equal(:,4,:))));

%transpose the 3D matrix to more easily compare with the original signals
dm_7bits_transpose = zeros(7, l/7, length(EsN0));
for pages = 1:length(EsN0)
    dm_7bits_transpose(:,:,pages) = dm_7bits(:,:,pages)';
end
%convert the corrected bits back to symbols to evaluate SER
for pages=1:length(EsN0)
    corrected_symbols(:,:,pages) = QAM_Map_Symbol(dm_7bits_transpose(4:7,:,pages), d, 4*l/7);
end
is_corrected_equal = (corrected_symbols == uncoded_symbols);
symbol_error =  reshape(is_corrected_equal(1,:,:) & is_corrected_equal(2,:,:),length(corrected_symbols),1,length(EsN0));
symbol_error_num = squeeze(length(corrected_symbols) - sum(symbol_error));
hold on 
semilogy(EsN0,symbol_error_num./(length(uncoded_symbols)),'*');

%Theoretical Symbol error rate when using hamming code
ser_ecc=zeros(1,length(EsN0));
for j=2:n
      ser_ecc = ser_ecc + nchoosek(n,j).*(noECC_Theoretical_BER.^j).*((1-noECC_Theoretical_BER).^(n-j));
end

ber_ecc=zeros(1,length(EsN0));
for j=2:n
    ber_ecc = ber_ecc + j*nchoosek(n,j).*(noECC_Theoretical_BER.^j).*((1-noECC_Theoretical_BER).^(n-j))/7;
end

hold on 
semilogy(EsN0,ser_ecc,'--');
xlabel('Es/N0 (dB)');
ylabel('Symbol Error Rate');
legend('SER(no ECC)', 'theoretical(no ECC)', 'SER(with ECC)', 'theoretical(with ECC)');
title('SER');

figure(2);
semilogy(EbN0_db,bit_error_num./l);
line = ones(length(EsN0),1)./10^4;
% semilogy(ECC_ber,line,'mx--');
hold on
semilogy(EbN0_db,line,'m--');

is_equal_ecc = reshape(dm_7bits_transpose(4:7,:,:) == pre_hamming,N/4,4,length(EsN0));
% size(is_equal_ecc)
bit_error_num_ecc = squeeze(N - (sum(is_equal_ecc(:,1,:))+sum(is_equal_ecc(:,2,:))+sum(is_equal_ecc(:,3,:))+sum(is_equal_ecc(:,4,:))));
hold on
semilogy(EbN0_db,bit_error_num_ecc./N);
hold on
semilogy(EbN0_db,ber_ecc,'x--');

xlabel('Eb/N0(dB)');
ylabel('Bit Error Rate');
legend('BER(no ECC)','10^-4','BER(with ECC)','theoretical BER with ECC');
title('BER');

% effective coding gain
% code_gain = 4*k/n - k*log(2)./EbN0;

function y=Q(x)
y=erfc(x./sqrt(2))/2;
end

function y=QAM_Map_Symbol(input_qua, d, l)
s1 = [-3*d/2; 3*d/2];s2 = [-d/2; 3*d/2];s3 = [d/2; 3*d/2];s4 = [3*d/2; 3*d/2];s5 = [-3*d/2; d/2];
s6 = [-d/2; d/2];s7 = [d/2; d/2];s8 = [3*d/2; d/2];s9 = [-3*d/2; -d/2];s10 = [-d/2; -d/2];s11 = [d/2; -d/2];
s12 = [3*d/2; -d/2];s13 = [-3*d/2; -3*d/2];s14 = [-d/2; -3*d/2];s15 = [d/2; -3*d/2];s16 = [3*d/2; -3*d/2];

bit1 = [0,0,1,0].'; bit2 = [0,1,1,0].'; bit3 = [1,1,1,0].'; bit4 = [1,0,1,0].'; bit5 = [0,0,1,1].'; bit6 = [0,1,1,1].';bit7 = [1,1,1,1].'; 
bit8 = [1,0,1,1].'; bit9 = [0,0,0,1].'; bit10 = [0,1,0,1].'; bit11 = [1,1,0,1].'; bit12 = [1,0,0,1].'; bit13 = [0,0,0,0].'; bit14 = [0,1,0,0].'; 
bit15 = [1,1,0,0].'; bit16 = [1,0,0,0].'; 
s_i = [s1, s2, s3, s4, s5, s5, s6, s7, s8, s9, s10, s11, s12, s13, s14, s15, s16];

b1 = repmat(bit1, 1, l/4); b2 = repmat(bit2, 1, l/4); b3 = repmat(bit3, 1, l/4); b4 = repmat(bit4, 1 ,l/4); b5 = repmat(bit5, 1, l/4); 
b6 = repmat(bit6, 1, l/4); b7 = repmat(bit7, 1, l/4); b8 = repmat(bit8, 1, l/4); b9 = repmat(bit9, 1, l/4); b10 = repmat(bit10, 1, l/4); 
b11 = repmat(bit11, 1, l/4); b12 = repmat(bit12, 1, l/4); b13 = repmat(bit13, 1, l/4); b14 = repmat(bit14, 1, l/4); b15 = repmat(bit15, 1, l/4); 
b16 = repmat(bit16, 1, l/4); 

% xor the target and the signal, return matrices 
% which have the location of the corresponding mapping signal
symbol1_mat = xor(input_qua,b1);symbol2_mat = xor(input_qua,b2);symbol3_mat = xor(input_qua,b3);symbol4_mat = xor(input_qua,b4);
symbol5_mat = xor(input_qua,b5);symbol6_mat = xor(input_qua,b6);symbol7_mat = xor(input_qua,b7);symbol8_mat = xor(input_qua,b8);
symbol9_mat = xor(input_qua,b9);symbol10_mat = xor(input_qua,b10);symbol11_mat = xor(input_qua,b11);symbol12_mat = xor(input_qua,b12);
symbol13_mat = xor(input_qua,b13);symbol14_mat = xor(input_qua,b14);symbol15_mat = xor(input_qua,b15);symbol16_mat = xor(input_qua,b16);
% logic index, to identify the mapping signals' location
% those who in the symbol_mat that equals to 0 after xor 
% is the corresponding mapping
location_s1 = ~logical(sum(symbol1_mat));location_s2 = ~logical(sum(symbol2_mat));location_s3 = ~logical(sum(symbol3_mat));
location_s4 = ~logical(sum(symbol4_mat));location_s5 = ~logical(sum(symbol5_mat));location_s6 = ~logical(sum(symbol6_mat));
location_s7 = ~logical(sum(symbol7_mat));location_s8 = ~logical(sum(symbol8_mat));location_s9 = ~logical(sum(symbol9_mat));
location_s10 = ~logical(sum(symbol10_mat));location_s11 = ~logical(sum(symbol11_mat));location_s12 = ~logical(sum(symbol12_mat));
location_s13 = ~logical(sum(symbol13_mat));location_s14 = ~logical(sum(symbol14_mat));location_s15 = ~logical(sum(symbol15_mat));
location_s16 = ~logical(sum(symbol16_mat));
% make it two dimensional to later point product the symbol
%i.e. generate the mapped signal
% using the location index to get the mapped signal value
signal1 = location_s1.*s1;signal2 = location_s2.*s2;signal3 = location_s3.*s3;signal4 = location_s4.*s4;signal5 = location_s5.*s5;
signal6 = location_s6.*s6;signal7 = location_s7.*s7;signal8 = location_s8.*s8;signal9 = location_s9.*s9;signal10 = location_s10.*s10;
signal11 = location_s11.*s11;signal12 = location_s12.*s12;signal13 = location_s13.*s13;signal14 = location_s14.*s14;signal15 = location_s15.*s15;
signal16 = location_s16.*s16;
signals = signal1+signal2+signal3+signal4+signal5+signal6+signal7+signal8+signal9+signal10+signal11+signal12+signal13+signal14+signal15+signal16;
y = signals;
end


