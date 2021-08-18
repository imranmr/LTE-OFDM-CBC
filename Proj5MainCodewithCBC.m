clc;
clear all;
close all;
load('Proj1PNSeq.mat');
load('Proj1TestData.mat');
load('Proj1TestEncryptedData.mat');
load('Proj1ModSymb.mat');
load('Proj1TransSymbStream.mat');
load('Proj1UpCnvrt2IFData.mat');
load('Proj1InputData.mat')

%Generates random 64 bits for both IV and Key
IV = randi([0,1], 1, 64);
%key is generated from PN seq


%Proj1 Input Data
Message = InputData;

x1 = [1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
n1 = length(x1);
len1 = 2^n1-1;
p1 = zeros(1,len1);
p1(1,1) = 0; %initialize to 0 due to x1 ending with 0
z1 = x1;
for y1 = 2:len1
    x1 = z1;
    for i = 1:n1 
        if (i ==2) %for x
            z1(1,i) = xor(x1(1,1),x1(1,19));
        elseif (i==3) %xor bit2 and bit5 and send into bit3 for x2
            z1(1,i) = xor (x1(1,2), x1(1,19));
        elseif (i==6) %for x5
            z1(1,i) = xor (x1(1,5),x1(1,19));
        elseif (i==1)
            z1(1,1) = x1(1,19);
        else %else just pass the value over
            z1(1,i) = x1(1,i-1);
        end
    end
    p1(1,y1) = z1(1,19);
end
%disp("M-sequence is p1")
%fprintf("%i ", p1);

%save PNSequence p1;

C = isequal(p1,PNSeq); %compares given PNSeq to generated PNSeq (1 == true)

%---------------------end of PNSeq-------------------------%
%----------------------------------------------------------%
%-----------------start of Encryption----------------------%

S = repmat(p1, 1, 39); %scale up PNSeq to size of TestData
modSymb = length(InputData)-length(S);
dif = p1(1:modSymb);
PNtotal = [S dif]; %Complete PNSeq
%TestEncryptedD = xor(PNtotal, InputData);
shift = 200000; %extracting key from PN sequence
key = PNtotal(shift:shift+63); 
TestEncryptedD = CBC_Encrypt(Message, IV, key);

save EncryptedData TestEncryptedD; %Encrypted bit stream (serial)

%---------------------end of Encryption-------------------------%
%----------------------------------------------------------%
%-----------------start of Mapping----------------------%

data = TestEncryptedD;

bw = 10000000; %10MHz, Starting bandwidth

qpsk_mapping = zeros(1,length(data));

%qpsk mapping 
for j=1:length(data)

    if data(j) == 0
        qpsk_mapping(j) = 1/sqrt(2);
    elseif data(j) == 1
        qpsk_mapping(j) = -1*(1/sqrt(2));
    end
    
end

%reformatting the structure of the QPSK mapping to match the following
    %1024M rows of data with 2 columns. First column is I, Second column is
    %Q

qpsk_mapping = reshape(qpsk_mapping, 2, length(qpsk_mapping)/2);
qpsk_mapping = qpsk_mapping';


%%%%CODE CHECKING%%%%%%%%%
modSymb = zeros(1,10240000);
for j=1:length(qpsk_mapping) 
    modSymb(j) = complex(qpsk_mapping(j,1),qpsk_mapping(j,2));
end
F = isequal (modSymb, ModSymb);
%%%%END OF CODE CHECK%%%%

save ModSymbTransmitter modSymb; %Bit to Modulation symbol mapping (serial)


% cyclic prefix 
% 
cp = zeros(1094,10000);

% need to fix code here for make it run faster
% storing temp to get the OFDMSymb is taking a long time
temp = zeros(1024,10000);
for i=1:length(cp)
    inverse_fft = ifft(modSymb(1024*(i-1)+1:1024*i));% get first 1024 and then next 1024 10240000
    temp2 = inverse_fft(955:1024); % last 70 
    temp3= [temp2 inverse_fft]; % this 1 symbol and cp 
    
    for j = 1:1094
        cp(j,i) = temp3(j);
    end
    
    for k = 1:1024
        temp(k,i) = inverse_fft(k);
    end
        
        
end

x =reshape(cp,1,[]);
inv = reshape(temp,1,[]);
check_if_equal = isequal (inv,TransSymbStream);

save OFDMsymb temp; % OFDM symbol generation (parallel) (post-IFFT)
save TransSymb x; % CP insertion (serial)

%%%%%%%%%%%%%%%%%%%%removed the upconversion%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%      END OF TRANSMITTER SIDE     %%%%%%%%%%%
amplitude = normrnd(0,0.25,[1,length(x)]);
phase = normrnd(0,sqrt(deg2rad(20)),[1,length(x)]);

complexForm = zeros(1,length(x));

for i=1:length(x)
    complexForm(i) = complex(amplitude(i)*cos(phase(i)),amplitude(i)*sin(phase(i)));
    
end

%do ifft here
complexForm = ifft(complexForm);

%add everything
noiseAddition = x + complexForm;

save NoiseAddition noiseAddition; % Additive channel noise (serial) 

%%%%%%%%%%%%      START OF RECEIVER SIDE     %%%%%%%%%%%

startReceiver = noiseAddition;
CPremove = zeros(1,10240000);
for i=1:10000
    chunk = startReceiver((1094*(i-1))+1:(1094*(i-1))+1094);
    rmCP = chunk(71:1094);
    for j=1:length(rmCP)
        CPremove((1024*(i-1))+j) = rmCP(j);
    end
end

save NOCP CPremove; % CP removal (serial) 

s2p = reshape(CPremove,1,[]);

%save OfdmSymStream s2p;

start = s2p;

noCP = zeros(1,length(start));

for i=1:length(start)/1024
    ft = fft(start(1024*(i-1)+1:1024*i));
    for j=1:1024
        noCP((1024*(i-1))+j) = ft(j);
    end
end

array = zeros(length(start), 2);
length(array);
for j=1:length(array)
    array(j,1) = real(noCP(j));
    array(j,2) = imag(noCP(j));
end

parallel_to_serial= zeros(1,length(array));
for j=1:length(array) 
    parallel_to_serial(j) = complex(array(j,1),array(j,2));
end

save ModSymbReceiver parallel_to_serial; % Modulation symbol recovery (serial) (post-FFT and P/S) 

%QPSK mapping========================================================

%for loop to iterate through parallel_to_serial
%the variable flag is to make sure that we put the bits in the correct
%place in new without using another loop 

encrypted_bits = zeros(1,length(parallel_to_serial)*2);
flag=0;
for i=1:length(array)
    if real(parallel_to_serial(i))>0 && imag(parallel_to_serial(i))>0 %top right
        encrypted_bits(i+flag:i+1+flag)=[0,0];
        
    elseif real(parallel_to_serial(i))>0 && imag(parallel_to_serial(i))<0 %bottom right
        encrypted_bits(i+flag:i+1+flag)=[0,1];
    
    elseif real(parallel_to_serial(i))<0 && imag(parallel_to_serial(i))>0 %top left
        encrypted_bits(i+flag:i+1+flag)=[1,0];
       
    elseif real(parallel_to_serial(i))<0 && imag(parallel_to_serial(i))<0 %bottom left
        encrypted_bits(i+flag:i+1+flag)=[1,1];
       
    elseif real(parallel_to_serial(i))==0 && imag(parallel_to_serial(i))<0
        %random bottom 2 quadrants
        temp = randi(2);
        if temp == 1
            encrypted_bits(i+flag:i+1+flag)=[1,1];
        elseif temp == 2
            encrypted_bits(i+flag:i+1+flag)=[0,1];
        end
    elseif real(parallel_to_serial(i))==0 && imag(parallel_to_serial(i))>0
        %random top 2 quadrants
        temp = randi(2);
        if temp == 1
            encrypted_bits(i+flag:i+1+flag)=[0,0];
        elseif temp == 2
            encrypted_bits(i+flag:i+1+flag)=[1,0];
        end
    elseif imag(parallel_to_serial(i))==0 && real(parallel_to_serial(i))>0
        %random right 2 quadrants
        temp = randi(2);
        if temp == 1
            encrypted_bits(i+flag:i+1+flag)=[0,0];
        elseif temp == 2
            encrypted_bits(i+flag:i+1+flag)=[0,1];
        end
    elseif imag(parallel_to_serial(i))==0 && real(parallel_to_serial(i))<0
        %random left 2 quadrants
        temp = randi(2);
        if temp == 1
            encrypted_bits(i+flag:i+1+flag)=[1,0];
        elseif temp == 2
            encrypted_bits(i+flag:i+1+flag)=[1,1];
        end
    end
    flag=flag+1;
end

save EncryptedBits encrypted_bits; %Modulation symbol to encrypted bit (serial)

%-------------------------------------------------------------------------------------------

% decryption 

%generate PN sequence 
x1 = [1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
n1 = length(x1);
len1 = 2^n1-1;
p1 = zeros(1,len1);
p1(1,1) = 0; %initialize to 0 due to x1 ending with 0
z1 = x1;
for y1 = 2:len1
    x1 = z1;
    for i = 1:n1 
        if (i ==2) %for x
            z1(1,i) = xor(x1(1,1),x1(1,19));
        elseif (i==3) %xor bit2 and bit5 and send into bit3 for x2
            z1(1,i) = xor (x1(1,2), x1(1,19));
        elseif (i==6) %for x5
            z1(1,i) = xor (x1(1,5),x1(1,19));
        elseif (i==1)
            z1(1,1) = x1(1,19);
        else %else just pass the value over
            z1(1,i) = x1(1,i-1);
        end
    end
    p1(1,y1) = z1(1,19);
end

%scale up PNSeq
S = repmat(p1, 1, 39);
modSymb = length(encrypted_bits)-length(S);
dif = p1(1:modSymb);

%The Complete PNSeq is PNtotal 
PNtotal = [S dif];
shift = 200000;
key = PNtotal(shift:shift+63);
%get the decrypted bits by performing xor on PNtotal and encrypted_bits 
%decrypted = xor(PNtotal,encrypted_bits);
decrypted = CBC_Decrypt(encrypted_bits, IV, key);

save De-encrypted decrypted; %De-encrypted bit stream (serial)

%----------------- End of Part 4 -----------------------------%

%%%%%%%%%%%%%%%%   CHECK FINAL OUTPUT   %%%%%%%%%%%%%%%5

finalchecker = isequal(InputData,decrypted)



function ciphertext = CBC_Encrypt(message, IV, Key)

    %condition: Mesage%64 != 0
    %If the data is not a multiple of 64, make it a multiple of 64
    if(mod(length(message), 64) ~= 0)
        array = zeros(1, 64 - mod(length(message), 64));
        message = [message array];
    end
    
    %Holding the initial value of IV for the first iteration
    tempIV = IV;
    Kfunc_output = zeros(1, 64);
    
    %Defining the structure of the output
    ciphertext = zeros(1,length(message));
    
    %Calculating the total number of loops
    loops = length(message)/64;
    
    %The for loop serves to xor the original message with IV, then
    %immediately xor with Key. IV is updated every loop as the output of
    %the Key function
    for i=1:loops
        %This variable (t) allows access to each 64 offset of the message
        t = (i - 1) * 64;
        for j=1:64
            ciphertext(1,t + j) = xor(message(1,t + j), tempIV(1,j));
            ciphertext(1, t+j) = xor(Key(1, j), ciphertext(1, t + j));
            Kfunc_output(1, j) = ciphertext(1, t + j); 
        end
        
        tempIV = Kfunc_output;
        
    end
   
end

function plaintext = CBC_Decrypt(ciphertext, IV, Key)
    
    temp = zeros(1,length(ciphertext));
    
    
    plaintext = zeros(1,length(ciphertext));

    loops = length(ciphertext)/64;
    
    %The for loop will run for all the blocks (64 bit) in the ciphertext
    %except for the first block
    for i=loops:-1:2
        %This variable (t) allows access to each 64 offset of the
        %ciphertext
        t = (i - 1) * 64;
        %IV_OS serves as the offset to obtain the temp IV from the
        %ciphertext
        IV_OS = (i - 2) * 64;
        
        %TempIV will take 64 bits of cipher text starting from the second 
        %to last least significant 64 bits and receeding to the second most
        %significant bit. I hope this part made sense... Just look at the
        %code if you are confused.
        tempIV = ciphertext(1, IV_OS + 1:IV_OS + 64);
        
        for k=1:64
            temp(1, t + k) = xor(Key(1,k), ciphertext(1, t + k));
            plaintext(1, t + k) = xor(temp(1, t + k), tempIV(1, k));
        end

    end
    
    % First block is taken care of in this last for loop
    for i=1:64
        temp(1,i) = xor(Key(1,i),ciphertext(1,i));
        plaintext(1,i) = xor(temp(1,i), IV(1,i));
    end

end


