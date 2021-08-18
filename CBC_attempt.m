clc;
clear all;
close all;
load("Proj1InputData.mat")

%Generates random 64 bits for both IV and Key
IV = randi([0,1], 1, 64);
k = randi([0,1], 1, 64);

%Proj1 Input Data
Message = InputData;

%Function calling CBC_Encrypt
CT = CBC_Encrypt(Message, IV, k);

%Function calling CBC_Decrypt
PT = CBC_Decrypt(CT, IV, k);

%Checking the output of the decryption with the original data
isequal(PT(1:length(Message)), Message)
eq(PT(1:length(Message)), Message);


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
