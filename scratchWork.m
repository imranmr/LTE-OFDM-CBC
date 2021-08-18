%this is for the encryption part
load('Proj1InputData.mat')

iv = ones(1,1024);
beforexor = iv;
key = toeplitz(mod(1:1,2),mod(1:1024,2));
encryptionVec = zeros(1,length(InputData));

for i=1:length(InputData)/1024
    afterxor = xor(beforexor,InputData(1024*(i-1)+1:1024*i));
    afterenc = xor(afterxor,key);
    encryptionVec(1024*(i-1)+1:1024*i) = afterenc;
    beforexor = afterenc;
end

%this is for the decryption part

iv2 = ones(1,1024);
input = encryptionVec;
key2 = toeplitz(mod(1:1,2),mod(1:1024,2));
decryptionVec = zeros(1,length(InputData));
len = length(input)/1024;

for j=1:length(input)/1024
    % get n and n-1 chunks
    current = input(1024*(len-j)+1:1024*(len-j+1));
    beforeDec = xor(key2,current);
    % n-1 chunk is the iv for the n chunk
    if j == length(input)/1024
        lastChunk = iv2;
    else
        lastChunk = input(1024*(len-j-1)+1:1024*(len-j));
    end
    plainText = xor(beforeDec,lastChunk);
    decryptionVec(1024*(len-j)+1:1024*(len-j+1)) = plainText;
end


checkifequal = isequal(decryptionVec, InputData);






array = ones(1, 3)
message = [message array]


