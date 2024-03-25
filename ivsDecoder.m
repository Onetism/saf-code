%clear
% clc
function spikeseq = ivsDecoder(filename, tnum)
close all

w=400;
h=250;

fid = fopen(filename, 'rb');
spikeseq = zeros(h, w, tnum, 'uint8');
fprintf('loading ivs data...\n');
para_frame = 1;

for t=1:tnum
    r=fread(fid,[(w+16)/8*para_frame h],'uint8');
    if size(r) == 0
        break
    end

    r_bin = dec2bin(r(1:w/8, :), 8); % 12500*8 char
    r_bin = r_bin.' - '0'; % 8*12500 double
    r_bin = flip(r_bin, 1); % 8*12500
    r_bin = reshape(r_bin, [], 1); %100000*1
    tmp_spkseq = reshape(r_bin, w, h, para_frame); % 400*250
    spikeseq(:, :, (t-1)*para_frame+1:para_frame*t) = permute(tmp_spkseq,[2,1,3]); 
 
end
if t~=tnum
    t = t -1;
end
fprintf('%d frames have loaded.\n', t);
spikeseq = spikeseq(:,:,1:t);
fclose(fid);
end