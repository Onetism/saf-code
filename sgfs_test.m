close all;
clear;
clc;

%%
O_MFDCT = [1, 1, -1, -1; 1, 1, -1, -1; -1, -1, 1, 1; -1, -1, 1, 1];
data_config = readtable('datasets/config.xlsx','VariableNamingRule','preserve');
disp(data_config(:,1:2));

for d = 1 : size(data_config,1)
data_index = d;
% data_index = input('Enter the datasets ID for test : ');
disp([data_config(data_index,2).(1){1} ' have been selected  for testing.']);

%%
maxframe = data_config(data_index,3).(1);
filepath = data_config(data_index,6).(1);
spikeseq = ivsDecoder(filepath{1}, maxframe);

%%
dt = data_config(data_index,5).(1);
sd_score = zeros(floor(size(spikeseq, 3)/dt), 1);
er_score = zeros(size(sd_score));
gd_score = zeros(size(sd_score));
mfdct_score = zeros(size(sd_score));
sgfs_score = zeros(size(sd_score));
for i = 1:dt:size(spikeseq, 3)-dt+1
    index = floor(i/dt)+1;
    power = sum(spikeseq(:, :, i:i+dt-1), 3);
    intensity = mean(power,'all');

    sd_score(index, 1) = sum(sum((power - intensity).^2))/(intensity^2) ;
    er_score(index, 1) = sum(sum( (power).^2 ));
    

    image = power;
    [gradient_magnitude, gradient_direction] = imgradient(image);
    gd_score(index, 1) = sum(gradient_magnitude, "all");
    
    image = image(1:248,:);
    blockRows = 4;
    blockCols = 4;

    blocks = mat2cell(image, repmat(blockRows, 1, size(image, 1) / blockRows), repmat(blockCols, 1, size(image, 2) / blockCols));
    for n = 1 : numel(blocks)
            m = blocks{n};
            mfdct_score(index, 1) = mfdct_score(index, 1) + sum(conv2(m, O_MFDCT).^2, 'all');
    end

    if (index < size(sd_score,1)/4)
        diff = 0;
    else 
        powerL = sum(spikeseq(:, :, 1:floor((i+dt)/2)-1), 3);
        powerH = sum(spikeseq(:, :, floor((i+dt)/2)-1:i+dt-1), 3);
    
        intensityL = mean(powerL,'all');
        intensityH = mean(powerH,'all');
        
        powerL_sd_score = sum(sum((powerL - intensityL).^2))/(intensityL^2);
        powerH_sd_score = sum(sum((powerH - intensityH).^2))/(intensityH^2);
        
        diff = powerH_sd_score - powerL_sd_score;
    end
    sgfs_score(index, 1) = diff;

end

%%
stop_index = find(sgfs_score == max(sgfs_score(floor(size(sd_score,1)/4):end)));
left_index = 1;
right_index = stop_index * dt;
left_index = floor(right_index - 0.618*(right_index - left_index));
poseA = [floor((left_index+right_index)/2)];
while(right_index - left_index > 10)
        L = sum(spikeseq(:, :, left_index:left_index+0.618*(right_index - left_index)), 3);
        R = sum(spikeseq(:, :, right_index-0.618*(right_index - left_index) : right_index), 3);
        
        meanL = mean(L,'all');
        meanR = mean(R,'all');

        loosA = sum(sum((L - meanL).^2))/(meanL^2);
        loosB = sum(sum((R - meanR).^2))/(meanR^2);

        if(loosA > loosB)
            right_index = left_index+0.618*(right_index - left_index);
        else
            left_index = right_index-0.618*(right_index - left_index);
        end
        poseA = [poseA, floor((left_index+right_index)/2)];
end

left_index = 1;
right_index = size(spikeseq,3);
poseB = [floor((left_index+right_index)/2)];
while(right_index - left_index > 10)
        L = sum(spikeseq(:, :, left_index:left_index+0.618*(right_index - left_index)), 3);
        R = sum(spikeseq(:, :, right_index-0.618*(right_index - left_index) : right_index), 3);
        
        meanL = mean(L,'all');
        meanR = mean(R,'all');

        loosA = sum(sum((L - meanL).^2))/(meanL^2);
        loosB = sum(sum((R - meanR).^2))/(meanR^2);

        if(loosA > loosB)
            right_index = left_index+0.618*(right_index - left_index);
        else
            left_index = right_index-0.618*(right_index - left_index);
        end
        poseB = [poseB, floor((left_index+right_index)/2)];
end

min_value = min(er_score);
max_value = max(er_score);
er_score = (er_score - min_value) / (max_value - min_value);
subplot(3,2,1);
plot(er_score);
title('Event Rate');

min_value = min(gd_score);
max_value = max(gd_score);
gd_score = (gd_score - min_value) / (max_value - min_value);
subplot(3,2,2);
plot(gd_score);
title('Gradient');

min_value = min(mfdct_score);
max_value = max(mfdct_score);
mfdct_score = (mfdct_score - min_value) / (max_value - min_value);
subplot(3,2,3);
plot(mfdct_score);
title('MF-DCT');

min_value = min(sd_score);
max_value = max(sd_score);

sd_score = (sd_score - min_value) / (max_value - min_value);
subplot(3,2,4);
plot(sd_score);
focal_positon_str = strjoin(cellstr(num2str(dt*find(sd_score==max(sd_score)).')), ', ');
title(['Spike Dispersion : ' focal_positon_str]);

subplot(3,2,5);
plot(poseA);
focal_positon_str = strjoin(cellstr(num2str(poseA(end).')), ', ');
title(['SGFS:' focal_positon_str]);

subplot(3,2,6);
plot(sgfs_score);
title('Real Time SGFS');

filefolder= fileparts(filepath{1});
saveas(gcf, [filefolder '/score_list.png']);
save([filefolder '/score_list.mat'], 'er_score', 'gd_score', 'mfdct_score', 'sd_score','sgfs_score','poseA','poseB');

scen_name = data_config(data_index,2).(1){1};
saveas(gcf, ['./datasets/' scen_name '.png']);
end