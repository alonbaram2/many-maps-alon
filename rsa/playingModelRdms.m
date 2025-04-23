scanFile1 = strtrim(ls(['/vols/Scratch/mgarvert/ManyMaps/scan_1.1/datafiles/Subj_1/*_session_1/data_1_1.mat']));
scanFile2 = strtrim(ls(['/vols/Scratch/mgarvert/ManyMaps/scan_1.1/datafiles/Subj_25/*_session_2/data_25_2.mat']));

load(scanFile1)
data1=data;
load(scanFile2)
data2=data;
clear data

clear ix ix1 ix2 order
data=data1;
for i = 1:17
    ix1(i,data.stimuli(2,:) == data.stimuli(1,i)) = 1;
    ix2(i,data.stimuli(1,:) == data.stimuli(2,i)) = 1;
    order(i) = find(data.stimuli(2,:) == data.stimuli(1,i));
    order2(i) = find(data.stimuli(1,:) == data.stimuli(2,i));
end
ix1 = 1-ix1;   ix2 = 1-ix2;

identity1 = [1-eye(17) ix1; ix2 1-eye(17)] ;

clear ix ix1 ix2 order
data=data2;
for i = 1:17
    ix1(i,data.stimuli(2,:) == data.stimuli(1,i)) = 1;
    ix2(i,data.stimuli(1,:) == data.stimuli(2,i)) = 1;
    order(i) = find(data.stimuli(2,:) == data.stimuli(1,i));
    order2(i) = find(data.stimuli(1,:) == data.stimuli(2,i));
end
ix1 = 1-ix1;   ix2 = 1-ix2;

identity2 = [1-eye(17) ix1; ix2 1-eye(17)] ;


