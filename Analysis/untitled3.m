cd D:\test\processed\catgt_test_8_g0

frames_start = readtable('test_8_g0_tcat.nidq.xia_7_0.txt')
size(frames_start)
A  = table2array(frames_start);
time = linspace(1,10,size(frames_start,1));

plot(time,A(:,1))

hist(A(:,1))

difference = table2array(diff(frames_start));

hist(difference(:,1))
