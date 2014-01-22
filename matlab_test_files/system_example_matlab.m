syst_fake=tf([1],[1 2 3])
syst_fake_dis=c2d(syst_fake,0.01)
[output,t]=step(syst_fake_dis)
plot(output)

out_len = length(output)
input=1:out_len;
input(:)=1;
[num,den]=stmcb_test(output,input,0,2);
sys_model=tf(num,den,0.01)
step(sys_model)
hold on
step(syst_fake)


% [b,a] = butter(6,0.2);              % Butterworth filter design
% h = filter(b,a,[1 zeros(1,100)]);   % Filter data using above filter
% [bb,aa] = stmcb_test(h,4,4)
