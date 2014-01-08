%% testing prony function versus python function
%   % Example:
%   %   Fit an IIR model to an impulse response of a lowpass filter.
%
[b,a] = butter(4,0.2);
% # b = 0.0048    0.0193    0.0289    0.0193    0.0048
% # a = 1.0000   -2.3695    2.3140   -1.0547    0.1874
b
a
impulseResp = impz(b,a);                % obtain impulse response
impulseResp
% impulseResp =
% 
%     0.0048
%     0.0307
%     0.0906
%     0.1679
%     0.2246
%     0.2335
%     0.1935
%     0.1238
%     0.0496
%    -0.0085
%    -0.0407
%    -0.0476
%    -0.0369
%    -0.0186
%    -0.0013
%     0.0100
%     0.0140
%     0.0121
%     0.0071
%     0.0017
%    -0.0022
%    -0.0040
%    -0.0039
%    -0.0026
%    -0.0010
%     0.0004
%     0.0011
%     0.0012
%     0.0009
%     0.0004
%     0.0000
%    -0.0003
%    -0.0004
%    -0.0003
%    -0.0002
%    -0.0000
%     0.0001
%     0.0001
%     0.0001
%     0.0001
%     0.0000
%    -0.0000
%    -0.0000
% python [0.00482434335771622,0.0307287177680857,0.0905946819548826,0.167944821844737,0.224641271344028]
denOrder=4; numOrder=4;                 % system function of order 4
[Num,Den]=prony_test(impulseResp,numOrder,denOrder);
% Num =0.0048    0.0193    0.0289    0.0193    0.0048
% Den =1.0000   -2.3695    2.3140   -1.0547    0.1874
display('Num')
Num
display('Den')
Den
% subplot(211);                           % impulse response and input
% stem(impz(Num,Den,length(impulseResp)));   
% title('Impulse Response with Prony Design');
% subplot(212);
% stem(impulseResp); title('Input Impulse Response');