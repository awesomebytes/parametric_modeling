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
% python [0.00482434335771622,0.0307287177680857,0.0905946819548826,0.167944821844737,0.224641271344028]
denOrder=4; numOrder=4;                 % system function of order 4
[Num,Den]=prony(impulseResp,numOrder,denOrder);
% Num =0.0048    0.0193    0.0289    0.0193    0.0048
% Den =1.0000   -2.3695    2.3140   -1.0547    0.1874
Num
Den
subplot(211);                           % impulse response and input
stem(impz(Num,Den,length(impulseResp)));   
title('Impulse Response with Prony Design');
subplot(212);
stem(impulseResp); title('Input Impulse Response');