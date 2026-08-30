clear, clc

%% Test of chargestate_analysis Function
clear, clc
% Inputs
qVec = [-3;-2;-1; 0; 1; 2; 3];
E0Vec = [3.3;3.2;2.4;2.0;2.5;3.0;2.8];
T = 300; 
numKbt = 3; 

% Function Call
[CTLVec,EVec,qPair,CTL,CTLMinVec,EMinVec,qMinPair,minCTL,qMinVec,E0MinVec] = chargestate_analysis(qVec,E0Vec,T,numKbt)

%% Test 2
clear, clc
% Inputs
qVec = [-1; 0; 1; 2; 3];
E0Vec = [3.3;3.2;2.5;3.0;2.8];
T = 300; 
numKbt = 3; 

% Function Call
[CTLVec,EVec,qPair,CTL,CTLMinVec,EMinVec,qMinPair,minCTL,qMinVec,E0MinVec] = chargestate_analysis(qVec,E0Vec,T,numKbt)

%% Test 3: Still can't figure out what is wrong here ?????
clear, clc
% Inputs
qVec = [1; 2; 3];
E0Vec = [3.3;3.2;2.4];
T = 300; 
numKbt = 3; 

% Function Call
[CTLVec,EVec,qPair,CTL,CTLMinVec,EMinVec,qMinPair,minCTL,qMinVec,E0MinVec] = chargestate_analysis(qVec,E0Vec,T,numKbt)