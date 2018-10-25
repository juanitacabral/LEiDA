%Entropy
load DEFINITIVE2711.mat P
probabilities= squeeze(P(9,:,1,1:10));

probabilitiesState8= squeeze(P(9,:,1,8));

%probabilities = a matrix of 86*10 containing probabilities of the ten
%states for each subject. s=1:51 is patients. s=52:86 is controls.

%our main result is that the probability of state 8 significantly differs
%between the groups (higher probabilities in controls). 


%entropy of the whole system for k=10

%E = -sum(P.*log2(P))

for s=1:86
 
H(s)=-sum(probabilities(s,:).*log(probabilities(s,:)));
 
end


%for k=10, corrected for values with 0 probability
for s=1:86
H(s)=0
for c=1:10
if probabilities(s,c)>0
H(s)=H(s)-sum(probabilities(s,c).*log2(probabilities(s,c)));
end
end
end





%note that some probabilities are 0 (for individual subjects), thus
%some entropy values are NaN. In neutral mood, replace 0 of subject 32 with lowest value 

%entropy of state 8 (for k=10) corrected for 0 probabilities


for s=1:86
H8(s)=0
if probabilitiesState8(s)>0
H8(s)=H8(s)-(probabilitiesState8(s).*log2(probabilitiesState8(s)));
end
end

p_H = ranksum(H(1:51), H(52:86))

p_H8 = ranksum(H8(1:51), H8(52:86))

p_P8 = ranksum(probabilitiesState8(1:51), probabilitiesState8(52:86))

%results when 0 into account

% p_H =
% 
%     0.6288
% 
% 
% p_H8 =
% 
%     0.0057


%results:

% p_H =
% 
%     0.4269
% 
% 
% p_H8 =
% 
%     0.0080
% 
% 
% p_P8 =
% 
%     0.0057

% mean(H8(1:51))

%     0.1869
% 
% mean(H8(52:86))
% 
% ans =
% 
%     0.2530

%Entropy SAD
load DEFINITIVE2711.mat P
probabilitiessad= squeeze(P(9,:,2,1:10));

probabilitiesState8sad= squeeze(P(9,:,2,8));

%probabilities = a matrix of 86*10 containing probabilities of the ten
%states for each subject. s=1:51 is patients. s=52:86 is controls.

%our main result is that the probability of state 8 significantly differs
%between the groups (higher probabilities in controls). 


%entropy of the whole system for k=10


for s=1:86

Hsad(s)=-sum(probabilitiessad(s,:).*log(probabilitiessad(s,:)));

end

for s=1:86
Hsad(s)=0
for c=1:10
if probabilitiessad(s,c)>0
Hsad(s)=Hsad(s)-sum(probabilitiessad(s,c).*log2(probabilitiessad(s,c)));
end
end
end


%note that some probabilities are 0 (for individual subjects), thus
%some entropy values are NaN

%entropy of state 8 (for k=10)


for s=1:86

H8sad(s)=-(probabilitiesState8sad(s).*log(probabilitiesState8sad(s)));

end

for s=1:86
H8sad(s)=0
if probabilitiesState8sad(s)>0
H8sad(s)=H8sad(s)-(probabilitiesState8sad(s).*log2(probabilitiesState8sad(s)));
end
end

p_Hsad = ranksum(Hsad(1:51), Hsad(52:86))

p_H8sad = ranksum(H8sad(1:51), H8sad(52:86))

p_P8sad = ranksum(probabilitiesState8sad(1:51), probabilitiesState8sad(52:86))

%results:

% p_H =
% 
%     0.4269
% 
% 
% p_H8 =
% 
%     0.0080
% 
% 
% p_P8 =
% 
%     0.0057
