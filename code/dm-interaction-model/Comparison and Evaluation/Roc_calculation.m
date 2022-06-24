% the Output is ROC, % AUC and ACC
% AUC : Area Under the Curve
% ACC : Accuracy [Roc_f AUC ACC]
function [Roc_f ,AUC, ACC] = Roc_calculation(f,Ground_truth)
save_TP_f = [];
save_FP_f = [];
condition_positive = sum(Ground_truth);%(norm(123-249)+1)+ (norm(301-360)+1) + (norm(670-770)+1);
condition_negative = length(Ground_truth) - condition_positive;
deltaStep = 0.03;

for i = 0:deltaStep:1
TP_f = 0;
FP_f = 0;
for j = 1: length(f)
if i == 0
if f(j)>= i
if Ground_truth(j)== 1
TP_f = TP_f + 1;
else
FP_f = FP_f +1;
end
end
else
if f(j)> i
if Ground_truth(j)== 1
TP_f = TP_f + 1;
else
FP_f = FP_f +1;
end
end 
end 
end
save_TP_f = [save_TP_f TP_f];
save_FP_f = [save_FP_f FP_f];
end
true_positive_f = save_TP_f/condition_positive;
false_positive_f = save_FP_f/condition_negative;
Roc_f = [false_positive_f;true_positive_f];
AUC = -trapz(false_positive_f,true_positive_f);
TN = condition_negative*(1-false_positive_f);
ACC = max(save_TP_f + TN)/(condition_positive + condition_negative);
