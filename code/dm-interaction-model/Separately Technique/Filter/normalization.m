function[distFinal] = normalization(currState, averageMat)
inThreshProb = 0.8;                                                         % If the probability of being in a superstate is 'inThreshProb', it is considered to be imposible to be outside the model
outThreshProb = 1 - inThreshProb;
nodesNumber = size(averageMat,1);                                           % Total of neurons
distances = pdist2(currState,averageMat);                                  % Calculation of norm 1 distance from current positions to each neuron
distFinal = distances./radius;

superStateProb = 1 - distFinal;                      % value 1 means mean and 0 is outside boundary

superStateProb = superStateProb .* (distFinal<1);                           % Sum elements
prob = max(superStateProb);                                                 % Calculation of maximum probability

superStateProb(1,nodesNumber+1) = max([1-(prob+outThreshProb), 0]);             % Probability of empty neuron is 1-probability of most likely neuron
superStateProb = superStateProb/sum(superStateProb);                        % Normalization
end
