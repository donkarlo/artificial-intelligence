function bag = bagDataRead(numbOfBag, dirData, numData)
label = num2str(numData);
codification{1,1} = [label '.bag'];
if numbOfBag > 4 || numbOfBag < 1
    display('Error: Bags IDs must go from 1 to 4')
else
    stringCod = codification{numbOfBag,1};
    cd(dirData)
    bag = rosbag(stringCod);
end
end