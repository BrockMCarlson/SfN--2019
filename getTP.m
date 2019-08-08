function TP = getTP(pEvC,pEvT)
count = 0;

for i = 1:size(pEvC,2)
	if ~any(pEvC{i} == 96) % This is not necessary on the evp trails
        % skip if trial was aborted and animal was not rewarded (event code 96)
        continue
    end   
    count= count+1;
    stimon   =  pEvC{i} == 23;
    stimoff  =  pEvC{i} == 24;
    
    for j = 1:size(pEvC{i},1)
        if pEvC{i}(stimon) 
            TP{count} = [pEvT{i}(stimon); pEvT{i}(stimoff)];

        end
    end
end
<<<<<<< Updated upstream
end
=======
end
>>>>>>> Stashed changes
