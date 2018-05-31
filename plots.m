data = load('measurements.txt','-ascii');
i=find(data(:,4)>0,1);

figure
stairs((data(i:end,5)-data(i,5))/1e3, data(i:end,1),'r')
title('Angle')

figure
stairs((data(i:end,5)-data(i,5))/1e3, data(i:end,3),'r')
title('Controller Input')
