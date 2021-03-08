%load the data for the graph
load('graph502.mat')
chosenGraph= graph502;
name = 'temp502';

%find the dimensions of the graph
a=length(chosenGraph);
[height, width, colour] = size(chosenGraph);

%loop to find all points of the graph
i=2;
plots = [0,16];

for n=1:a
    graph = [];
    
    %finds all red pixels in the column
    graph(:,n)=chosenGraph(:,n,1)>160 & chosenGraph(:,n,2)<150 & chosenGraph(:,n,3)<150;
    
    %find 
    points=find(graph(:,n),1);
    
    %Checks whether the column is empty
    noPoints = isempty(points);
    
    
    %if pixels are found
    if noPoints == false
        %finds mean posisiotn of all pixels in the column
        xValues = mean(points);
        %scales according the image
        x = (n)*2000/width;
        y = (height - xValues)*2000/height;
        %plots(n,:)=[n,xValues];
        %saves all values into an appropriately sized matrix
        plots(i,:) = [x,y];
        i=i+1;
    else
        
    end
end

%save the data to the correct variable names
tempdata = plots(:,2);
timedata = plots(:,1);

%scale the data
tempdata(i-4:i-2) = tempdata(i-1);


[timedata, index] = unique(timedata);
tempdata = tempdata(index);

save(name, 'timedata', 'tempdata')

%plot to ensure correct graph
plot(timedata,tempdata)
