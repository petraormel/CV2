function [R, t] = icp(source, target)
    R = eye(3);

    t = 0;

    mins = zeros(size(source,2),1);
    size(target,2)
    minj = 0;
    source(:,1)
    target(:,1)
    
    for i=1:size(source,2)
        minsource = inf;
        for j=1:size(target,2)
            if norm(source(:,i)-target(:,j)) < minsource
                
                minsource = norm(source(:,i)-target(:,j));
                minj = j;
               
            end
        end
        mins(i) = minj;    
    end
    mins

end