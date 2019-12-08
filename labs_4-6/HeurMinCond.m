function [result] = HeurMinCond(A_inf, A_sup, iter_number)
  
m = size(A_inf,1); 
n = size(A_inf,2);   
 
if(nargin >= 3)
    NN = iter_number; 
else
    NN = 10;
end

Matr1 = ones(m,n); 
Matr2 = ones(m,n);   

MinCond = Inf; 
  
  
for k = 1:NN 

    EPM = randi([0,1],m,n); 
          
    for i = 1:m
        for j = 1:n
            if EPM(i,j) == 0 
                Matr1(i,j) = A_inf(i,j); 
                Matr2(i,j) = A_sup(i,j); 
            else 
                Matr1(i,j) = A_sup(i,j);
                Matr2(i,j) = A_inf(i,j);                 
            end 
        end
    end 

    c1 = cond(Matr1,2); 
    c2 = cond(Matr2,2); 
    if MinCond > c1 
        MinCond = c1; 
    end
    if MinCond > c2 
        MinCond = c2; 
    end
    
end

result = MinCond;

end