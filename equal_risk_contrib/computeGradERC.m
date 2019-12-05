function gval = computeGradERC (x)

global Q
  
  n = size(Q,1) ;
  k = ((n*n)-n)/2; 

  if(size(x,1)==1)
     x = x';
  end
  
  % Insert your gradiant computations here
  % You can use finite differences to check the gradient
  
  %calculating y 
  y = x.*Q*x;
  
  %calculating all the possible partial derivative y index =i, x index = j
  derv_y = zeros(size(Q));
  
  for i = 1:n
      for j = 1:n
          if j==i
              derv_y(i,j) = Q(i,j)*x(i) + dot(Q(i,:),x);
          else
              derv_y(i,j) = Q(i,j) * x(i);
          end                   
      end
  end
  
  %calculating yi - yj
  y_diff = zeros(1,k);
  
  count=1;
  for i = 1:n-1
      for j = i+1:n
          y_diff(count) = y(i)-y(j);  
          count = count +1;
      end
  end
  
  %calculating partial derivative difference 
  derv_y_diff = zeros(n,k);
  
  count=1;
  for index = 1:n
      for i = 1:n-1
          for j = i+1:n
              derv_y_diff(index,count) = derv_y(i,index) - derv_y(j,index);
              count = count +1;
          end
      end
      count = 1;
  end

  
  %calculating the gradient for obj function
  gval = zeros(n,1); 
  
  for i = 1:n
      gval(i) = 4*y_diff*derv_y_diff(i,:)';
  end

end
