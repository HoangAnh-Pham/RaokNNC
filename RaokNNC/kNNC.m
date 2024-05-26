% Apply the k-nearest neighbor comparison
function eval = kNNC(p,x_new,X,Fitness,Const,en,k)  
  Xc=X; Fc=Fitness; Cc=Const;
  
  % Eliminate the p-solution in the population
  % Xc(p,:)=[]; Fc(p)=[]; Cc(p)=[];
  
  id = knnchk(x_new,Xc,k);
  for i=1:length(id)
    ecom(i) = ebetter(Fc(id(i)),Cc(id(i)),Fitness(p),Const(p),en); 
  end
  ecom0 = ecom==0; 
  ecom1 = ecom==1;
  if sum(ecom1)>=sum(ecom0), 
     eval=1; 
  else eval=0; 
  end
end
