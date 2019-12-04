alpha=1/3; 
beta = .99;  
delta=.025; 
sigma = 2;
roe=0.5; 
sge=0.2; 


num_grid=5;
[z_grid, Pi]=TAUCHEN(num_grid,roe,sge,3);
z_grid=exp(z_grid');


%dist
PI_inv=Pi^1000;
PI_inv=PI_inv(1,:)';
N_s=z_grid*PI_inv; 
a_low = 0; 
a_high = 100;
numa = 50;
na=500;

a = linspace(a_low, a_high, numa); 
iva = linspace(a_low, a_high, na);




Kmin=20;
Kmax=50;
tolerk=1;
while abs(tolerk)>.01
    if Kmax-Kmin<0.00001
       break
   end
    
    K_guess=(Kmin+Kmax)/2;
    Interests= alpha*K_guess^(alpha-1)*N_s^(1-alpha)+(1-delta);
    Wage=(1-alpha)*K_guess^alpha*N_s^(-alpha);
    
   	cons = bsxfun(@minus, Interests* a', a);
    cons = bsxfun(@plus, cons, permute(z_grid, [1 3 2])*Wage);
    re = (cons .^ (1-sigma)) ./ (1 - sigma); 
    re(cons<0)=-Inf;
    
    v_guess = zeros(num_grid, numa);
    
    
    v_tol = 1;
    while v_tol >10e-5
        
	   value_mat=re+beta*repmat(permute((Pi*v_guess),[3 2 1]), [numa 1 1]);
       [vfn, pol_index] = max(value_mat, [], 2);
       vfn=permute(vfn, [3 1 2]);
       v_tol = max(abs(vfn-v_guess));
       v_tol = max(v_tol(:));
       
       v_guess = vfn;
    end
    
    pol_index=permute(pol_index, [3 1 2]);
    pol_fn = a(pol_index);
    
    
    V=zeros(num_grid,na);
    for i=1:num_grid
        V(i,:)=interp1(a,v_guess(i,:),iva);
    end
   
    
    
    
    cons = bsxfun(@minus, Interests* iva', iva);
    cons = bsxfun(@plus, cons, permute(z_grid, [1 3 2])*Wage);
   
    
    re = (cons .^ (1-sigma)) ./ (1 - sigma); % current period utility
    re(cons<0)=-Inf;
    
    value_mat=re+beta*repmat(permute((Pi*V),[3 2 1]), [na 1 1]);
    [vfn, pol_index] = max(value_mat, [], 2);
    vfn=permute(vfn, [3 1 2]);
    pol_index=permute(pol_index, [3 1 2]);
    pol_fn = iva(pol_index);
    
    mu=zeros(num_grid,na);
    mu(:)=1/(na*num_grid);
    
    dist=1;
    
    %iteration
  while dist>10e-6 
      muprime = zeros(size(mu));
     [z_ind, a_ind, mass] = find(mu); 
    
    for ii = 1:length(z_ind)
        apr_ind = pol_index(z_ind(ii),a_ind(ii)); 
        
        muprime(:, apr_ind) = muprime(:, apr_ind) +(Pi(z_ind(ii), :)*mass(ii))';
        
    end
    dist = max(max(abs(mu-muprime)));
    mu=muprime;
  end
   
   K=sum(sum(mu.*pol_fn));
   tolerk=K-K_guess;
   if tolerk>0;
       Kmin=K_guess;
   else Kmax=K_guess;
   end
end

% policy function
figure(1)
title(['Policy Function'])
plot(iva,pol_fn)
z_name=cellstr(num2str(z_grid'));


popul=reshape(mu',[num_grid*na,1]);
wealth=reshape(repmat(iva,num_grid,1)',[num_grid*na,1]);

mu=sum(mu);
figure(2)
title('Assets Dist')
bar(iva,mu)



%lorenz curve%%%%
Wealth=sortrows([wealth,popul,popul.*wealth]);
Wealth=cumsum(Wealth);
pw=Wealth(:,2);
pw=pw(end);
Wealth(:,2)=Wealth(:,2)/pw;
w=Wealth(:,3);
w=w(end);
Wealth(:,3)=Wealth(:,3)/w;
giniwealth2 = 1 - sum((Wealth(1:end-1,3)+Wealth(2:end,3)) .* diff(Wealth(:,2)));

figure(3)
title('Lorenz Curve' )
area(Wealth(:,2),Wealth(:,3),[0.5,0.5,1.0])
hold on
plot([0,1],[0,1])
axis square
hold off


y=Wage*z_grid;
A=repmat(a,[num_grid,1]);
Y=repmat(y',[1,na]);

C=Y+Interests*A-pol_fn;
CF=C(:,pol_index');
cfprime=reshape(CF,[num_grid na num_grid]);

i=0;

while i < num_grid+1

    c1(i,:)=Pi(i,:)*cfprime(:,:,i);
i=i+1;

end




EE=sum(sum(abs(C.^(-sigma)-beta*c1.^(-sigma)*Interests).*mu));