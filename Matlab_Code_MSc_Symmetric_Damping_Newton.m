function Diode_Slides_Symmetric_Damping
%% Program starts

%% Parameters list
 para.epsilon = 11.7 * 8.85 * 1e-14;
       para.k = 1.38e-23;
       para.q = 1.6021e-19;
       para.T = 300;
   para.n_int = 1.22*1e10;

%% Grid and Variables
%  Define the mesh
% 
% $$-5\times10^{-4}=x_0<x_1<...<x_m<x_{m+1}=5\times10^{-4}$$
%  
%

%grid
para.M = 199; % #points
para.a = -5e-4; % left boundary
para.b =  5e-4; % right boundary
para.h = (para.b-para.a)/(para.M+1); % 
     
%variables
    d_psi = 100 + ones(para.M,1);
      psi = zeros(para.M,1); 
        f = zeros(para.M,1);
psi_left  = -para.k*para.T/para.q*log(1e18/para.n_int);
psi_right =  para.k*para.T/para.q*log(1e18/para.n_int);



%% Initial guess
% 
% Charge netrual guess: 
%
% $$\psi(x) = \frac{kT}{q}\textnormal{arsinh}(\frac{D(x)}{2n_i}).$$
%
% linear initial guess:
% 
% $$\psi(x) =  \frac{x}{x_N - x_0}(\psi_{M+1}-psi_0) $$
%


for i = 1 : para.M

    %Charge neutral guess
%    psi(i) = para.k*para.T/para.q*asinh(0.5*Dop(para.a+i*para.h)/para.n_int);
    
    %Charge neutral guess with a perturbation 
    %psi(i) = para.k*para.T/para.q*asinh(0.5*Dop(para.a+i*para.h)/para.n_int)-0.1;
    
    %x=-0.5 if x<0 else x=0.5
%     if(i<para.M/2)
%         psi(i) = para.k*para.T/para.q*asinh(0.5*Dop(para.a+i*para.h)/para.n_int)-0.5;
%     else
%         psi(i) = para.k*para.T/para.q*asinh(0.5*Dop(para.a+i*para.h)/para.n_int)+0.5;
%     end


    %Linear initial guess
    psi(i) = (psi_right - psi_left)/(para.b-para.a)*(i*para.h+para.a);
    
    %Linear initial guess with a fliping 
    %psi(i) = (psi_left - psi_right)/(para.b-para.a)*(i*para.h+para.a);
    
    %zero guess
    %psi(i) = 0;
    
    %sin(i*pi)
    %psi(i) = 0.5*sin(pi/para.b*(i*para.h+para.a)); 
        
    %exp(x)
    %psi(i) = exp(i*para.h+para.a); %almost linear
    
    %tan(x)
    %psi(i) = tan(i*para.h+para.a); %almost linear
    
end
%plot(psi);

%% General damping Newton's method
%
% For the damped Newton's method, the update setp is replaced by
%
% $$\psi^{k+1} = \psi^{k} + \lambda d\psi^k,$$
%
% where the damping parameter $\lambda$ is needed to be determined. The
% above expression implies that, after we have calculated the Newton
% correction, we determine the scalar factor $\lambda$ such that 
% the $\lambda d\psi^k$ is the best possible update. One way to do that is 
% introducing the following functional:
% 
% $$g(\lambda) = \frac{1}{2}(F(\psi^{k} + \lambda d\psi^k),
%   F(\psi^{k} + \lambda d\psi^k)) = \frac{1}{2} 
%   ||F(\psi^k+\lambda d\psi^k)||^2, $$
%
% where the $(\cdot,\cdot)$ is the normal inner product in n-dimensional space.
%
% The minimum $g$ with respect to $\lambda$ means $F(\psi^{k} + \lambda d\psi^k)$ close to zero as
% much as possible, i.e., $F(\psi^{k+1})\rightarrow 0$. Remeber that $\psi$ satisfying
% $F(\psi)=0$ is the solution of original system.
% 
% For obtaining the $\lambda$, we take the derivative of this functional with
% respect to the damping parameter:
%
% $$g'(\lambda) = (F'(\psi^{k} + \lambda d\psi^k)d\psi^k,
%   F(\psi^{k} + \lambda d\psi^k)),$$
%
% and sovling $\lambda$ by let it equal to zero.
%
% In addition, we can check that
%
% $$g'(0) = (F'(\psi^{k})d\psi^k,
%   F(\psi^{k}))=-(F(\psi^{k}),F(\psi^{k}))\leq 0,$$
%
% where the final equality is obtained by using the definition of the
% Newton correction. This calculation shows that the Newton correction
% provides a descent direction for the functional g, at least starting from
% 0. *So for small values of the damping parapmeter, the value of g will
% decrease*.
%
% In conclusion, for small values of the damping parameter,
%
% $$F(\psi^{k} + \lambda d\psi^k), F(\psi^{k} + \lambda d\psi^k) 
%   <F(\psi^{k}), F(\psi^{k}).$$
%
% A solution of the original problem F(u)=0 necessarily leads to a minimum
% of $||F(\psi^k+\lambda d\psi^k)||^2$. Vice versa, if we find an iterate,
% a corresponding Newton correction and a damping parameter such that the
% value of $g$ is zero, then we know the following vector is a solution of
% the non-linear problem:
%
% $$\psi^k+\lambda d\psi^k$$
%
% In total, the main idea of damping is that in each step, we determine a damping
% parameter such that the value of the residual decreases.
% 
% 1. First solve $F'(\psi^k)d\psi^k=-F(\psi^k)$,
%
% 2. Find $\lambda$: solve $\min g(\lambda) =
% \frac{1}{2}||F(\psi^k+\lambda d\psi^k)||^2$,
%
% 3. updata $\psi^{k+1} = \psi^{k} + \lambda d\psi^k$.
%
% Remark: the second step may use the Newton's method again.
%

%% Symmetric damped Newton's method
%
% When the original problem is symmetric, then we can derivate a better
% damping method. 
%
% Symmetry means that the elements of the Jacobian matrix satisfy  
% $\frac{\partial f_i}{\partial \psi_j} = \frac{\partial f_j}{\partial \psi_i}$.
%
% Followings are two importmant lemmas and a corollary which are the
% foundation of symmetric cases.
%
% *Lemma*  Let the functions $f_1, \cdots, f_N: R^N \rightarrow R$ are
% symmetic. Assume further that $f_i$, considered as a function of $x_i$
% only, is in $L^1(R)$. Then there exists a real-valued function $V_N$ such
% that 
%
% $$f_i(\psi_i, \cdots, \psi_N)=\frac{\partial V_N(\psi_i, \cdots, \psi_N)}{\partial \psi_i}$$
%
% for all $i = 1, 2, \cdots, N$.
%
% *Corollary* Let $F: R^N\rightarrow R^N$ have components $f_1, \cdots,
% f_N$ which satisfy the requirement in above Lemma. Then there exists a
% real-valued function $V_N: R^N \rightarrow R$ such that 
%
% $$F^T(\psi) = \nabla V_N(\psi).$$
%
% From this corollary it follows that, if $F(\psi^*)=0$, we also have that
% $\nabla V_N(\psi^*)=0$. Thus, $\psi^*$ is a stationary point of $V_N$.
%
% Hence, finding the zero of the function $F$ is equaivalent to finding a
% minimum of the constructed function. However, the solution may also be a
% maximum or some type of stationary point, such as a saddle point. But the
% following lemma tells us that with extra conditions of $F$, $V_N$ is a
% convex function and has a unique minimum at $\psi^*$.
% 
% *Lemma* Let $F: R^N\rightarrow R^N$ be such that $F'(\psi)$ is symmetic
% and 
%
% $$v^TF'(\psi)v\geq kv^Tv \quad \forall \psi v \in R^N, $$
%
% for some constant $k$. Then the function $V_N$ provided by Lemma above is
% strictly convex on $R^N$, i.e.,
%
% $$V_N(\lambda \psi + (1-\lambda)v)<\lambda V_N(\psi)+(1-\lambda)V_N(v), 
% \quad \forall \lambda \in (0,1), \forall \psi,v \in R^N.$$
%
% Furthermore, 
%
% $$V_N(\psi)\rightarrow +\infty \quad as \quad ||\psi||\rightarrow \infty.$$
%
% If $F$ satisfies the conditions in the corollary, and if $F'(\psi)$ is
% continuous, $F$ is termed *uniformly monotone*.
%
% Hence, for symmetric case, we replace 
%
% $$g(\lambda) =\frac{1}{2}||F(\psi^k+\lambda d\psi^k)||^2,$
%
% by 
%
% $$g(\lambda)=V_N(\psi^k+\lambda d\psi^k).$$
%
% We have:
%
% $$g'(\lambda) = (F(\psi^k+\lambda d\psi^k), d\psi^k).$$
%
% Hence,
%
% $$g'(0) = (F(\psi^k), d\psi^k) $$
%
% showing that indeed that Newton direction provides descent.
%
% Another nice property is the following:
%
% $$g''(\lambda) = (F'(\psi^k+\lambda d\psi^k)d\psi^k, d\psi^k)>0$$
%
% due to the positive definiteness of $F'.$ This implies that the 
% functional is *strictly convex*! This can then be used to find the *unique 
% minimum* of the functional.
%
% An important question is whether the damped Newton procedure discussed above will
% actually lead to the solution of the original problem $F(u)=0$.
%
% In the general case, there ws some doubt about the convergence since 
% we could not guarantee that
%
% $$g^{k+1}(0)\leq g^{k}(\lambda^k),$$
%
% where the supscript $k$ means the iteration step.
%
% This problem does not occur in the present case, since
%
% $$g^{k+1}(0)=V_N(\psi^{k+1}).$$
%
% so that certainly
%
% $$g^{k+1}(\lambda^{k+1})<g^{k}(\lambda^{k}).$$
%
% In other words: a sequence $\psi^0, \psi^1, \psi^2, \cdots$ is generated
% which satisfies
%
% $$V_N(\psi^*)\leq V_N(\psi^{k+1}) \leq V_N(\psi^{k}),$$
%
% for all $k\geq0$. From the bounded monotone convergence theorem it then
% follows that the sequence of $V_N(\psi^k)$ has a limit.
%
%% Algorithm: Symmetric damped Newton's method
%% 1. First get $d\psi^k$
%
% by $F'(\psi^k)d\psi^k=-F(\psi^k)$
%
Mat = diag(ones(para.M-1,1)*(-para.epsilon/para.h),1) + ...
      diag(ones(para.M-1,1)*(-para.epsilon/para.h),-1);
flag = 0;
max_err =[];
while max(abs(d_psi)) > 1e-16
    %iterative matrix
    for i = 1 : para.M
        Mat(i,i) = 2*para.epsilon/para.h + ...
            para.q*para.n_int*( para.q/(para.k*para.T)*exp(-para.q*psi(i)/(para.k*para.T))+...
            para.q/(para.k*para.T)*exp(para.q*psi(i)/(para.k*para.T)) );
    end
    %RHS of iterative formula
    for i = 1 : para.M
        if i == 1
            f(i) = - ( para.epsilon/para.h*(-psi_left + 2*psi(i) - psi(i+1)) - ...
                para.q*(  para.n_int*exp(-para.q*psi(i)/(para.k*para.T)) - ...
                para.n_int*exp(para.q*psi(i)/(para.k*para.T)) + Dop(para.a+i*para.h) ) ); 
        elseif i == para.M
            f(i) = - ( para.epsilon/para.h*(-psi(i-1) + 2*psi(i) - psi_right) - ...
                para.q*(  para.n_int*exp(-para.q*psi(i)/(para.k*para.T)) - ...
                para.n_int*exp(para.q*psi(i)/(para.k*para.T)) + Dop(para.a+i*para.h) ) );
        else     
            f(i) = - ( para.epsilon/para.h*(-psi(i-1) + 2*psi(i) - psi(i+1)) - ...
                para.q*(  para.n_int*exp(-para.q*psi(i)/(para.k*para.T)) - ...
                para.n_int*exp(para.q*psi(i)/(para.k*para.T)) + Dop(para.a+i*para.h) ) );
        end
    end
    
    %1. obtaining the d\psi^k
    d_psi = Mat\f;

%% 2. Find $\lambda$ 
%
% solve $\min g(\lambda)=V_N(\psi^k+\lambda d\psi^k),$
% 
% by $g'(\lambda)=0,$
%
% where $g'(\lambda) = (F(\psi^k+\lambda d\psi^k), d\psi^k),$ and 
% $F=[F_1, F_2, \cdots, F_M]^T$ for $M$ gird points, namely $g'(\lambda)$
% is the dot product of $F=[F_1, F_2, \cdots, F_M]^T$ and $d\psi^k=[d\psi^k_1, d\psi^k_2, \cdots, d\psi^k_M].$
%
%
% Some guidelines:
% 1. It does not seem wise to determine the damping parameters in a Newton 
% iteration to extreme accuracy; finding the minimum approximately is sufficient
%
% 2. As we wish to make use, or benefit, from the quadratic convergence of Newton’s method 
% when near to the solution, it seems wise to check whether $\lambda=1$ is a suitable candidate.
%
% This is supported by the fact that the quadratic polynomial
% 
% $$p(\lambda) = g(0) + (\lambda-\frac{1}{2}\lambda^2)g'(0), $$
%
% agrees with the damping functional and its first derivative at 0. Since
%
% $$p'(\lambda) = (1-\lambda)g'(0),$$
%
% the minimum of $p^l$ is found for $\lambda=1$. Thus, we have the third
% guideline;
%
% 3. When starting the search, first calculate $g'(1)$.
%
% Indeed: we know that the derivative at 0 is negative. Also, that we have a convex functional.
% Hence, if the derivative at 1 is still negative, we can start from there, the minimum will be beyond 1.
% If the derivative at 1 is positive, then we know that the minimum is somewhere between 0 and 1.
%

    %2. solve the lambda 
      %sum_g is the dot product of (g^k)'(lambda) on page 27 slide 07
      %we want it equal to 0, then using Newton iteration to get lambda.
    if abs(sum_g(1, psi, d_psi, para)) < 1e-6
        lambda = 1;
        %fprintf("directly 1\n");
    else
        if sum_g(1, psi, d_psi, para) < 0
            lambda = 1;
        else
            lambda = 0;
        end
%        err_lambda = sum_g(lambda, psi, d_psi, para);
        G = @(lambda) sum_g(lambda, psi, d_psi, para);
        options = optimoptions('fsolve','TolX',1e-15);
        lambda = fsolve(G, lambda, options);
%         while abs(err_lambda) > 1e-6
%             %Nested Newton's method
%             % λ = λ - g'(λ)/g''(λ);
%             %sum_dg:= g''(λ)
% 
%             %lambda = lambda - sum_g(lambda, psi, d_psi, para)/sum_dg(lambda, psi, d_psi, para);
%             %err_lambda = sum_g(lambda, psi, d_psi, para);
%         end
    end
    
    %draw the damping functional graph
%     figure(flag+1)
%     points = 0:0.01:5;
%     plot_g = sum_g(points, psi, d_psi);
%     plot(points, plot_g);
%     fprintf('lambda is %f\n', lambda);

%% 3. updata $\psi^{k+1} = \psi^{k} + \lambda d\psi^k$.  
    %3. update
    psi = psi + lambda * d_psi;
    
    flag = flag + 1;
    max_err(flag) = max(abs(d_psi));
end

%% figure
%% Draw figures

% plot(a:h:b,[psi_left; psi; psi_right] ,'b', 'LineWidth', 1.5);
% ylabel('$\psi$', 'Interpreter', 'latex', 'FontSize', 20);
% xlabel('$x$', 'Interpreter', 'latex', 'FontSize', 20);
% set(gcf,'color','w');
% set(gca,'linewidth',1)
% print( '-dpng', '-r300','~\Desktop\Master_Thesis\Thesis_Material\DiSol_CT');

%title('Diode-FVM-CorrTrans');
fprintf('It takes %d iterative step(s)\n', flag);
figure
plot(max_err, 'b', 'LineWidth', 1.5);
set(gca, 'YScale', 'log')
set(gcf,'color','w');
set(gca,'linewidth',1)
xlabel('iteration steps', 'Interpreter', 'latex', 'FontSize', 20);
ylabel('max correction $|d\psi|$','Interpreter', 'latex', 'FontSize', 20);
%print( '-dpng', '-r300','~\Desktop\Master_Thesis\Thesis_Material\DiIter_SD_linear');
%title('Diode-FVM-CorrTrans');
end

%%
function D = Dop(x)

  if x < 0
      D = -10^18;
  elseif x == 0
      D = 0;
  else
      D = 10^18;
  end
  
end
%%
function g = sum_g(lambda, psi, d_psi, para)
    
    psi_left  = -para.k*para.T/para.q*log(1e18/para.n_int);
    psi_right =  para.k*para.T/para.q*log(1e18/para.n_int); 
  
    g = 0;
    i = 1;
    g = g + (para.epsilon/para.h*(-psi_left+...
            2*(psi(i)+lambda*d_psi(i))-(psi(i+1)+lambda*d_psi(i+1)))-...
            para.q*( para.n_int*exp(-para.q*(psi(i) + lambda*d_psi(i))/(para.k*para.T))-...
            para.n_int*exp(para.q*(psi(i) + lambda*d_psi(i))/(para.k*para.T))+Dop(para.a+i*para.h)))...
            *d_psi(1);
        
    for i = 2 : para.M-1
        g = g + (para.epsilon/para.h*(-(psi(i-1)+lambda*d_psi(i-1))+...
            2*(psi(i)+lambda*d_psi(i))-(psi(i+1)+lambda*d_psi(i+1)))-...
            para.q*( para.n_int*exp(-para.q*(psi(i) + lambda*d_psi(i))/(para.k*para.T))-...
            para.n_int*exp(para.q*(psi(i) + lambda*d_psi(i))/(para.k*para.T))+Dop(para.a+i*para.h)))...
            *d_psi(i);
    end
    
    i = para.M;
    g = g + (para.epsilon/para.h*(-(psi(i-1)+lambda*d_psi(i-1))+...
            2*(psi(i)+lambda*d_psi(i))-psi_right)-...
            para.q*( para.n_int*exp(-para.q*(psi(i) + lambda*d_psi(i))/(para.k*para.T))-...
            para.n_int*exp(para.q*(psi(i) + lambda*d_psi(i))/(para.k*para.T))+Dop(para.a+i*para.h)))...
            *d_psi(length(psi));
end
%%
function d_g = sum_dg(lambda, psi, d_psi, para)
    d_g = 0;
    
    i = 1;
    d_g = d_g + (para.epsilon/para.h*(0 + 2*d_psi(i) - d_psi(i+1))...
            - para.q*(-para.n_int*para.q/(para.k*para.T)*d_psi(i)*exp(-para.q*(psi(i)+lambda*d_psi(i))/(para.k*para.T))...
            -para.n_int*para.q/(para.k*para.T)*d_psi(i)*exp(para.q*(psi(i)+lambda*d_psi(i))/(para.k*para.T))))*d_psi(i);
        
    for i = 2 : para.M-1
        d_g = d_g + (para.epsilon/para.h*(-d_psi(i-1) + 2*d_psi(i) - d_psi(i+1))...
            - para.q*(-para.n_int*para.q/(para.k*para.T)*d_psi(i)*exp(-para.q*(psi(i)+lambda*d_psi(i))/(para.k*para.T))...
            -para.n_int*para.q/(para.k*para.T)*d_psi(i)*exp(para.q*(psi(i)+lambda*d_psi(i))/(para.k*para.T))))*d_psi(i);
    end
    
    i = para.M;
    d_g = d_g + (para.epsilon/para.h*(-d_psi(i-1) + 2*d_psi(i) - 0)...
            - para.q*(-para.n_int*para.q/(para.k*para.T)*d_psi(i)*exp(-para.q*(psi(i)+lambda*d_psi(i))/(para.k*para.T))...
            -para.n_int*para.q/(para.k*para.T)*d_psi(i)*exp(para.q*(psi(i)+lambda*d_psi(i))/(para.k*para.T))))*d_psi(i);
end










