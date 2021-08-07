function Diode_Correction_Trans
%% Problem description
% In simulating semiconductor devices, a famous drift-diffusion model reads
%
% $$-\nabla \cdot (\epsilon\nabla\psi) = q(p-n+D) $$
%
% Rewrite the carrier concentrations in terms of the quasi-Fermilevels, the
% equation becomes non-linear:
%
% $$-\nabla \cdot (\epsilon\nabla\psi) = q[n_i \exp(\frac{q(\phi_p-\psi)}{kT})
%                                          -n_i
%                                          \exp(\frac{q(\psi-\phi_n)}{kT})+D],
%                                          $$
%
% with zero applied bias boundary conditions
%
% $$\psi(x=-5\cdot10^{-4}) = -\frac{kT}{q}\log\frac{10^{18}}{n_i} $$ 
%
% $$\psi(x=5\cdot10^{-4}) = \frac{kT}{q}\log\frac{10^{18}}{n_i} $$
%
% Denoting the right hand side by
%
% $$ h(\psi) = q[n_i \exp(\frac{q(\phi_p-\psi)}{kT})
%                                          -n_i
%                                          \exp(\frac{q(\psi-\phi_n)}{kT})+D].
%                                          $$
%
% we can check that $h'(\psi)>0$ is trivial. Hence, the non-linear Poisson
% euqation is eplliptic.
% Next step, one should scale the equations, since the problem area is
% extremely small(micrometers), and the size of coefficients differs a lot
% in size
%
% After that, we discretize the non-linear Poisson equation using the finite
% volume method, and use the midpoint rule for the zero-order term,
% obtaining the following discrete equation:
%
% $$ \epsilon_{i-1/2}\frac{\psi_{i} - \psi_{i-1}}{x_{i} - x_{i-1}}
%   -\epsilon_{i+1/2}\frac{\psi_{i+1} - \psi_{i}}{x_{i+1} - x_{i}}    
%   +q[-n_i \exp(\frac{q(\phi_{p,i}-\psi_i)}{kT})
%      +n_i \exp(\frac{q(\psi_i-\phi_{n,i})}{kT})+D(x_i)]=0. $$
%                
% For solving this equation, we focus on generic Newton's iteration at first.
% Assume $\phi_p$ and $\phi_n$ are 0, and $\epsilon$ is constant$.
%
% Let 
%
% $$ f_i(\psi) = \frac{\epsilon}{h}(-\psi_{i-1}+2\psi_{i}-\psi_{i+1}) + q(n(\psi_i)-p(\psi_i)+D(x_i)), $$
%
% where 
%
% $$ n(\psi_i) = n_i\exp(\frac{-q\psi_i}{kT}), \\\ p(\psi_i) =
% n_i\exp(\frac{q\psi_i}{kT}). $$
%
% Using the Taylor expansion, we get the Newton's iteration, that
%
% $$f_i(\psi^l) + \sum_{j=1}^{m}\frac{\partial f_i}{\partial \psi_j}\bigg|_{\psi^{(l)}}\psi_j^{(l)}=0 $$
%
% where $l=0,1,2,\cdots$ refers to the number of iteration. It can be
% written as 
%
% $$ \frac{\epsilon}{h}(-\Delta \psi_{j-1}^{(l)} + 2\Delta \psi_{j}^{(l)} - 
%     \Delta \psi_{j+1}^{(l)}) + \frac{\partial q(n(\psi_i)-p(\psi_i)+D(x_i))}{\partial \psi_{j}}\Delta \psi_{j}^{(l)} = - f_i(\psi^l)$$
%
% Let $M_j^l = \frac{2\epsilon}{h} + \frac{\partial
% q(n(\psi_i)-p(\psi_i)+D(x_i))}{\partial \psi_{j}}$, then the above
% equation can be written as
% 
% <<MATRIX.PNG>>
% 
% Finally, $\psi^{l+1} = \psi^{l} + \Delta \psi^{l}$.
%
%
% Next step, we consider the *symmetric damping Nemton's method*
% and the *correction transformation method*.
%
%% Program starts

%% Parameters list
epsilon = 11.7 * 8.85 * 1e-14;
      k = 1.38e-23;
      q = 1.6021e-19;
      T = 300;
  n_int = 1.22*1e10;

  
%% Grid and Variables
%  Define the mesh
% 
% $$-5\times10^{-4}=x_0<x_1<...<x_m<x_{m+1}=5\times10^{-4}$$
%  
%

%grid
M = 199; % #points
a = -5e-4; % left boundary
b =  5e-4; % right boundary
h = (b-a)/(M+1); % 
     
%variables
    d_psi = 100 + ones(M,1);
      psi = zeros(M,1); 
  d_v_psi = zeros(M,M);
        v = zeros(M,1);
psi_left  = -k*T/q*log(1e18/n_int);
psi_right =  k*T/q*log(1e18/n_int);

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


for i = 1 : M
    
    %Charge neutral guess
    %psi(i) = k*T/q*asinh(0.5*Dop(a+i*h)/n_int);  
    
    %Charge neutral guess with a perturbation 
    %psi(i) = k*T/q*asinh(0.5*Dop(a+i*h)/n_int)-0.1;
    
    %x=-0.5 if x<0 else x=0.5
%     if(i<M/2)
%         psi(i) = k*T/q*asinh(0.5*Dop(a+i*h)/n_int)-0.5;
%     else
%         psi(i) = k*T/q*asinh(0.5*Dop(a+i*h)/n_int)+0.5;
%     end
    
    %Linear initial guess
    psi(i) = (psi_right - psi_left)/(b-a)*(i*h+a);
    
    %Linear initial guess with a fliping
    %psi(i) = (psi_left - psi_right)/(b-a)*(i*h+a);
    
    %zero guess
    %psi(i) = 0;
    
    %sin(i*pi)
    %psi(i) = 0.5*sin(pi/b*(i*h+a));
    
    %exp(x)
    %psi(i) = exp(i*h+a); %almost linear
        
    %tan(x)
    %psi(i) = tan(i*h+a); %almost linear
    
end
%Initial guess graph
% plot(linspace(a,b,M),  psi ,'b', 'LineWidth', 1.5);
% %title('Boundary layer')
% ylabel('$\psi^0(x)$', 'Interpreter', 'latex', 'FontSize', 20);
% xlabel('$x$', 'Interpreter', 'latex', 'FontSize', 20);
% set(gcf,'color','w');
% set(gca,'linewidth',1)
% print( '-dpng', '-r300','~\Desktop\Master_Thesis\Thesis_Material\Ini_Diode_tan');

%% Correction transformation
% 1. Choose an initial guess $\psi^0$ and put $k=0$,
%
% 2. Evaluate $du^k$ by solving $\nabla_{\psi}F(\psi^k)du^k=-F(\psi^k)$,
%
% 3. Transform the correction $d\psi^k$ to $dv^k$ and do $dv^k =
% \nabla_{\psi}v(\psi^k)d\psi^k$,
%
% 4. Construct a new approximation $v^{k+1} = v^k+dv^k$,
%
% 5. Determine $\psi^{k+1}$ such that $v(\psi^{k+1}) = v^{k+1}$,
%
% 6. Increase $k$ by 1 and go to step 2 until convergence.
%
% In this case
%
% $$ f_i(\psi) = \frac{\epsilon}{h}(-\psi_{i-1}+2\psi_{i}-\psi_{i+1}) + q(n(\psi_i)-p(\psi_i)+D(x_i)). $$
%
% Then, we take 
%
% $$v=q(n(\psi_i)-p(\psi_i)+D(x_i)),$$ 
%
% such that we do not need to solve a nonlinear system again, and we can rewrite it as 
%
% $$v = q(2n_{i}\textnormal{sinh}(\frac{q}{kT}\psi_i)-D(x))$$
%
%
Mat = diag(ones(M-1,1)*(-epsilon/h),1) + diag(ones(M-1,1)*(-epsilon/h),-1);
flag = 0;
max_err =[];
%initial v = v(psi_0) 
for i = 1 : M   
    %v(i) = 2*q*n_int*sinh(q*psi(i)/(k*T)) - q * Dop(a+i*h) ;
    v(i) = 2*q*n_int*sinh(q*psi(i)/(k*T));
end
while max(abs(d_psi)) > 1e-12
    %iterative matrix
    for i = 1 : M
        Mat(i,i) = 2*epsilon/h + q*n_int*( q/(k*T)*exp(-q*psi(i)/(k*T)) + q/(k*T)*exp(q*psi(i)/(k*T)) );
    end
    
    %RHS of the iterative formula
    f = fun(psi, n_int, q, k, T, epsilon, h, psi_left, psi_right, a);
    %CT on page 24 slides 08
    
    %2. Evaluate d_psi
    d_psi = Mat\f;
    
    %3. Transform the correction d_psi to d_v
    for i = 1 : M
        %dv/dpsi
        d_v_psi(i,i) = 2*q^2*n_int/(k*T)*cosh(q*psi(i)/(k*T));
    end
    d_v = d_v_psi * d_psi;
    
    %4. Construct a new approximation
    v = v + d_v;
    
    %5. Determine psi(k+1) such that v(psi_k+1) = v_k+1
    for i = 1 : M    
        %In this case, the inverse function is explicit. 
        %We do not need a nested Newton's iteration
        %psi(i) = asinh((v(i) + q*Dop(a+i*h))/(2*q*n_int))*k*T/q;
        psi(i) = asinh(v(i)/(2*q*n_int))*k*T/q;
    end
    
    %psi = psi + d_psi;
    flag = flag + 1;
    max_err(flag) = max(abs(d_psi));
end



%% Draw figures

plot(a:h:b,[psi_left; psi; psi_right] ,'b', 'LineWidth', 1.5);
ylabel('$\psi(x)$', 'Interpreter', 'latex', 'FontSize', 20);
xlabel('$x$', 'Interpreter', 'latex', 'FontSize', 20);
set(gcf,'color','w');
set(gca,'linewidth',1)
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
%print( '-dpng', '-r300','~\Desktop\Master_Thesis\Thesis_Material\DiIter_CT_linear');
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

function f = fun(psi, n_int, q, k, T, epsilon, h, psi_left, psi_right, a)
    M = length(psi);
    f = zeros(M,1);
        i = 1;
            f(i) = - ( epsilon/h*(-psi_left + 2*psi(i) - psi(i+1)) - ...
                q*(  n_int*exp(-q*psi(i)/(k*T)) - n_int*exp(q*psi(i)/(k*T)) + Dop(a+i*h) ) ); 
        for i = 2:M-1    
            f(i) = - ( epsilon/h*(-psi(i-1) + 2*psi(i) - psi(i+1)) - ...
                q*(  n_int*exp(-q*psi(i)/(k*T)) - n_int*exp(q*psi(i)/(k*T)) + Dop(a+i*h) ) );
        end
        i = M;
            f(i) = - ( epsilon/h*(-psi(i-1) + 2*psi(i) - psi_right) - ...
                q*(  n_int*exp(-q*psi(i)/(k*T)) - n_int*exp(q*psi(i)/(k*T)) + Dop(a+i*h) ) );
end





