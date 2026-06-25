% 
%
% MAIN_1_EX2: SIGNAL PROCESSING (1D SIGNALS)
%
% 

% clear environment
close all; clear all; clc; format short;

addpath("GENERATE_DATA_1D\")

% initialize pseudo-random numbers generators
rng('default'); rng(1);

% SHOW WHAT?
  SHOW_RESTORED_SIGNALS   = 1;

% set algorithms to be tested 
  TEST_DENOISING   = 0;  
  TEST_RESTORATION = 1;
  TEST_INPAINTING  = 0;

% For Tikhonov regularization: set values of the regularization parameter lambda to be tested
  lambdas_min     = 10^-4;     %10;
  lambdas_max     = 30;      %10^6;
  lambdas_n       = 500;      %500;
  lambdas         = linspace(lambdas_min,lambdas_max,lambdas_n);

% ---------------INITIALIZE----------------------

% Generate - and eventually visualize - the original 
% and corrupted (that is, blurred/noisy/missing) signal b
GENERATE_DATA_1D
N = length(b);

%--------------------------------------------------------------------------
% DENOISING %
%--------------------------------------------------------------------------
if ( TEST_DENOISING == 1 ) %signal 1D DO NOT ADD EXTRA NOISE, DO NOT ADD BLUR
    fprintf('\n\n\nSIGNAL DENOISING by TIKHONOV:\n');
    %
    % Smoothing (degree = 2)
    % D is the second-order difference matrix.
      e = ones(N, 1);
      
      % No boundary conditions assumption
        D = spdiags([e -2*e e], 0:2, N-2, N);

      % Reflexive boundary condition assumption
        %D = spdiags([e -2*e e],-1:1, N , N ); D(1,2) = 2; D(end,end-1) = 2;

    DT = D'; DTD = DT * D;

    % Solve the least square problem
      error = [];
      I = eye(size(D,2));
      for i  = 1 : length(lambdas)
          %
          X_TIKH = (I + lambdas(i) * DTD) \ b;
          rel_err = compute_rel_err(xtr,X_TIKH);
          error = [error;lambdas(i),rel_err];
          fprintf('\tThe error for lambda equal to %.5f is: %.5f\n',lambdas(i),rel_err)
      end
      [~,i] = min(error(:,2));
      lambda_optimal = (lambdas(i));
      fprintf('\nThe best value of lambda is %.5f\n',lambda_optimal);
      fprintf('The error in correspondence of %.5f is %.5f\n',lambda_optimal,min(error(:,2)));

      X_TIKH = (I + lambda_optimal * DTD) \ b;

    if ( SHOW_RESTORED_SIGNALS == 1 )
        
            picturewidth = 20;
            hw_ratio = 0.6; %Set for the height

            fig_1 = figure(200);
            clf(200)
            n = 1:N;
            plot(n, X_TIKH,'b','LineWidth', 1.5);hold on 
            plot(xtr,'r','LineWidth', 1.5)
            legend('Restored signal ','Original signal x','Location','northwest')
            title('\textbf{Smoothed signal}','Interpreter','latex'); axis tight;

            set(findall(fig_1,'-property','FontSize'),'FontSize',21);
            set(findall(fig_1,'-property','Box'),'Box','off');
            set(findall(fig_1,'-property','Interpreter'),'interpreter','latex');
            set(findall(fig_1,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex');
            set(fig_1,'units','centimeters','Position',[3 3 picturewidth hw_ratio * picturewidth]);
            pos = get(fig_1,'Position');
            set( fig_1,'PaperPositionMode','Auto','PaperUnits','centimeters','papersize',[pos(3),pos(4)] );
            print(fig_1,'pdf_figure','-dpdf','-vector','-fillpage');    
    end
end

%--------------------------------------------------------------------------
% RESTORATION %
%--------------------------------------------------------------------------
if ( TEST_RESTORATION == 1 )
    fprintf('\n\n\n SIGNAL RESTORATION:\n');

    % Derivative regularization (noisy)
      e = ones(N, 1);
      
      % No boundary conditions assumption
        D = spdiags([e -2*e e], 0:2, N-2, N);

      % Reflexive boundary condition assumption
        %D = spdiags([e -2*e e],-1:1, N , N ); D(1,2) = 2; D(end,end-1) = 2;


    % Solve least squares problem
      error = [];
     
      AT = A';  ATA = AT * A;
      DT = D';  DTD = DT * D;

      tau = 0.98;

      for i  = 1 : length(lambdas)
          %
          X_TIKH = (ATA + lambdas(i) * DTD) \ (AT*b);
          rel_err = compute_rel_err(xtr,X_TIKH);
          %gg = norm(A*X_TIKH - b,2)^2 / d;
          %error = [error;lambdas(i),gg, tau * n_sigma^2];
          error = [error;lambdas(i),rel_err];
          fprintf('\tThe error for lambda equal to %.5f is: %.5f\n',lambdas(i),rel_err)
      end
      [~,i] = min(error(:,2));
      lambda_optimal = lambdas(i);
      fprintf('\nThe best value of lambda is %.5f\n',lambda_optimal);
      fprintf('The error in correspondence of %.5f is %.5f\n',lambda_optimal,min(error(:,2)));

      X_TIKH = (ATA + lambda_optimal * DTD) \ (AT*b);

    %%%
    if ( SHOW_RESTORED_SIGNALS == 1 )
        %
        picturewidth = 20;
        hw_ratio     = 0.7;     %Set for the height

        fig_2 = figure(200);
        clf(200)
        plot(X_TIKH,'LineWidth',1.5); hold on; 
        plot(xtr,'r','LineWidth',1.5);
        legend('Restored ','Original signal','Location','northwest')
        xlabel('i'); ylabel('x(i)');
        title('\textbf{Smoothed signal}','Interpreter','latex' );
        axis tight;

        set(findall(fig_2,'-property','FontSize'),'FontSize',21);
        set(findall(fig_2,'-property','Box'),'Box','off');
        set(findall(fig_2,'-property','Interpreter'),'interpreter','latex');
        set(findall(fig_2,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex');
        set(fig_2,'units','centimeters','Position',[4 4 picturewidth hw_ratio * picturewidth])
        pos = get(fig_2,'Position');
        set( fig_2,'PaperPositionMode','Auto','PaperUnits','centimeters','papersize',[pos(3),pos(4)] )
        print(fig_2,'pdf_figure','-dpdf','-vector','-bestfit');
    end
end

%-----------------------------------------------------------------------
% INPAINTING ( %signal 1D s_id=11, DO NOT ADD NOISE, DO NOT ADD BLUR )
%-----------------------------------------------------------------------

if ( TEST_INPAINTING == 1 )
    % xtr original signal
    % b   corrupted signal (Missing data appear as NaN's in vector b)

    % Load of the signal for the inpainting, no blur and noise to be
    % added
      s_id = 11;
      s_file     = 'input/inpainting_data.mat';
      s_name     = 'SIGNAL_INPAINTING';
      load(s_file,'-mat');

      b = Axtr;
    
    % Signal dimension
      d = numel(xtr);

    % Plot of the original signal and the known samples

      picturewidth = 20;
      hw_ratio = 0.6; %Set for the height

      hfig_3 = figure(101);
      b_color = 'b'; xtr_color = 'r';

      x_min = min( [ min(b) min(b) min(b) ] );
      x_max = max( [ max(b) max(b) max(b) ] );
      x_range = x_max - x_min;
      x_min = x_min - 0.1 * x_range - 0.1;
      x_max = x_max + 0.1 * x_range + 0.1;

      hold on;
      plot(1:d,xtr,xtr_color,'LineWidth',1.5);
      plot(1:d,b,b_color,'LineWidth',1.5);
      axis ([1 d x_min x_max]);
      xlabel('i');ylabel('x(i)');grid on;
      title('\textbf{v (reconstructed samples)}, \textbf{b (input samples)}','Interpreter','latex');
      legend('v','b');

      set(findall(hfig_3,'-property','FontSize'),'FontSize',21);
      set(findall(hfig_3,'-property','Box'),'Box','off');
      set(findall(hfig_3,'-property','Interpreter'),'interpreter','latex');
      set(findall(hfig_3,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex');
      set(hfig_3,'units','centimeters','Position',[3 3 picturewidth hw_ratio * picturewidth]);
      pos = get(hfig_3,'Position');
      set( hfig_3,'PaperPositionMode','Auto','PaperUnits','centimeters','papersize',[pos(3),pos(4)] );
      print(hfig_3,'pdf_figure','-dpdf','-vector','-fillpage');     


      

    fprintf('\n\n\n SIGNAL INPAINTING:\n');
   
    % Estimate missing data by least squares:
    % Minimize the energy of second-order derivative subject to the data consistency constraint.

    % Define matrix D, D represents the second-order derivative (2nd-order difference).
    % D is defined as a sparse matrix so that Matlab
    % subsequently uses fast solvers for banded systems.
      e = ones(N, 1);
      D = spdiags([e -2*e e], 0:2, N-2, N);

    % Matrix D (only for visualization) 
    % Fist corner of D:
    % full(D(1:5, 1:5))
    % Last corner of D:
    % full(D(end-4:end, end-4:end))

    % Define matrices S and Sc
    k = isfinite(b);                    % k : logical vector, indexes known values

    S = speye(N);
    S(~k, :) = [];                      % S : sampling matrix

    Sc = speye(N);                      % Sc : complement of S
    Sc(k, :) = [];

    L = sum(~k)                         % L : number of missing values

    % Estimate missing data

    % Least square estimation of missing data.
    % Note that the system matrix is banded so the system
    % equations can be solved very efficiently with a fast banded system solver.
    % By defining S and D as sparse matrices, Matlab calls a fast
    % banded system solver by default.
    
     x = b;

    % Let's define b_1 the signal of only the known values
      b_1 = b(k);

    % We want to recover the signal v of the missing samples from the
    % original signal xtr
      DT  = D';
      ScT = Sc'; ST = S';

      v   = - inv(Sc*DT*D*ScT) * (Sc*DT*D*ST*b_1);

      X_INPAINTNG = ST * b_1 + ScT * v;
   
    if ( SHOW_RESTORED_SIGNALS == 1 )

        picturewidth = 20;
        hw_ratio = 0.6; %Set for the height

        hfig_2=figure(200);
        clf(200);
        n = 1:N;
        plot(n, x, 'k*', n(~k), v ,'r.','LineWidth',1);
        legend('Known data', 'Estimated data')

        set(findall(hfig_2,'-property','FontSize'),'FontSize',21);
        set(findall(hfig_2,'-property','Box'),'Box','off');
        set(findall(hfig_2,'-property','Interpreter'),'interpreter','latex');
        set(findall(hfig_2,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex');
        set(hfig_2,'units','centimeters','Position',[3 3 picturewidth hw_ratio * picturewidth]);
        pos = get(hfig_2,'Position');
        set( hfig_2,'PaperPositionMode','Auto','PaperUnits','centimeters','papersize',[pos(3),pos(4)] );
        print(hfig_2,'pdf_figure','-dpdf','-vector','-fillpage');     


        hfig_1 = figure(201);
        clf(201)
        plot(n, X_INPAINTNG,'b','LineWidth',1.5)
        title('\textbf{Final signal} x','interpreter','latex')
        legend('x')

        set(findall(hfig_1,'-property','FontSize'),'FontSize',21);
        set(findall(hfig_1,'-property','Box'),'Box','off');
        set(findall(hfig_1,'-property','Interpreter'),'interpreter','latex');
        set(findall(hfig_1,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex');
        set(hfig_1,'units','centimeters','Position',[3 3 picturewidth hw_ratio * picturewidth]);
        pos = get(hfig_1,'Position');
        set( hfig_1,'PaperPositionMode','Auto','PaperUnits','centimeters','papersize',[pos(3),pos(4)] );
        print(hfig_1,'pdf_figure','-dpdf','-vector','-fillpage');       

        hfig = figure(202);
        clf(202)
        n = 1:N;
        plot(n, x,'b','LineWidth',1.5)
        legend('b')
        title('\textbf{INPUT SAMPLES} b','interpreter','latex');

        set(findall(hfig,'-property','FontSize'),'FontSize',21);
        set(findall(hfig,'-property','Box'),'Box','off');
        set(findall(hfig,'-property','Interpreter'),'interpreter','latex');
        set(findall(hfig,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex');
        set(hfig,'units','centimeters','Position',[3 3 picturewidth hw_ratio * picturewidth]);
        pos = get(hfig,'Position');
        set( hfig,'PaperPositionMode','Auto','PaperUnits','centimeters','papersize',[pos(3),pos(4)] );
        print(hfig,'pdf_figure','-dpdf','-vector','-fillpage');
        
    end
end
