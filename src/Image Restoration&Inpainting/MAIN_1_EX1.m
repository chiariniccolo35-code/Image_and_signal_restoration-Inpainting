% 
%
% MAIN_2_EX1: TEST IMAGE RESTORATION 
%             DENOISE + DEBLUR  by TV-L2 AND TIK-L2
%             INPAINTING        by TIk-L2
% 

% clear environment
close all; clear all; clc; format short;

addpath("GENERATE_DATA_2D\")

% initialize pseudo-random numbers generators
  rng('default'); rng(2);

% Set algorithms to be tested
  TEST_TV_L2_GRAD_DESC    = 1;
  TEST_TIK_L2_U_REST      = 0;
  TEST_TIK_L2_U_INP       = 0;

% SHOW WHAT?
  SHOW_RESTORED_IMAGES = 1;

  SAVE_RESTORED_IMAGES = 0;

% SET ALGORITHMS PARAMETERS for TV_L2 with Gradient Descent (GD) and TIKHONOV: 
  % set values of the regularization parameter to be tested
    lambdas_min     = 0.01;
    lambdas_max     = 0.3;
    lambdas_n       = 10;
    lambdas         = linspace(lambdas_min,lambdas_max,lambdas_n);

% allocate arrays where storing the obtained 
% relative errors, discrepancies

  % restored images: associated relative errors
    TIK_L2_U_REL_ERRS   = zeros(1,lambdas_n);

    %TIK_L2_U_RES_NORM   = zeros(1,lambdas_n);
    %TIK_L2_REG_NORM     = zeros(1,lambdas_n);

    GD_REL_ERRS         = zeros(1,lambdas_n);

  % Restored images: associated discrepancies
    TIK_L2_U_DISCREPS   = zeros(1,lambdas_n);
    GD_DISCREPS         = zeros(1,lambdas_n);

  % ISNR for the restored image
    ISNR = zeros(1,lambdas_n);


AEI_VIS_F = 2;

% -------------------------------------------------------------------------
% INITIALIZE
% -------------------------------------------------------------------------

% generate - and eventually visualize - the original 
% and corrupted (that is, blurred and noisy) image b
  GENERATE_DATA_2D;

% extract image dimensions
  [h,w] = size(b); % height and width, in pixels
  d     = h * w;   % total number of pixels

% compute relative error of b, error between original image and corrupted(
% blur + noise) one
  b_REL_ERR = compute_rel_err(xtr,b);

% -------------------------------------------------------------------------
% TIK-L2 RESTORATION %
%--------------------------------------------------------------------------

% Compute approximate solutions (that is, minimizers) of the
% unconstrained TIK-L2 variational model: 
%
% u*(mu) = argmin J(u;lambda),   J(u;lambda) = ||D u||_2^2 / 2 + lambda ||A u - b||_2^2 / 2  
%
if ( TEST_TIK_L2_U_REST == 1 ) & (exist('b_k') )
    fprintf('\n\n\nIMAGE RESTORATION by TIK-L2-U:\n');
    % for each given lambda, closed-form solution by solving the linear system:
    % (D^TD/lambda+A^TA)x = A^T b
    % we can solve the linear system by fast transforms...
    % we assume periodic boundary conditions -> solution by FFT
    % we pre-compute all terms in the system which does not depend on lambda
    K_DFT           = psf2otf(b_k,size(b));
    Dh_DFT          = psf2otf([1,-1],size(b));
    Dv_DFT          = psf2otf([1;-1],size(b));
    KT_DFT          = conj(K_DFT);
    DhT_DFT         = conj(Dh_DFT);
    DvT_DFT         = conj(Dv_DFT);
    DTD_DFT         = DhT_DFT .* Dh_DFT + DvT_DFT .* Dv_DFT;
    KTK_DFT         = KT_DFT .* K_DFT;
    KTb_DFT         = KT_DFT .* fft2(b);    
    for lambda_i = 1:lambdas_n
        %
        lambda = lambdas(lambda_i); 

        % Estimate of the solution reconstructed by TIK_L2
          TIK_L2_U_xest = real( ifft2( KTb_DFT ./ (KTK_DFT + lambda * DTD_DFT) ) );

        % Estimate of the Relative Error
          TIK_L2_U_ITRS(lambda_i)     = 1;
          TIK_L2_U_REL_ERRS(lambda_i) = compute_rel_err(xtr,TIK_L2_U_xest);

        % Estimate of the residual variance
          TIK_L2_U_RES                = real( ifft2( K_DFT .* fft2(TIK_L2_U_xest) - fft2(b) ) );
          TIK_L2_U_DISCREPS(lambda_i) = 255* norm( TIK_L2_U_RES(:) ) / sqrt(N);

        % Estimate of the residual and regularized norm
          %TIK_L2_U_RES_NORM(lambda_i) = 255 * norm( TIK_L2_U_RES(:), 2) ;
          %TIK_L2_REG_NORM(lambda_i)   = 255 * norm( TIK_L2_U_xest(:), 2) ;

        % ISNR
          ISNR(lambda_i) = compute_snr(TIK_L2_U_xest,xtr) - compute_snr(b,xtr);
              
        % show results
        fprintf('\nlambda(%02d) = %5.2f: rel-err = %6.4f (%6.4f), discr = %5.2f (%5.2f), its = %d\n', ...
            lambda_i,lambda,TIK_L2_U_REL_ERRS(lambda_i),b_REL_ERR,TIK_L2_U_DISCREPS(lambda_i),n_sigma,TIK_L2_U_ITRS(lambda_i));    
        
        if ( SAVE_RESTORED_IMAGES == 1 )
            imwrite(uint8(255*xtr),sprintf('../GENERATE_DATA_2D/output/%s_A_ORIG.png',im_name));
            imwrite(uint8(255*Kxtr),sprintf('../GENERATE_DATA_2D/output/%s_B_BLUR.png',im_name));
            imwrite(uint8(255*b),sprintf('../GENERATE_DATA_2D/output/%s_C_CORR.png',im_name));
            imwrite(uint8(255*TIK_L2_U_xest),sprintf('../GENERATE_DATA_2D/output/%s_D_REST_TIK-L2-U.png',im_name));
            imwrite(uint8(AEI_VIS_F*255*abs(xtr-TIK_L2_U_xest)),sprintf('../GENERATE_DATA_2D/output/%s_E_ERR_TIK-L2-U.png',im_name));
        end
    end

    [~,i] = min(TIK_L2_U_REL_ERRS);
    lambda_optimal = lambdas(i);
    fprintf('\nThe best value of lambda is %.4f\n',lambda_optimal);
    fprintf( 'The error in correspondence of %.4f is %.5f\n',lambda_optimal,min(TIK_L2_U_REL_ERRS) )

    if ( SHOW_RESTORED_IMAGES == 1 )
         figure(20)
         set(gcf,'Position',get(0,'ScreenSize'));
         subplot(2,2,1)
         imshow(uint8(255*xtr));
         title(sprintf('\\textbf{ORIGINAL} ( %d x %d )',h,w),'Interpreter','latex');
         subplot(2,2,2)
         imshow(uint8(255*b));
         title(sprintf('\\textbf{CORRUPTED}: $rel_{err} = %7.5f$',b_REL_ERR),'Interpreter','latex');
         subplot(2,2,3)
         imshow(uint8(255*TIK_L2_U_xest));
         title(sprintf('\\textbf{RESTORED by TIK-L2-U} ($\\lambda = %6.3f$): $rel_{err} = %7.5f$',lambda_optimal,TIK_L2_U_REL_ERRS(i)),'interpreter','latex');
         subplot(2,2,4)
         imshow(uint8( AEI_VIS_F*255*abs(xtr-TIK_L2_U_xest) )); 
         title('\textbf{ABSOLUTE ERROR IMAGE}','Interpreter','latex');
         pause(0.1);

         figure(22)
         %set(gcf,'Position',get(0,'ScreenSize'));
         hold on;
         plot(lambdas,ISNR,'b');plot(lambda_optimal,ISNR(i),'r*')
         xlabel('$\lambda$','Interpreter','latex');
         ylabel('ISNR($\lambda$) [dB]','Interpreter','latex');
         axis tight;grid on;
         title('TIK-L2-U: ACCURACY MEASURES versus $\lambda$','Interpreter','latex')
         legend('ISNR function','max value ISNR','location','northwest','Interpreter','latex')

    end
    
    % show multiple results
    if ( lambdas_n > 1 )
        figure(21)

        subplot(1,2,1)
        plot([lambdas(1),lambdas(end)],[b_REL_ERR,b_REL_ERR],'b'); hold all;
        plot(lambdas,TIK_L2_U_REL_ERRS,'r');
        xlabel('regularization parameter $\lambda$','Interpreter','latex');
        ylabel('$rel_{err}(\lambda)$','Interpreter','latex');
        legend('corrupted image $rel_{err}$','TIK-L2 image $rel_{err}$','interpreter','latex');
        title('TIK-L2-U:  $rel_{err}(\lambda)$','Interpreter','latex'); axis tight;

        subplot(1,2,2)
        plot([lambdas(1),lambdas(end)],[n_sigma,n_sigma],'b'); hold all;
        plot(lambdas,TIK_L2_U_DISCREPS,'r');
        xlabel('regularization parameter $\lambda$','Interpreter','latex'); 
        ylabel('discr($\lambda$)','Interpreter','latex');
        legend('Noise variance','Residual variance');
        title('TIK-L2-U:  discr($\lambda$)','Interpreter','latex');axis tight;
    end
end

%--------------------------------------------------------------------------
% TV-L2 RESTORATION %
%--------------------------------------------------------------------------
if ( TEST_TV_L2_GRAD_DESC == 1 ) & (exist('b_k'))
    % set parameters
    eps = 0.001;
    tau = eps/6;
    itrs_max = 2000;
    ch_th = 1e-5;
    if ( lambdas_n == 1 )
        debug_cw = 1;
    else
        debug_cw = 0;
    end
    % set iterations initial guess
    x0 = b;
    %x0 = mean(b(:)) * ones(size(b));
    fprintf('\n\n\nIMAGE RESTORATION by TV-L2 and GRADIENT DESCENT:\n');
    for lambda_i = 1 : lambdas_n
        %
        lambda = lambdas(lambda_i);
        
        %%TODO implement
        [TV_L2_U_xest,itrs,r] = TV_L2_U_GD(b,b_k,lambda,eps,x0,tau,itrs_max,ch_th,debug_cw,xtr);    
        %xest  = xtr;  %TOBE CHANGED
        %itrs  = 1;    %TOBE CHANGED
         
        GD_REL_ERRS(lambda_i) = compute_rel_err(xtr,TV_L2_U_xest);
        GD_RES                = real( ifft2( K_DFT .* fft2(TV_L2_U_xest) - fft2(b) ) );
        GD_DISCREPS(lambda_i) = 255*norm( GD_RES(:) ) / sqrt(N);

        ISNR(lambda_i) = compute_snr(TV_L2_U_xest,xtr) - compute_snr(b,xtr);
        
        fprintf('\nlambda(%02d) = %5.4f:  rel-err = %6.6f (%6.4f), discr = %5.2f (%5.2f), itrs = %d',lambda_i,lambda,GD_REL_ERRS(lambda_i),b_REL_ERR,GD_DISCREPS(lambda_i),n_sigma,itrs);
     
    end

    [~,i] = min(GD_REL_ERRS);
    lambda_optimal = lambdas(i);

    if ( SHOW_RESTORED_IMAGES == 1 )
            figure(20)
            set(gcf,'Position',get(0,'ScreenSize'));
            subplot(1,3,1)
            imshow(uint8(255*xtr));
            title(sprintf('\\textbf{ORIGINAL} (%d x %d)',h,w),'interpreter','latex');
            subplot(1,3,2)
            imshow(uint8(255*b));
            title(sprintf('\\textbf{CORRUPTED}: $rel_{err}$ = %7.5f',b_REL_ERR),'Interpreter','latex');
            subplot(1,3,3)
            imshow(uint8(255*TV_L2_U_xest));
            title(sprintf('\\textbf{RESTORED by TV-L2 GD} $\\lambda = %6.5f$: $rel_{err}$ = %7.5f',lambda_optimal,GD_REL_ERRS(i)),'Interpreter','latex');
            pause(1);

            figure(22)
            hold on;
            plot(lambdas,ISNR,'b');plot(lambda_optimal,ISNR(i),'r*')
            xlabel('$\lambda$','Interpreter','latex');
            ylabel('ISNR($\lambda$) [dB]','Interpreter','latex');
            axis tight;grid on;
            title('TV-L2-U: ACCURACY MEASURES versus $\lambda$','Interpreter','latex')
            legend('ISNR function','max value ISNR','Interpreter','latex')
    end
    
    if ( lambdas_n > 1 )
        figure(21)

        subplot(1,2,1)
        plot([lambdas(1),lambdas(end)],[b_REL_ERR,b_REL_ERR],'b'); hold all;
        plot(lambdas,GD_REL_ERRS,'r');
        xlabel('regularization parameter $\lambda$','interpreter','latex');
        ylabel('$rel_{err}(\lambda)$','interpreter','latex');
        legend('corrupted image $rel_{err}$','interpreter','latex');
        title('TV-L2-U: $rel_{err}(\lambda)$','Interpreter','latex'); axis tight;
        axis tight;

        subplot(1,2,2)
        plot([lambdas(1),lambdas(end)],[n_sigma,n_sigma],'b'); hold all;
        plot(lambdas,GD_DISCREPS,'r');
        xlabel('regularization parameter $\lambda$','Interpreter','latex');
        ylabel('discr($\lambda$)','Interpreter','latex');
        legend('noise stdv','interpreter','latex');
        title('TV-L2-U:  discr($\lambda$)','Interpreter','latex');axis tight;
        axis tight;
    end
end

%--------------------------------------------------------------------------
% TIK-L2 INPAINTING %
%--------------------------------------------------------------------------

% Compute approximate solutions (that is, minimizers) of the
% unconstrained TIK-L2 variational model: 
%
% u*(mu) = argmin J(u;lambda),   J(u;lambda) = ||D u||_2^2 / 2 + lambda ||M u - b||_2^2 / 2  
% 
if ( TEST_TIK_L2_U_INP == 1 )
    error = [];

    Dh_1d = spdiags([-ones(w,1),ones(w,1)], 0:1, w ,w); Dh_1d(end,1) = 1;
    Dv_1d = spdiags([-ones(h,1),ones(h,1)], 0:1, h ,h); Dv_1d(end,1) = 1; 

    Ih = speye(h); Iv = speye(w);

    Dh = kron(Dh_1d,Ih); Dv = kron(Iv,Dv_1d);

    DTD = Dh' * Dh + Dv' * Dv;

    for lambda_i = 1 : lambdas_n
        %
        lambda = lambdas(lambda_i);
      
        % TODO solve for xest with  M inpainting matrix
        N = h*w;
        S = spdiags(M(:),0,N,N);

        TIK_L2_U_INP = reshape( (S + lambda * DTD) \ ( S*b(:)), h, w) ;

        rel_err_f(lambda_i) = compute_rel_err(xtr,TIK_L2_U_INP);

        error = [error;lambda, rel_err_f(lambda_i)];

        ISNR(lambda_i) = compute_snr(TIK_L2_U_INP,xtr) - compute_snr(b,xtr);

        TIK_L2_INP_RES                =  M .* TIK_L2_U_INP - b;
        INP_DISCREPS(lambda_i) = 255*norm( TIK_L2_INP_RES(:) ) / sqrt(N);

        fprintf('\nlambda(%02d) = %5.4f: rel-err = %6.5f (%6.5f)\n', ...
             lambda_i,lambda,rel_err_f(lambda_i),b_REL_ERR);   

    end

    [~,i] = min(error(:,2));
    lambda_optimal = lambdas(i);

    if ( SHOW_RESTORED_IMAGES == 1 )

            figure(20)
            set(gcf,'Position',get(0,'ScreenSize'));
            subplot(1,3,1)
            imshow(uint8(255*xtr));
            title(sprintf('\\textbf{ORIGINAL} (%d x %d)',h,w),'Interpreter','latex');
            subplot(1,3,2)
            imshow(uint8(255*b));
            title(sprintf('\\textbf{CORRUPTED}: $rel_{err}$ = %7.5f',b_REL_ERR),'Interpreter','latex');
            subplot(1,3,3)
            imshow(uint8(255*TIK_L2_U_INP));
            title(sprintf('\\textbf{RESTORED by TIK-L2} ($\\lambda$ = %6.5f): $rel_{err}$ = %7.5f',lambda_optimal,rel_err_f(i)),'Interpreter','latex');
            pause(1);

            figure(22)
            hold on;
            plot(lambdas,ISNR,'b');plot(lambda_optimal,ISNR(i),'r*')
            xlabel('$\lambda$','Interpreter','latex');
            ylabel('ISNR($\lambda$) [dB]','Interpreter','latex');
            axis tight;grid on;
            title('TIK-L2-INPAINTING: ACCURACY MEASURES versus $\lambda$','Interpreter','latex')
            legend('ISNR function','max value ISNR','location','northwest','Interpreter','latex')
    end

    if ( (n_sigma ~= 0) &&  (lambdas_n > 1) )
        figure(21)

        subplot(1,2,1)
        plot([lambdas(1),lambdas(end)],[b_REL_ERR,b_REL_ERR],'b'); hold all;
        plot(lambdas,rel_err_f,'r');
        xlabel('regularization parameter $\lambda$','interpreter','latex');
        ylabel('$rel_{err}(\lambda)$','interpreter','latex');
        legend('corrupted image $rel_{err}$','interpreter','latex');
        title('TIK-L2-INPAINTING: $rel_{err}(\lambda)$','Interpreter','latex'); axis tight;
        axis tight;

        subplot(1,2,2)
        plot([lambdas(1),lambdas(end)],[n_sigma,n_sigma],'b'); hold all;
        plot(lambdas,INP_DISCREPS,'r'); % TO BE CHANGEd
        xlabel('regularization parameter $\lambda$','Interpreter','latex');
        ylabel('discr($\lambda$)','Interpreter','latex');
        legend('noise stdv','interpreter','latex');
        title('TV-L2-INPAITING:  discr($\lambda$)','Interpreter','latex');axis tight;
        axis tight;
    end

end