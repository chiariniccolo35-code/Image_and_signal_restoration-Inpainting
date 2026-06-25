%clear all; clc; close all;

rng('default'); rng(1);

addpath('GENERATE_DATA_1D/input')

%-----------------------------------------------------------------------
% RESTORATION (NOISE + BLUR)
%-----------------------------------------------------------------------
  % Forward degradation model: given xtr (signal) compute b = A*xtr+err
  % OUTPUT:
    %   xtr         original uncorrupted signal
    %   Axtr A*xtr  blurred signal
    %   b           blurred and noisy signal

%-----------------------------------------------------------------------
% INPAINTNG ( directly made inside the MAIN_1_EX2 )
%-----------------------------------------------------------------------
  % Reconstruction of the missing samples of the signal
  % OUTPUT:
    % xtr           original uncorrupted signal
    % b             known samples of the original signal xtr
  % OBS => For the inpainting case b is equal to Axtr

SHOW_CORRUPTIONS_DENOISE     = 0;
SHOW_CORRUPTIONS_RESTORATION = 1;

% -------------------------------------------------------------------------
% SELECT THE UNCORRUPTED SIGNAL TO BE TESTED

s_id = 7; % index of the uncorrupted signal to be tested (see the switch-case below)
switch s_id    
    case 0 % 
        s_name     = 'SIGNAL_0';
        f           = 1;
        xtr         = 20 * ones(200,1) / f;
    case 1 % 
        s_file     = 'input/signal_1.mat';
        s_name     = 'SIGNAL_1';
        xtr         = load(s_file);
        xtr         = xtr.xtr;
    case 2 % 
        s_file     = 'input/signal_2.mat';
        s_name     = 'SIGNAL_2';
        xtr         = load(s_file);
        xtr         = xtr.xtr;
    case 3 % 
        s_file     = 'input/signal_3.mat';
        s_name     = 'SIGNAL_3';
        xtr         = load(s_file);
        xtr         = xtr.xtr;
    case 4 % 
        s_file     = 'input/signal_4.mat';
        s_name     = 'SIGNAL_4';
        xtr         = load(s_file);
        xtr         = xtr.xtr;
    case 5 % 
        s_file     = 'input/signal_5.mat';
        s_name     = 'SIGNAL_5';
        xtr         = load(s_file);
        xtr         = xtr.xtr;
    case 6 % 
        s_file     = 'input/signal_6.mat';
        s_name     = 'SIGNAL_6';
        xtr         = load(s_file);
        xtr         = xtr.xtr;
    case 7 % 
        s_name     = 'SIGNAL_7';
        xtr=zeros(200,1);
        for k=1:200
            if (k<50)
                xtr(k)=0;
            elseif (k<100)
                xtr(k)=1;
            elseif (k<150)
                xtr(k)=4;
            else
                xtr(k)=0;
            end
        end
    case 8 % from Regularization toolbox 
        s_name     = 'PHILLIPS';
        xtr        = zeros(128,1); %The order n must be a multiple of 4 
        [A,Axtr,xtr]  = phillips(128);
        if SHOW_CORRUPTIONS_DENOISE, Axtr = xtr; end
    case 9 % from Regularization toolbox 
        s_name     = 'BAART';
        xtr        = zeros(256,1); 
        [A,Axtr,xtr]  = baart(256);
        if SHOW_CORRUPTIONS_DENOISE, Axtr = xtr; end
    case 10 % from Regularization toolbox 
        s_name     = 'SHAW';
        xtr        = zeros(256,1); 
        [A,Axtr,xtr]  = shaw(256);
        if SHOW_CORRUPTIONS_DENOISE, Axtr = xtr; end
    case 11 % smooth noisy data (for denoising exercise) 
        % do not add extra noise or blur
        % A noisy ECG signal b (ECG waveform generator ECGSYN)
        s_file     = 'input/ECG_data.mat';
        s_name     = 'SIGNAL_SMOOTHING';
        load(s_file,'-mat');
        xtr  = b;
        Axtr = b;
    case 12 % smooth noisy data (Fourier domain) 
        % do not add extra noise or blur
        % A noisy wav file signal b 
        % Load speech waveform data
        s_file     = 'input/sp1.wav';
        [sp1, fs] = wavread(s_file);
        fprintf('Sampling rate: %d samples/second \n', fs)
        M = 500;                        % M : length of signal
        xtr = sp1(5500+(1:M));            % s : signal (without noise)
        s_name     = 'SIGNAL_Speech_waveform';
        % Make noisy signal by adding white Gaussian noise
        w = 0.1 * randn(M,1);      % w : zero-mean Gaussian noise
        b = xtr + w;                 % y : noisy speech signal
        Axtr = b;
end

% extract signal dimension
d = numel(xtr);

% -------------------------------------------------------------------------
% SET THE BLUR PARAMETERS, THEN GENERATE AND STORE THE BLUR PSF (KERNEL) b_k, 
% FINALLY BLUR THE ORIGINAL IMAGE xtr --> A xtr
% For (signal_id <  8),  CONSTRUCT the BLUR matrix A
% For (signal_id >= 8),  matrix A, if needed, is previously defined.
% -------------------------------------------------------------------------
if s_id < 8

    % set the blur parameters then generate and store the blur kernel, b_k
    if(SHOW_CORRUPTIONS_DENOISE == 1)
        %
        fprintf('\nPRESS ANY KEY TO PROCEDE\n');
        pause    
        Axtr = xtr;
    else
        b_k_type = 2; % index of the blur type (see the switch-case below)
        switch b_k_type        
               case 1 % average
                   b_r             = 5; % radius (pixels)
                   b_bc_type       = 0; % boundary conditions type (1->Dirilichet; 2->Periodic; 3->Reflective)
                   [b_k,b_k_c]     = generate_b_k_1D(b_k_type,b_r); % kernel
                   b_type_descr    = 'AVER';
               case 2 % Gaussian
                    b_r             = 5; % radius (pixels) ...the band is 2 * b_r + 1
                    b_s             = 0.5; % standard deviation
                    b_bc_type       = 3; % boundary conditions type (1->Dirilichet; 2->Periodic; 3->Reflective)
                    [b_k,b_k_c]     = generate_b_k_1D(b_k_type,[b_r,b_s]); % kernel
                    b_type_descr    = 'GAUSS';
               case 3 % motion
                    b_l             = 10; % length (pixels) 
                    b_c             = 2; % center
                    b_bc_type       = 0; % boundary conditions type (1->Dirilichet; 2->Periodic; 3->Reflective)
                    [b_k,b_k_c]     = generate_b_k_1D(b_k_type,[b_l,b_c]); % kernel
                    b_type_descr    = 'MOTION';
       otherwise disp('ERROR: blur does not exist!')
       end
 
       % generate the blur (convolution) matrix (according to the blur kernel)
         A = generate_A_from_b_k_1D (b_k,b_k_c,b_bc_type,d);
       % compute the blurred image
         fprintf('\nPRESS ANY KEY TO PROCEDE\n');
         pause
         Axtr  = A * xtr;
    end
end

% -------------------------------------------------------------------------
% SET THE NOISE PARAMETERS, THEN ADD NOISE TO THE BLURRED IMAGE A xtr --> b
% -------------------------------------------------------------------------
n_type_id = 1; % index of the noise type (see the switch-case below)
switch n_type_id
    case 0 % no noise
        n_mean          = 0.0;  % mean 
        n_sigma         = 0.0;  % stdv 
        n_type_descr    = 'NO NOISE';
        b               = Axtr;
    case 1 % additive white Gaussian noise (AWGN)
        n_mean          = 0.0;    % mean 
        n_sigma         = 0.05;   % stdv 
        n_type_descr    = 'AWGN';
        n_realiz        = n_mean + n_sigma * randn(size(Axtr));
        b               = Axtr + n_realiz;
        norm2_error     = norm(n_realiz)^2;  
        % n_sigma^2*d = squared norm of the noise (norm2_error)
    case 2 % additive white uniform noise (AWUN) 
           % (given n_mean and n_sigma -> uniform pdf in [n_mean-sqrt(3)n_sigma,n_mean+sqrt(3)n_sigma])
        n_type_descr    = 'AWUN';
        n_mean          = 0.0;  % mean 
        n_sigma         = 0.05;   % stdv 
        n_hs            = n_sigma * sqrt(3); % half-size of the support        
        n_realiz        = n_mean + n_hs * ( 2 * (rand(size(Axtr)) - 0.5) );
        b               = Axtr + n_realiz;
    case 3 % additive white Laplace noise (AWLN)
        n_type_descr    = 'AWLN';
        n_mean          = 0.0;  % mean 
        n_sigma         = 20;   % stdv 
        n_sp            = n_sigma / sqrt(2); % scale parameter
        U               = rand(size(Axtr)) - 0.5;
        n_realiz        = n_mean - n_sp * ( sign(U) .* log(1 - 2 * abs(U)) );
        b               = Axtr + n_realiz;
    case 4 % multiplicative white Gaussian noise (MWGN)
        n_mean          = 0.0;  % mean 
        n_sigma         = 0.5;  % stdv 
        n_type_descr    = 'MWGN';
        n_realiz        = 1 + ( (n_sigma) * randn(size(Axtr)) );
        b               = Axtr .* n_realiz;    
    case 7 % impulsive: salt & pepper
        n_type_descr    = 'ISPN';
        n_sigma = 11;
        p   = 0.1;
        ph  = p/2;
        P   = rand(size(Axtr));
        b   = Axtr;
        b(P < ph)       = 0;
        b(P > (1-ph))   = 1;
    case 8 % impulsive: random-valued
        n_type_descr    = 'IRVN';
        p   = 0.1;
        P   = rand(size(Axtr));
        P1  = rand(size(Axtr));
        b   = Axtr;
        b(P < p) = P1(P < p);
end

if (SHOW_CORRUPTIONS_RESTORATION == 1)
    % Plot signals
    picturewidth = 20;
    hw_ratio     = 0.7;     %Set for the height

    fig_1 = figure(101);

    Axtr_color = 'm';
    b_color = 'b';
    xtr_color = 'r';

    x_min = min( [ min(Axtr) min(b) min(xtr) ] );
    x_max = max( [ max(Axtr) max(b) max(xtr) ] );
    x_range = x_max - x_min;
    x_min = x_min - 0.1 * x_range - 0.1;
    x_max = x_max + 0.1 * x_range + 0.1;

    hold on;
    plot(1:d,xtr,xtr_color,'LineWidth',1.5);
    plot(1:d,Axtr,Axtr_color,'LineWidth',1.5);
    plot(1:d,b,b_color,'LineWidth',1.5);
    axis ([1 d x_min x_max]);
    xlabel('i'); ylabel('x(i)');
    title('\textbf{ x (original), Ax (blurred) and b (blur+noise) }','Interpreter','latex');
    legend('x','Ax','b','Location','northwest');grid off;

    set(findall(fig_1,'-property','FontSize'),'FontSize',21);
    set(findall(fig_1,'-property','Box'),'Box','off');
    set(findall(fig_1,'-property','Interpreter'),'interpreter','latex');
    set(findall(fig_1,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex');
    set(fig_1,'units','centimeters','Position',[4 4 picturewidth hw_ratio * picturewidth])
    pos = get(fig_1,'Position');
    set( fig_1,'PaperPositionMode','Auto','PaperUnits','centimeters','papersize',[pos(3),pos(4)] )
    print(fig_1,'pdf_figure','-dpdf','-vector','-bestfit ');
end

if (SHOW_CORRUPTIONS_DENOISE == 1)
    %
    picturewidth = 20;
    hw_ratio     = 0.7;     %Set for the height

    fig = figure(101);

    b_color = 'b';
    xtr_color = 'r';

    x_min = min( [ min(b) min(b) min(b) ] );
    x_max = max( [ max(b) max(b) max(b) ] );
    x_range = x_max - x_min;
    x_min = x_min - 0.1 * x_range - 0.1;
    x_max = x_max + 0.1 * x_range + 0.1;

    hold on;
    plot(1:d,xtr,xtr_color,'LineWidth',1.5);plot(1:d,b,b_color,'LineWidth',1.5);
    axis ([1 d x_min x_max]);
    xlabel('i'); ylabel('x(i)');
    title('\textbf{x (original signal)}, \textbf{b (noisy signal)}','Interpreter','latex');
    legend('x','b');grid off;

    set(findall(fig,'-property','FontSize'),'FontSize',21);
    set(findall(fig,'-property','Box'),'Box','off');
    set(findall(fig,'-property','Interpreter'),'interpreter','latex');
    set(findall(fig,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex');
    set(fig,'units','centimeters','Position',[3 3 picturewidth hw_ratio * picturewidth])
    pos = get(fig,'Position');
    set( fig,'PaperPositionMode','Auto','PaperUnits','centimeters','papersize',[pos(3),pos(4)] )
    print(fig,'pdf_figure','-dpdf','-vector','-bestfit');

end

