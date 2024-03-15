function dA = model_rhs_treatment(t,A,p,ct,cR,rho1,rho2,effect)
    % Function for calculating the right-hand sides of the ODEs of the
    % model. The effects of RUX are taken into account by the ct, cR and rho 
    % variables.
    % Input:
    %   t: Time
    %   A: Vector containing the variables of the Cancitis model
    %       A(1): x0, haematopoietic stem cells (HSC)
    %       A(2): x1, haematopoietic progenitor cells (HPC)
    %       A(3): x2, mature blood cells (MBC)
    %       A(4): y0, malignant haematopoietic stem cells (mHSC)
    %       A(5): y1, malignant haematopoietic progenitor cells (mHPC)
    %       A(6): y2, malignant mature blood cells (mMBC)
    %       A(7): a, cellular debris
    %       A(8): s, cytokine signal
    %   p: Struct containing the parameters of the model
    %       p.alphax0: Proliferation rate for HSCs
    %       p.px0: Fraction of self-renewal for HSCs
    %       p.alphay0: Proliferation rate for mHSCs
    %       p.py0: Fraction of self-renewal for mHSCs
    %       p.Ax0: Multiplication factor for differentiation of HSC to HPC
    %       p.Ay0: Multiplication factor for differentiation of mHSC to mHPC
    %       p.alphax1: Proliferation rate for mHSCs
    %       p.px1: Probability of self-renewal for mHSCs
    %       p.dx1: Death rate for mHSCs
    %       p.alphay1: Proliferation rate for mHSCs
    %       p.py1: Probability of self-renewal for mHSCs
    %       p.dy1: Death rate for mHSCs
    %       p.Ax1: Multiplication factor for differentiation of HPC to MBC
    %       p.Ay1: Multiplication factor for differentiation of mHPC to mMBC
    %       p.ea: Elimination rate for cellular debris
    %       p.rs: Rate for stimulation of immune system due to dead cells
    %       p.es: Elimination rate for cytokines (stimulation of immune
    %       system)
    %       p.I: (Time-dependent) exogene stimulation of the immune system,
    %       here taken to be piece-wise constant and thus a parameter
    %       p.cxx: Inhibitory strength of the signaling bone marrow niche
    %       feedback from HSC to HSC
    %       p.cxy: Inhibitory strength of the signaling bone marrow niche
    %       feedback from HSC to mHSC
    %       p.cyx: Inhibitory strength of the signaling bone marrow niche
    %       feedback from mHSC to HSC
    %       p.cyy: Inhibitory strength of the signaling bone marrow niche
    %       feedback from mHSC to mHSC
    %       p.sx0: Parameter for controlling Michaelis-Menten factor on the
    %       self-renewal terms of the HSCs
    %       p.sy0: Parameter for controlling Michaelis-Menten factor on the
    %       self-renewal terms of the mHSCs
    %   ct: Vector with times for the drug concentrations
    %   cR: Vector with RUX concentrations corresponding to the times in ct
    %   rho1: Parameter that determines the strength of the patient's 
    %   response to RUX treatment in terms of the first parameter (sy0).
    %   rho2: Parameter that determines the strength of the patient's 
    %   response to RUX treatment in terms of the second parameter (dy1).
    %   effect: String determining which effect of the drug that is used
    %   for fitting. See the code for the available options.
    %       
    % Output:
    %   dA: Vector with right-hand sides of the ODEs in the Cancitis model

    % Storage
    dA = zeros(8,1);

    % Pick out each variable
    x0 = A(1);
    x1 = A(2);
    x2 = A(3);
    y0 = A(4);
    y1 = A(5);
    y2 = A(6);
    a = A(7);
    s = A(8);

    % Pick out each parameter
    alphax0 = p.alphax0;
    px0 = p.px0;
    alphay0 = p.alphay0;
    py0 = p.py0;
    Ax0 = p.Ax0;
    Ay0 = p.Ay0;
    alphax1 = p.alphax1;
    px1 = p.px1;
    dx1 = p.dx1;
    alphay1 = p.alphay1;
    py1 = p.py1;
    dy1 = p.dy1;
    Ax1 = p.Ax1;
    Ay1 = p.Ay1;
    dx2 = p.dx2;
    dy2 = p.dy2;
    ea = p.ea;
    rs = p.rs;
    es = p.es;
    I = p.I;
    cxx = p.cxx;
    cxy = p.cxy;
    cyx = p.cyx;
    cyy = p.cyy;
    sx0 = p.sx0;
    sy0 = p.sy0;

    % Find time in ct closest to current time
    [~,ID] = min(abs(ct-t));

    % Extract corresponding drug doses
    cnowR = cR(ID);

    % Add drug effect to the relevant parameter
    if strcmp(effect,'sy0dy1')
        sy0 = (1+rho1*cnowR)*sy0;
        dy1 = (1+rho2*cnowR)*dy1;
    elseif strcmp(effect,'sy0')
        sy0 = (1+rho1*cnowR)*sy0;
    elseif strcmp(effect,'dy1')
        dy1 = (1+rho2*cnowR)*dy1;
    end

    % Calculate phi-functions
    phix = 1/(1+cxx*x0+cxy*y0);
    phiy = 1/(1+cyx*x0+cyy*y0);

    % Rhs of ODEs
    x0dot = alphax0*(2*px0*phix*s/(sx0+s)-1)*x0;
    x1dot = alphax1*(2*px1-1)*x1+2*Ax0*alphax0*(1-px0*phix*s/(sx0+s))*x0-dx1*x1; 
    x2dot = 2*Ax1*alphax1*(1-px1)*x1-dx2*x2;
    y0dot = alphay0*(2*py0*phiy*s/(sy0+s)-1)*y0;
    y1dot = alphay1*(2*py1-1)*y1+2*Ay0*alphay0*(1-py0*phiy*s/(sy0+s))*y0-dy1*y1; 
    y2dot = 2*Ay1*alphay1*(1-py1)*y1-dy2*y2;
    adot = dx1*x1+dy1*y1+dx2*x2+dy2*y2-ea*a*s;
    sdot = rs*a-es*s+I;

    % Collect rhs in vector
    dA(1) = x0dot;
    dA(2) = x1dot;
    dA(3) = x2dot;
    dA(4) = y0dot;
    dA(5) = y1dot;
    dA(6) = y2dot;
    dA(7) = adot;
    dA(8) = sdot;
end

