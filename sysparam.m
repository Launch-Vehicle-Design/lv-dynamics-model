function param = sysparam(filename)

    param.m0 = 2267.962;
    param.L = 240*0.0254;
    param.dia = 24*0.0254;
    
    % 1st stage
    param.ms1 = 120
    param.mp1 = 
    
    % 2nd stage
    param.ms2 = 
    param.mp2 = 

    if exist("filename","var")
        load(filename)
    end