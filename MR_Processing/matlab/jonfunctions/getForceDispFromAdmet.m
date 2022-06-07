function P = getForceDispFromAdmet(admetfile,disps)
    fstruct = load(admetfile);
    admetdisps = fstruct.loaddisp.Positionmm;
    admetloads = fstruct.loaddisp.LoadN;
    
    P = interpn(smooth(admetdisps),smooth(admetloads),disps);
end