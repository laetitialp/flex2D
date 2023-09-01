function [BCDOF,BCVAL] =  bcapply(nodesx,nodesy,nodeso,valx,valy,valo,BCDOF,BCVAL)
thetaxdof = (nodesx-1)*3+1;
thetaxval = ones(size(thetaxdof)).*valx;

thetaydof = (nodesy-1)*3+2;
thetayval = ones(size(thetaydof)).*valy;

omegadof  = (nodeso-1)*3+3;
omegaval = ones(size(omegadof)).*valo;

BCDOF = [BCDOF,thetaxdof,thetaydof,omegadof];
BCVAL = [BCVAL,thetaxval,thetayval,omegaval];
