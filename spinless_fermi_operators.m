function siteops = spinless_fermi_operators()
%SPINLESS_FERMI_OPERATORS Summary of this function goes here
%   Detailed explanation goes here
  siteops = struct();
  siteops.cdag = [0 0; 1 0];
  siteops.c = (siteops.cdag)';
  siteops.id = eye(2);
  siteops.ph = [1 0; 0 -1];
  siteops.cdag_ph = siteops.cdag*siteops.ph;
  siteops.c_ph = siteops.c*siteops.ph;
  siteops.n = siteops.cdag*siteops.c;
end

