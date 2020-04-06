function initPClockFromPEph(obj)
 % Use the clock data from the .sp3 to populate the precise clock field
 
 
 Clck.Cepochs = obj.PEph.epochs';
 Clck.Cclk    = obj.PEph.clock_bias';
 Clck.Cclk_sig = nan(size(Clck.Cclk));
 
 Clck.PRNs = obj.PEph.PRN;
 Clck.constInds = obj.PEph.constellation;


obj.PClock = Clck;

end