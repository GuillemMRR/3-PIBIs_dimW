monos = {'npa', 4}  
qubits_NV = -nvOptimize(PIBIout3quBits, monos, 'irreps' , NVSettings)

% We confirm the variational bound, -0.25. Note that the irreps are of
% dimenson 1, is a linear program.

