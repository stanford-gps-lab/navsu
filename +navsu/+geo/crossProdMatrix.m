function Omega = crossProdMatrix(w)

% creating the cross product matrix (w x)

Omega = [  0   -w(3)  w(2);
          w(3)   0   -w(1);
         -w(2)  w(1)   0  ];


end