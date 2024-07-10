classdef PIBIout3quBits < NVProblem
    properties
       forceReal = true;
    end
    methods
        function X = sampleOperators(self)
            U1 = qdimsum.Random.unitary(3);
            U2 = qdimsum.Random.unitary(3);
            P = diag([1 1 0]);  % Change to P = diag([1 1 1]) to obtain the bound for qutrits -0.5
            X = {P, U1(1,:)'*U1(1,:), U1(2,:)'*U1(2,:), U2(1,:)'*U2(1,:), U2(2,:)'*U2(2,:)};
        end
        function rho = sampleStateKraus(self)
            rho = qdimsum.Random.pureNormalizedDensityMatrix(3);
        end
        function obj = computeObjective(self, X, rho)
            bellOp1 = X{1}*X{2}*X{1} + X{1}*X{5}*X{1} - X{1}*(X{2}-X{5})*X{1}*(X{2}-X{5})*X{1} + X{1}*X{3}*X{1} + X{1}*X{4}*X{1} - X{1}*(X{3}-X{4})*X{1}*(X{3}-X{4})*X{1};
            obj = -real(trace(bellOp1*rho));
        end   
        function generators = symmetryGroupGenerators(self)
              generators = [1 2 3 4 5; 1 3 2 5 4; 1 4 5 2 3]; 
         end        
    end
end

