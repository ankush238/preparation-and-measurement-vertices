d=2; % d is the dimension of the communicating system
ops=sdpsettings('solver','mosek','verbose',0);
D1=2/3; D2=1; % D1, D2 are the distinguishability of first and second sender respectively

A1v=sdpvar(d,d,'hermitian','complex'); sig=sdpvar(d,d,'hermitian','complex'); % sig is the auxiliary variable to bound the distinguishability of the first sender
A2v=sdpvar(d,d,'hermitian','complex');
A3v=sdpvar(d,d,'hermitian','complex');
conA=[A1v>=0;A2v>=0;A3v>=0;trace(A1v)==1;trace(A2v)==1;trace(A3v)==1;sig>=(1/3)*A1v;sig>=(1/3)*A2v;sig>=(1/3)*A3v;trace(sig)<=D1];

B1v=sdpvar(d,d,'hermitian','complex'); sgi=sdpvar(d,d,'hermitian','complex'); % sgi is the auxiliary variable to bound the distinguishability of the second sender
B2v=sdpvar(d,d,'hermitian','complex');
conB=[B1v>=0;B2v>=0;trace(B1v)==1;trace(B2v)==1;sgi>=(1/2)*B1v;sgi>=(1/2)*B2v;trace(sgi)<=D2];

M1v=sdpvar(d^2,d^2,'hermitian','complex');
M2v=sdpvar(d^2,d^2,'hermitian','complex');
conM=[M1v>=0;M2v>=0;M1v+M2v==eye(d^2)];

A1 = RandomDensityMatrix(d); B1 = RandomDensityMatrix(d);
A2 = RandomDensityMatrix(d); B2 = RandomDensityMatrix(d);
A3 = RandomDensityMatrix(d);

Sq=zeros(1,3);
for r=1:3
    %  p(2|1, 1) − 3p(2|2, 1) + p(2|3, 1) − p(2|1, 2) + p(2|2, 2) + p(2|3, 2) ⩽ 6D1 + 2D2 − 3
    for iter = 1:25
        obj_state = trace(0*M1v * (-kron(A1, B1) + kron(A2, B1) + 0*kron(A3, B1) -kron(A1, B2) - kron(A2, B2) + kron(A3, B2)) +...
                            M2v * (kron(A1, B1) -3* kron(A2, B1) + kron(A3, B1)-kron(A1, B2) + kron(A2, B2)+ kron(A3, B2)));
        game_M = optimize(conM, -real(obj_state), ops);
        M1 = value(M1v);
        M2 = value(M2v);

        obj_A = trace(0*M1 * (-kron(A1v, B1) + kron(A2v, B1) + 0*kron(A3v, B1) -kron(A1v, B2) - kron(A2v, B2) + kron(A3v, B2)) +...
                        M2 * (kron(A1v, B1) -3* kron(A2v, B1) + kron(A3v, B1)-kron(A1v, B2) + kron(A2v, B2)+ kron(A3v, B2)));
        game_A = optimize(conA, -real(obj_A), ops);
        A1 = value(A1v);
        A2 = value(A2v);
        A3 = value(A3v);

        obj_B = trace(0*M1 * (-kron(A1, B1v) + kron(A2, B1v) + 0*kron(A3, B1v) -kron(A1, B2v) - kron(A2, B2v) + kron(A3, B2v)) +...
                       M2 * (kron(A1, B1v) -3* kron(A2, B1v) + kron(A3, B1v)-kron(A1, B2v) + kron(A2, B2v)+ kron(A3, B2v)));
        game_B = optimize(conB, -real(obj_B), ops);
        B1 = value(B1v);
        B2 = value(B2v);
    end
     Sq(1,r)=value(obj_B);
end

bestSq=max((Sq));
bestSq1=real(bestSq);
disp(bestSq1);
