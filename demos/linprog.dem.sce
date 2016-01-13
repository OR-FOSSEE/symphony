mode(1)
//
// Demo of linprog.sci
//

Optimal problems
//Linear program, linear inequality constraints
f=[-1,-1/3];
A=[1,1;1,1/4;1,-1;-1/4,-1;-1,-1;-1,1];
b=[2,1,2,1,-1,2]
[xopt,fopt,exitflag,output,lambda]=linprog(f, A, b);
//Linear program with Linear Inequalities and Equalities
f=[-1,-1/3];
A=[1,1;1,1/4;1,-1;-1/4,-1;-1,-1;-1,1];
b=[2,1,2,1,-1,2]
Aeq=[1,1/4]
beq=[1/2]
[xopt,fopt,exitflag,output,lambda]=linprog(f, A, b, Aeq, beq);
//Linear program with all constraint types
f=[-1,-1/3];
A=[1,1;1,1/4;1,-1;-1/4,-1;-1,-1;-1,1];
b=[2,1,2,1,-1,2]
Aeq=[1,1/4]
beq=[1/2]
lb=[-1,-0.5]
ub=[1.5,1.25]
[xopt,fopt,exitflag,output,lambda]=linprog(f, A, b, Aeq, beq, lb, ub);
//Linear program with all constraint types, and writing the given problem in MPS/LP file format
f=[-1,-1/3];
A=[1,1;1,1/4;1,-1;-1/4,-1;-1,-1;-1,1];
b=[2,1,2,1,-1,2]
Aeq=[1,1/4]
beq=[1/2]
lb=[-1,-0.5]
ub=[1.5,1.25]
options=list("MaxIter",[200],"WriteMps","ON","WriteLp","ON");
mps="../path to the mps file../.."
lp="../path to the lp file../.."
[xopt,fopt,exitflag,output,lambda]=linprog(f, A, b, Aeq, beq, lb, ub,options,mps,lp);
Primal Infeasible Problem
f=[-1,-1,-1];
A=[1,2,-1];
b=[-4];
Aeq=[1,5,3;1,1,0];
beq=[10,100];
lb=[0,0,0];
ub=[%inf,%inf,%inf];
[xopt,fopt,exitflag,output,lambda]= linprog(f,A,b,Aeq,beq,lb,ub);
Dual Infeasible Problem
f=[3,5,-7];
A=[-1,-1,4;1,1,4];
b=[-8,5];
Aeq=[];
beq=[];
lb=[-%inf,-%inf,-%inf];
ub=[%inf,%inf,%inf];
[xopt,fopt,exitflag,output,lambda]= linprog(f,A,b,Aeq,beq,lb,ub);
//========= E N D === O F === D E M O =========//
