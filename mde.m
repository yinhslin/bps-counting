(* ::Package:: *)

(* ::Section:: *)
(*MDE*)


(*
	Execute "math -script mde.m -N " followed by the desired value of N. 
	"math" is the path to Mathematica executable.
*)
options = $CommandLine;

(*
	Uncomment below to specify N directly, suitable when using GUI.
*)
(*options = {"-N", "2"}*)

param[flag_] := Module[
		{position, flagList}
	, 
		flagList = StringSplit[options];
		position = Position[flagList, "-" <> flag];
		Switch[Length[position], 
			1, If[position[[1, 1]] >= Length[flagList], 
				True, 
				flagList[[position[[1, 1]]+1]][[1]]
				], 
			0, Null, 
			_, WriteString["stdout", "flag -" <> flag <> " duplicated\n"];
				Abort[]
		]
	];

NN = param["N"] // ToExpression;


(* ::Subsection:: *)
(*Parameters*)


cc=-3(NN^2-1);
(* degree of MDE *)
n=((NN+1)/2)^2;

(* basis of weigh d holomorphic modular forms as monomials in E4 E6 *)
EisParitions[d_]:=EisParitions[d]=Product[e[n],{n,#}]&/@IntegerPartitions[d,Infinity,{4,6}];
NumVars[d_]:=Sum[Length[EisParitions[w]],{w,2,d,2}];

(* can choose to overconstrain as sanity check *)
cut=-cc/24+NumVars[2n]+1;


(* ::Subsection:: *)
(*Grand odd Schur*)


(* ::Text:: *)
(*2205.00818 (2.14)*)


Timing[F=I EllipticTheta[1,z,Sqrt[q+O[q]^(cut+2)]];
F=F//TrigToExp;
yx=1/2 (x+Sqrt[4+x^2])+O[x]^cut;
pre=1/DedekindEta[Log[q+O[q]^(cut+2)]/(2\[Pi] I)]^3;
F=F/.z->Log[y]/I/.y->yx;
index[NN_]:=Series[(-1)^((NN-1)/2)CoefficientList[Normal[F],x][[NN+1]]*pre,{q,0,cut}]/;OddQ[NN];][[1]]//Print;


Timing[f=index[NN];][[1]]//Print;


(* ::Subsection:: *)
(*MDE*)


Clear[ES];
ES[n_]:=ES[n]=1-Series[2n/BernoulliB[n]Sum[m^(n-1)q^(m d),{d,1,cut},{m,1,cut}],{q,0,cut}];
repl=e[n_]:>ES[n];

Clear[qD,covD]
qD[expr_,n_]:=q D[qD[expr,n-1],q]/;n>0;
qD[expr_,0]:=expr;
covD[expr_,n_]:=(covD[expr,n]=qD[covD[expr,n-1],1]-2(n-1)/12ES[2]covD[expr,n-1])/;n>=1;
covD[expr_,0] := expr;


\[Alpha][0,1]:=1;
Timing[mde=q^(cc/24)*Sum[(covD[f,n-t]/.repl)Sum[\[Alpha][t,i]EisParitions[2t][[i]]/.repl,{i,1,Length[EisParitions[2t]]}],{t,0,n}]//Expand;][[1]]//Print;


Timing[eqns=CoefficientList[mde,q];
vars=Variables[eqns];
sol=Solve[#==0&/@eqns,vars][[1]];][[1]]//Print;


(* ::Subsection:: *)
(*Indicial*)


covDInd[n_]:=Product[(h-cc/24-2(i-1)/12),{i,1,n}];
Timing[ind=Sum[(covDInd[n-t])Sum[(\[Alpha][t,i]/.sol),{i,1,Length[EisParitions[2t]]}],{t,0,n}]//Expand;][[1]]//Print;
Solve[ind==0,h]//Print;
