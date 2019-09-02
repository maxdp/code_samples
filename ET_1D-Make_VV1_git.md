```Mathematica
(* Section: Define potential *)

Clear[ee, a, LL, V0]
prec = MachinePrecision; (* increase for higher accuracy *)

(* define potential; will need x and z separate *)
Vx[x_] := 
 EllipticTheta[3, (pi x)/(2 a), 
  E^(-((pi^2 ee^2)/(2 a^2)))]
Vz[z_] := Exp[-z^2/(2 sigma^2)]
V1[x_, z_] := Evaluate@(Vx[x] Vz[z] - V0)

(* set main parameters *)
ee = SetPrecision[0.4, prec];
sigma = SetPrecision[0.4, prec];
a = SetPrecision[1.0, prec];
LL = SetPrecision[3.0, prec];

V0 = 0; (* if this changes, need V0 term in code below *)

(* check Vmax is same as desired (see above) and Vmin is zero *)
Print["Vmax=", Vmax = Chop@Maximize[V[x, z], {x, z}][[1]]]
Print["Vmin=", 
 Vmin = Chop@Minimize[V[x, z], {x, z}][[1]](*," (=",Vmin,"?)"*)]
Print["Using V0=", V0]

(* get saddle potentials, depends on getting correct saddlepts first *)
saddlePts = {{a, 0}};
Print["Vsaddle1=", Vsaddle1 = V[saddlePts[[1]]]]



(* Section: Create H matrix, dimension (2Nt+1) x (Nt+1) on each side *)

(* Using H_QM=-d^2/dx^2+-d^2/dz^2+V(x,z) and following derivation in Hsu and Reichl (2005) gives a matrix equation for amplitudes A_nx(nu,k) *)

Clear[Vxmat, Vymat, Vxymat, Vzmat, VmatS, kx, ky]
Nt = 6;

(* calculate x integrals *)
Xi[nz_, z_] := 
 If[nz == 0, 1/Sqrt[2 LL], 
  1/Sqrt[LL] Cos[(nz pi (z + LL))/(2 LL)]]
Vxmat = Table[
    If[Abs[dnx] <= 6,
     NIntegrate[
      Vx[x] 1/(2 a) Exp[I (dnx) pi x/a], {x, -a, a}]
     (* else *), 0],
    {dnx, -2 Nt, 
     2 Nt}]; // AbsoluteTiming (* for either x or y *)
Vxmat //MatrixForm (* print *)

(* calculate z integrals *)
(* Vzmat is built as a half-matrix to get 2x efficiency *)
Vzmat = Table[
    If[EvenQ@Abs[nz-nzp],
     NIntegrate[Vz[z] Xi[nz, z] Xi[nzp, z], {z, -LL, LL}(*,
      AccuracyGoal->16*)]
     (* else *), 0],
    {nz, 0, Nt}, {nzp, 0, Nt}]; // AbsoluteTiming
Vzmat //MatrixForm (* print *)

(* combine to make whole V-matrix *)
(* using V0=0 to skip doing the V0 term -- if V0!=0 I'll need to add that term *)
VmatS = SparseArray@Flatten[ParallelTable[Flatten@Table[
         Vxmat[[nx - nxp + 2 Nt + 1]] Vzmat[[nz + 1, nzp + 1]]
         , {nx, -Nt, Nt}, {nz, 0, Nt}], {nxp, -Nt, Nt}, {nzp, 0, Nt}],
       1]; // Timing // AbsoluteTiming



(* Section: Saving and retrieving data *)

(* save data *)

SetDirectory[
  "/Users/mollyporter/Documents/UT Austin Stuff/Research \
(Austin)/Quantum Chaos/data"];

Export["1D_VV1_Nt6_L3_Vmat_sparse.mx", VmatS];

(* retrieve data *)

SetDirectory[
  "/Users/mollyporter/Documents/UT Austin Stuff/Research \
(Austin)/Quantum Chaos/data"];

VmatS0 = Import["1D_VV1_Nt6_L3_Vmat_sparse.mx"]; // AbsoluteTiming


