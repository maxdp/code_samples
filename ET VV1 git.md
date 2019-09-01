```Mathematica
Clear[\[Epsilon], a, LL, V0]
prec = MachinePrecision; (* increase for higher accuracy *)

Vx[x_] := 
 EllipticTheta[3, (\[Pi] (x(*-(1/Sqrt[17])a*)))/(2 a), 
  E^(-((\[Pi]^2 \[Epsilon]^2)/(2 a^2)))]
Vz[z_] := Exp[-(z(*+(1/Sqrt[13])a*))^2/(2 \[Sigma]^2)]
V1[x_, z_] := Evaluate@(Vx[x] Vz[z] - V0)

(* set main parameters *)
\[Epsilon] = SetPrecision[0.4, prec];
\[Sigma] = SetPrecision[0.4, prec];
a = SetPrecision[1.0, prec];
LL = SetPrecision[3.0, prec];

V0 = 0; (* if this changes, need V0 term in code below *)

(* check Vmax is same as desired (see above) and Vmin is zero *)
\
Print["Vmax=", Vmax = Chop@Maximize[V[x, z], {x, z}][[1]]]
Print["Vmin=", 
 Vmin = Chop@Minimize[V[x, z], {x, z}][[1]](*," (=",Vmin,"?)"*)]
Print["Using V0=", V0]

(* get saddle potentials, depends on getting correct saddlepts first *)

saddlePts = {{a, 0}};
Print["Vsaddle1=", Vsaddle1 = V[saddlePts[[1]]]]