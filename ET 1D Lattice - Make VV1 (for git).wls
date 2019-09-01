(* ::Package:: *)

(* ::Subsection:: *)
(*1D Elliptic Theta Lattice: V matrix*)


(* ::Subsubsection:: *)
(*Define potential*)


(* ::Input:: *)
(*Clear[\[Epsilon],a,LL,V0]*)
(*prec=MachinePrecision; (* increase for higher accuracy *)*)
(**)
(*Vx[x_]:=EllipticTheta[3,(\[Pi] (x(*-(1/Sqrt[17])a*)))/(2 a),E^(-((\[Pi]^2 \[Epsilon]^2)/(2 a^2)))]*)
(*Vz[z_]:=Exp[-(z(*+(1/Sqrt[13])a*))^2/(2 \[Sigma]^2)]*)
(*V1[x_,z_]:=Evaluate@(Vx[x]Vz[z]-V0)*)
(**)
(*(* set main parameters *)*)
(*\[Epsilon]=SetPrecision[0.4,prec];*)
(*\[Sigma]=SetPrecision[0.4,prec];*)
(*a=SetPrecision[1.0,prec];*)
(*LL=SetPrecision[3.0,prec];*)
(**)
(*V0=0; (* if this changes, need V0 term in code below *)*)
(**)
(*(* check Vmax is same as desired (see above) and Vmin is zero *)*)
(*Print["Vmax=",Vmax=Chop@Maximize[V[x,z],{x,z}][[1]]]*)
(*Print["Vmin=",Vmin=Chop@Minimize[V[x,z],{x,z}][[1]](*," (=",Vmin,"?)"*)]*)
(*Print["Using V0=",V0]*)
(**)
(*(* get saddle potentials, depends on getting correct saddlepts first *)*)
(*saddlePts={{a,0}};*)
(*Print["Vsaddle1=",Vsaddle1=V[saddlePts[[1]]]]*)


(* ::Input:: *)
(**)


(* ::Input:: *)
(*cp=ContourPlot[-14.9603*(V1[x,z]+0.5 V2[x,z]),{x,-3a,3a},{z,-LL,LL},PlotRange->All,LabelStyle->36,Contours->10,FrameLabel->{"x",Rotate["z",-90 Degree]},AspectRatio->1]*)


(* ::Input:: *)
(*linep=Plot[{LL,-LL},{x,-3,3},PlotStyle->Black]*)


(* ::Input:: *)
(*Show[cp,linep,ImageSize->400]*)


(* ::Input:: *)
(**)


(* ::Input:: *)
(*Plot3D[-14.9603*(V1[x,z]+0.5 V2[x,z]),{x,-3a,3a},{z,-LL,LL},PlotRange->All,ColorFunction->(ColorData[{"TemperatureMap",{0,1.2}}][#3]&),PlotStyle->Opacity[0.8],LabelStyle->36,AxesLabel->{"x","z"}]*)


(* ::Subsubsection:: *)
(*Create H matrix, dimension (2Nt+1) x (Nt+1) on each side*)


(* ::Text:: *)
(*Using H_QM=-d^2/dx^2+-d^2/dz^2+V(x,z) and following derivation in Hsu (2005) gives a matrix equation for amplitudes A_nx(\[Nu],k):*)


(* ::Input:: *)
(*VmatS1=VmatS;*)


(* ::Input:: *)
(*Clear[Vxmat,Vymat,Vxymat,Vzmat,VmatS,kx,ky]*)
(*Nt=36;*)
(**)
(*\[Xi][nz_,z_]:=If[nz==0,1/Sqrt[2LL],1/Sqrt[LL] Cos[(nz \[Pi] (z+LL))/(2LL)]]*)
(*Vxmat=Table[*)
(*If[Abs[\[CapitalDelta]nx]<=6,*)
(*NIntegrate[Vx[x] 1/(2a) Exp[I (\[CapitalDelta]nx)\[Pi] x/a],{x,-a,a}]*)
(*(* else *),0],*)
(*{\[CapitalDelta]nx,-2Nt,2Nt}]; //AbsoluteTiming(* for either x or y *)*)
(*Vxmat*)
(**)
(*(* Vzmat is built as a half-matrix to get 2x efficiency *)*)
(*Vzmat=Table[*)
(*If[(*EvenQ@Abs[nz-nzp]*)True,(* for ASYM can't skip any -- CAREFUL *)*)
(*NIntegrate[Vz[z]\[Xi][nz,z]\[Xi][nzp,z],{z,-LL,LL}(*,AccuracyGoal\[Rule]16*)]*)
(*(* else *),0],*)
(*{nz,0,Nt},{nzp,0,Nt}];//AbsoluteTiming*)
(**)
(*(* using V0=0 to skip doing the V0 term -- if V0\[NotEqual]0 I'll need to add that term *)*)
(*VmatS=SparseArray@Flatten[ParallelTable[Flatten@Table[*)
(*(* note no U0 here, it's in "S and Heff" *)*)
(*Vxmat[[nx-nxp+2Nt+1]]Vzmat[[nz+1,nzp+1]]*)
(*,{nx,-Nt,Nt},{nz,0,Nt}],{nxp,-Nt,Nt},{nzp,0,Nt}],1];//Timing//AbsoluteTiming*)
(*(*Vmat//MatrixForm*)*)
(**)
(*Vzmat//MatrixForm*)


(* ::Input:: *)
(*MatrixPlot[(VmatS0-VmatS)[[1;;20,1;;20]],PlotLegends->Automatic]*)


(* ::Input:: *)
(**)


(* ::Input:: *)
(*(* hidden: full H if I want it for some reason *)*)


(* ::Input:: *)
(*KmatS=SparseArray[Table[{n,n},{n,(2Nt+1)(Nt+1)}]->Flatten@Table[*)
(*1/2 ((kx+nx \[Pi]/a)^2+(\[Pi]/(2LL) nz)^2)*)
(*,{nx,-Nt,Nt},{nz,0,Nt}]];//AbsoluteTiming*)
(*(*Kmat//MatrixForm;*)*)
(*U=14.960335515053726;*)
(*HmatS=KmatS-U VmatS;//AbsoluteTiming*)
(*(*Hmat//MatrixForm*)*)


(* ::Input:: *)
(**)


(* ::Text:: *)
(*Saving and Retrieving data*)


(* ::Input:: *)
(*foo=VmatS;*)


(* ::Input:: *)
(*(* save data *)*)
(*SetDirectory["/Users/mollyporter/Documents/UT Austin Stuff/Research (Austin)/Quantum Chaos/data"];*)
(*(*SetDirectory["/home/mporter/Downloads"];*)*)
(**)
(*Export["1D_VV2a_Nt36_L5_Vmat_sparse.mx",VmatS];*)


(* ::Input:: *)
(**)


(* ::Input:: *)
(*(* retrieve data *)*)
(*SetDirectory["/Users/mollyporter/Documents/UT Austin Stuff/Research (Austin)/Quantum Chaos/data"];*)
(*(*SetDirectory["/home/mporter/Downloads"];*)*)
(**)
(*VmatS0=Import["1D_VV2a_Nt36_L5_Vmat_sparse.mx"];//AbsoluteTiming*)
