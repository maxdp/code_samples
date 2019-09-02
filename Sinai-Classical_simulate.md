```Mathematica
(* Section: Define potential and parameters *)

prec = MachinePrecision;

(* Change bumpQ for bumps vs holes; True means bumps, False means holes *)
bumpQ = True;
numBumpOrHole = 1;

(* repeat V definition for convenience *)
(* will get shifted below so that Vmin=0.0 *)
If[numBumpOrHole == 1,

  V[x_, y_] := 
   If[bumpQ, 1, -1]*U0 ((EllipticTheta[3, -((pi x)/(2 a)), E^(-((pi^2 ee^2)/(2 a^2)))]
        EllipticTheta[3, -((pi y)/(2 a)), E^(-((pi^2 ee^2)/(2 a^2)))])/(4 a^2) - V0)
        
  ,(*else ==3 presumably*)
  
  V[x_, y_] :=
   If[bumpQ, 1, -1]*
    U0 (( EllipticTheta[3, -((pi x)/(2 a)), E^(-((pi^2 ee^2)/(2 a^2)))]
        EllipticTheta[3, -((pi y)/(2 a)), E^(-((pi^2 ee^2)/(2 a^2)))])/(4 a^2) + 
        4 (EllipticTheta[3, -((pi (x - 1/Sqrt[7]))/(2 a)), E^(-((pi^2 ee^2)/(2 a^2)))]
        EllipticTheta[3, -((pi (y - 1/Sqrt[11]))/(2 a)), E^(-((pi^2 ee^2)/(2 a^2)))])/(4 a^2) +
        (EllipticTheta[3, -((pi (x - 1/Sqrt[13]))/(2 a)), E^(-((pi^2 ee^2)/(2 a^2)))]
        EllipticTheta[3, -((pi (y + 1/Sqrt[5]))/(2 a)), E^(-((pi^2 ee^2)/(2 a^2)))])/(4 a^2) - V0)
  ];
  
  V[{x_, y_}] := V[x, y]

(* set main parameters *)
ee = SetPrecision[0.7, prec];
a = SetPrecision[1.0, prec];

U0 = 1.;(*placeholder, see below*)
V0 = 0;(*placeholder, see below*)

(* find correct V0 and U0 *)
Vmax = 100;(*desired, determines U0*)
V0 = Minimize[V[x, y]/U0, {x, y}][[1]];
U0 = Vmax/Maximize[V[x, y], {x, y}][[1]];

(* set energy and time *)
En = SetPrecision[1, prec]; (* should often be overidden by xmax method below or by EnList *)
dt = Rationalize[0.2]; (* should sometimes be overridden by lyapunov stuff below *)



(* Section: Simulation *)

Clear[En]
EnList = SetPrecision[Table[en, {en, {100}}], prec];
EnModList = 
  Table[EnList[[j]] - If[bumpQ, 0, Vmin], {j, Length@EnList}];

skippedICs = {};

tmax = 10;

Clear[x0, px0, y0, py0, xi0, pxi0, yi0, pyi0, x, y, px, py, x1, y1, t,
   t0, xSample, pxSample];
pointsArray =(*Parallel*)Table[

     EnMod = If[bumpQ, En, En - Vmin];
            
       Clear[sampled, testAbort, abortNow];
       testAbort = False; (* if True, will ABORT after first IC *)
   
           sampled = False; (* do not change *)
       
       abortNow = False; (* do not change *)
       
       Clear[xSample, pxSample];
       
       (* a set of ICs *)
       xpxList=Flatten[Join[
       Table[{j,k},{j,{-0.22,-0.24}},{k,{-4,4}}],
       {}],1];
       
       (* only for a set of ICs *)
       numICs = Length@xpxList;
       {xSample, pxSample} = 
        xpxList[[2]];(*for config space viewing*)
       
       
       If[En == EnList[[1]],
        Print["numICs=", numICs = Length@xpxList](* 
        this numICs won't count if parallel *)
        ];
       
       points = ParallelTable[
         {xi, pxi} = xpxpair;(* 
         necessary for using xpxpair *)
         
         Clear[x, px, y, py, sol, s, ps, ptrans];
         yi = trenchYofX[xi];(*starting on SOS*)
                  
         If[En - Abs[pxi]^2 < V[xi, yi], 
          Null, (*throw out unphysical cases*)(*TESTING!!*)
          
          Clear[pyi];
          pyirules = Solve[pxi^2 + pyi^2 + V[xi, yi] == En, pyi];(*,
          WorkingPrecision\[Rule]highPrec*)
          
          pyi = Chop[Evaluate[pyi /. pyirules[[2]]], 
            10^-15]; (*Using [[2]] because it has pyi>0*)
          
          eqnx = x'[t] == 2 px[t];
          eqnpx = px'[t] == -D[V[x[t], y[t]], x[t]];
          eqny = y'[t] == 2 py[t];
          eqnpy = py'[t] == -D[V[x[t], y[t]], y[t]];
          
          
          sol1 = NDSolve[{eqnx, eqnpx, eqny, eqnpy, x[0] == xi, 
             px[0] == pxi, y[0] == yi, py[0] == pyi}, {x, px, y, 
             py}, {t, 0, tmax}, MaxSteps -> 10000000, 
            MaxStepSize -> 0.001(*,WorkingPrecision\[Rule]prec-1*)];(* 
          SMALL ee or E=Vmax => use MaxStepSize *)
          
          sol = sol1[[1]];
          
          x1[t_] = Evaluate[x[t] /. sol];
          y1[t_] = Evaluate[y[t] /. sol];
          
          x[t_] = Mod[Evaluate[x[t] /. sol], 2 a, -a];
          px[t_] = Evaluate[px[t] /. sol];
          y[t_] = Mod[Evaluate[y[t] /. sol], 2 a, -a];
          py[t_] = Evaluate[py[t] /. sol];
          
          roots = {};
          For[t0 = 0, t0 < tmax, t0 = t0 + dt,
           
           Quiet[  (*CHANGE FOR 3-PEAK VS 1-PEAK AND BUMP VS WELL*)
  
                      root = FindRoot[
              
              Abs@(* for bump *)
                 y[t] - 
                Abs@Mod[trenchYofX[x[t]], 2 a, -a] == 0 (* 
              repeat below *)
              , {t, t0}(*,
              WorkingPrecision\[Rule]Precision@y[t0]*)];
            t1 = Evaluate[t /. root];
            
            
            If[(*t0-dt/2\[LessEqual]t1<t0+dt/2 &&*)
             t1 >= 0 && t1 <= tmax
              && 
              Round[Abs@y[t1] - Abs@Mod[trenchYofX[x[t1]], 2 a, -a], 
                10^-2 a] == 0,(* copy from above, add Round, 
             and change t\[Rule]t1 *)
             
             AppendTo[roots, Chop[t1]]]
            ]
           ];
          roots = Sort[DeleteDuplicates[roots, #1 == #2 &]];
          
          (*sample for debugging*)
          
          If[{xi, pxi} == {xSample, pxSample},

           x0[t_] = x[t]; px0[t_] = px[t]; y0[t_] = y[t]; 
           py0[t_] = py[t];
           xi0 = xi; pxi0 = pxi; yi0 = yi; pyi0 = pyi; 
           roots0 = roots;
           sampled = True;
           abortNow = testAbort;
           ];
          If[abortNow, Abort[]];
          
          pointsPart = {};
          For[i = 1, i <= Length[roots], i++,
           (*If[Chop[Re[py[roots[[i]]]]]>0,(*only use py>
           0 to avoid a mirror image*)
           AppendTo[
           pointsPart,{Chop[Re[x[roots[[i]]]]],Chop[Re[px[roots[[
           i]]]]]/Sqrt[En]}];
           ];*)
           (*borrowed from my Horsley-
           reproducing code*)
           t1 = roots[[i]];
           
           s = trenchSofXmod[
             x[t1]];(*CAREFUL about what SofX function is here*)
      theta1 = ArcTan[dtrench[x[t1]]];
           theta2 = 
            ArcTan[ If[px[t1] == 0, \[Infinity], py[t1]/px[t1]]]
             + 
             If[px[t1] < 0, pi, 
              0]; (*second term corrects for Arctan limited range*)
  
           p = Sqrt[px[t1]^2 + py[t1]^2];
           ps = p Cos[theta2 - theta1];
           ptrans = p Sin[theta2 - theta1];
           (*Print[{"1:",ptrans,s,t1}];*)
           If[ptrans > 0,
            (*Print[{"2:",ptrans,s,t1}];*)
            
            AppendTo[pointsPart, {s, ps(*/Sqrt[EnMod]*)}];
            ]
           ];
          pointsPart
          
          , Null]
          
         , {xpxpair, 
          xpxList}]; // Timing // AbsoluteTiming;
     
     Print["Finished E=", En];
     
     pointsOrig = points;
     points = Flatten[points, 0];(* flatten 0 for xpxpair, 1 for xi,
     pxi *)
     points
     , {En, EnList}]; // Timing // AbsoluteTiming
     
pointsArray = 
  Table[Select[pointsArray[[j]], # =!= Null &], {j, 1, 
    Length[pointsArray]}];



(* Section: Plotting surfaces of section *)

f = 1;
{xf, pxf} = {0, 0 Sqrt[EnModList[[1]]]};
Table[ListPlot[pointsArray[[j]], 
  PlotRange -> {{xf - a/f, xf + a/f},{pxf - 1.1 Sqrt[EnModList[[j]]]/f, 
     pxf + 1.1 Sqrt[EnModList[[j]]]/f}},
  PlotStyle -> {{PointSize[0.01], 
     ColorData[97, "ColorList"][[1]](* enforces blue *)}}, 
  AxesOrigin -> {-a, -Sqrt[EnModList[[j]]]}, 
  PlotLabel -> 
   SetPrecision[
    Row@{"prec=", Round[prec, 0.1], ", E=", EnList[[j]], ", t=", tmax,
       " dt=", dt, ", ICs=", numICs - skippedICs[[j]], ", U0=", 
      U0,", ee=", ee, 
      ", a=", a}, MachinePrecision], Ticks -> {True, True}, 
  FrameTicksStyle -> {Directive[16], Directive[16]}, 
  Frame -> {True, True, True, True}, 
  FrameLabel -> {Style["x", 16], 
    Style["px", 16]}, ImageSize -> 700],
    {j, 1, Length[pointsArray]}]



