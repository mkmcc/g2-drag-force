(* ============================================================================== *)
(* import astrometry data from the wiki *)
g2posdata = Select[Import["../obs-data/wiki-pos.dat", "Table"],
                 VectorQ[#, NumericQ] &];

g1posdata = Select[Import["../obs-data/pfuhl-g1-pos-lband.dat", "Table"],
                   VectorQ[#, NumericQ] &];

(* convert from arc-sec to pc *)
g2posdata = Map[# * {1, 0.04, 0.04, 0.04, 0.04} &, g2posdata];
g1posdata = Map[# * {1, 0.04, 0.04, 0.04, 0.04} &, g1posdata];



(* repeat for velocity *)
g2veldata = Select[Import["../obs-data/wiki-vr.dat", "Table"],
                   VectorQ[#, NumericQ] &];

g1veldata = Select[Import["../obs-data/pfuhl-g1-vr.dat", "Table"],
                   VectorQ[#, NumericQ] &];

(* convert from km/s to pc/yr *)
g2veldata = Map[# * {1, 1.02269032*10^-6, 1.02269032*10^-6} &, g2veldata];
g1veldata = Map[# * {1, 1.02269032*10^-6, 1.02269032*10^-6} &, g1veldata];



(* ============================================================================== *)
(* fit the observational data with 2nd order polynomials in time, so
   we can interpolate and take derivatives *)
RAfit = With[{time    = g2posdata[[All, 1]], 
              ra      = g2posdata[[All, 2]], 
              raerror = g2posdata[[All, 3]]},
             NonlinearModelFit[Transpose[{time, ra}], 
                               a * t^2 + b * t + c, {a, b, c}, t, 
                               Weights -> 1/raerror^2]];

DECfit = With[{time     = g2posdata[[All, 1]], 
               dec      = g2posdata[[All, 4]], 
               decerror = g2posdata[[All, 5]]},
              NonlinearModelFit[Transpose[{time, dec}], 
                                a * t^2 + b * t + c, {a, b, c}, t, 
                                Weights -> 1/decerror^2]];

VLOSfit = With[{time = g2veldata[[All, 1]], 
                vlos = g2veldata[[All, 2]], 
                err  = g2veldata[[All, 3]]},
               NonlinearModelFit[Transpose[{time, vlos}], 
                                 a * t^2 + b * t + c, {a, b, c}, t, 
                                 Weights -> 1/err^2]];



(* ============================================================================== *)
(* generate an initial condition *)
$t0 = 2013.33;

$icfit = NMinimize[{(With[{r = {RAfit[$t0]  + dra,  DECfit[$t0]  + ddec,  z},
                           v = {RAfit'[$t0] + dvra, DECfit'[$t0] + dvdec, VLOSfit[$t0] + dvlos}},
                          (ecc[r, v]       - 0.9762)^2     / (0.0074)^2
                          + (inc[r, v]     - 118.1)^2      / (2.0)^2
                          + (ascnode[r, v] - 81.9)^2       / (4.3)^2
                          + (omega[r, v]   - 97.2)^2       / (2.0)^2
                          + (asemi[r, v]   - 1.048*0.04)^2 / (0.247*0.04)^2]
                     + (dra)^2   / (0.0003)^2
                     + (ddec)^2  / (0.0003)^2
                     + (dvra)^2  / (3*0.0003)^2
                     + (dvdec)^2 / (3*0.0003)^2) / 5,
                    (*  *)
                    (* constraints *)
                    -0.04 < z    < 0.0,
                    -0.10 < dra  < 0.1,
                    -0.10 < ddec < 0.1},
                   (*  *)
                   (* search space *)
                   {{z,     -0.04,      0.0},
                    {dvra,  -0.003,     0.003},
                    {dvdec, -0.003,     0.003},
                    {dra,   -0.003,     0.003},
                    {ddec,  -0.003,     0.003},
                    {dvlos, -0.0000511, 0.0000511}},
                   Method -> {"SimulatedAnnealing", SearchPoints -> 30}];

(* output some diagnostic information *)
Print["Initial Condition:"];
Print[$icfit];
Print[""];

With[{r = {RAfit[$t0]  + dra,  DECfit[$t0]  + ddec,  z},
      v = {RAfit'[$t0] + dvra, DECfit'[$t0] + dvdec, VLOSfit[$t0] + dvlos}},
     Print["e: "     <> ToString[ecc[r,v]     /. Last[$icfit]]];
     Print["i: "     <> ToString[inc[r,v]     /. Last[$icfit]]];
     Print["a: "     <> ToString[asemi[r,v]   /. Last[$icfit]]];
     Print["Omega: " <> ToString[ascnode[r,v] /. Last[$icfit]]];
     Print["omega: " <> ToString[omega[r,v]   /. Last[$icfit]]];
     Print["r: " <> ToString[r /. Last[$icfit]]];
     Print["v: " <> ToString[v /. Last[$icfit]]];
    ];


Print["loaded initial-cond.m"];

(* Local Variables: *)
(* mode: mathematica *)
(* coding: utf-8-unix *)
(* End: *)
