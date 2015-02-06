(* ==============================================================================*)
(* modelling function.  this integrates the equation of motion *)

mkorbit[modelParams_?VectorQ, tlim_: 2100] :=
    Module[{sol, tbeg, tend,
            (* helper functions for equation of motion *)
            vbg, rho, M, vrel, dragForce,

            (* physical parameters for the cloud and background *)
            J = {Sin[theta] * Cos[phi],
                 Sin[theta] * Sin[phi],
                 Cos[theta]} /. modelParams,

            r0 = 0.04},

           (* background velocity as a function of 3D radius *)
           vbg[r_?ThreeVectQ] := (fkep * Sqrt[gm/Norm[r]]
                                  * Cross[J, r]/(Norm[r]*Norm[J])) /. modelParams;

           vrel[r_?ThreeVectQ, v_?ThreeVectQ] := (v - vbg[r]);


           (* background density as a function of 3D radius *)
           rho[r_?ThreeVectQ] := (1.65*10^5)*(Norm[r]/r0)^(-alpha) /. modelParams;

           (* return v^2 * area/R^2 as an {x,y,z} vector *)
           (* - assume the cloud is tidally elongated in its direction of motion,
                so we caclculate the force in a coordinate system aligned with v. 
                (not v_rel).  then transform back *)
           dragForce[r_?ThreeVectQ, v_?ThreeVectQ] :=
           With[{rot = RotationMatrix[{v, {1,0,0}}],
                 vr  = vrel[r, v]},
                Module[{ram = Sign[rot.vr]*(rot.vr)^2,
                        area = {1, 1+gf, 1+gf}},
                      Inverse[rot].(ram * area)]] /. modelParams;

           (* cloud Mach # as a function of 3D radius and 3D velocity *)
           M[r_?ThreeVectQ, v_?ThreeVectQ] := Norm[v - vbg[r]]/Sqrt[(5/3)*gm/Norm[r]];

           sol = NDSolve[{r'[t] == v[t],
                          (*  *)
                          v'[t] == -(gm/Norm[r[t]]^3)*r[t]
                                   - (rho[r[t]])/(10^(-10) + sigmag2) * (1 + 2/(0.001 + beta*M[r[t], v[t]]^2))
                                   * (dragForce[r[t], v[t]]),
                          (*  *)
                          v[$t0] == {RAfit'[$t0]  + dvra  + vraShift,
                                     DECfit'[$t0] + dvdec + vdecShift,
                                     VLOSfit[$t0] + dvlos},
                          (*  *)
                          r[$t0] == {RAfit[$t0]  + dra  + raShift,
                                     DECfit[$t0] + ddec + decShift,
                                     z}} /. Last[$icfit] /. modelParams,
                         {r, v},
                         {t, 1000, tlim},
                         (* *)
                         (* cut off integration at small radii so NDSolve
                           doesn't return an error, and also at apocenter,
                           so we don't waste time on past behavior *)
                         Method -> {"EventLocator",
                                    "Event" -> {(Norm[r[t]] - 10^-4),
                                                UnitStep[$t0 - t] (r[t].r'[t])},
                                    "Direction" -> {1, 1}}] // First;

           {tbeg, tend} = First[InterpolatingFunctionDomain[r /. sol]];
           {sol, tbeg, tend}];

Print["loaded mk-orbit.m"];

(* Local Variables: *)
(* mode: mathematica *)
(* coding: utf-8-unix *)
(* End: *)
