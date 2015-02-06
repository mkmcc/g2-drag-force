(* ============================================================================== *)
(* preliminary definitions *)

gm = 1.9393 * 10^-8; (* GM in pc^3/year^2 *)

Needs["DifferentialEquations`InterpolatingFunctionAnatomy`"];

(* pattern to match only numeric 3-vectors *)
(* - need this to prevent mathematica from trying to pre-evaluate
     expressions with symbolic arguments. *)
ThreeVectQ = Function[{v}, Length[v] == 3 && VectorQ[v, NumericQ]];



(* ============================================================================== *)
(* analysis functions *)
(* - taken from appendix B of Lu et al. 2009 *)

asemi[r_?ThreeVectQ, v_?ThreeVectQ] := 
    (2/Norm[r] - v.v/gm)^(-1);

ecc[r_?ThreeVectQ, v_?ThreeVectQ] := 
    (1 - Norm[Cross[r, v]]^2/(gm asemi[r, v]))^(1/2);

inc[r_?ThreeVectQ, v_?ThreeVectQ] := 
    Block[{j = Cross[r, v]}, 
          ArcCos[-j.{0, 0, 1}/Norm[j]] * 180/Pi];

ascnode[r_?ThreeVectQ, v_?ThreeVectQ] :=
    Block[{j,n},
          j = Evaluate[Cross[r, v]];
          n = Evaluate[Cross[j,{0, 0, 1}]];
          ArcTan[n.{1,0,0}/n.{0,1,0}]*180/Pi];
 
omega[r_?ThreeVectQ, v_?ThreeVectQ] :=
    Block[{j = Cross[r, v], zcrossj, eccvector},
          zcrossj = Cross[{0, 0, 1}, j];
          eccvector = Cross[v,j]/gm - r/Norm[r];
          If[eccvector.{0,0,1} < 0, 
             360 - ArcCos[(zcrossj.eccvector)/Norm[zcrossj]/Norm[eccvector]]*180/Pi, 
             ArcCos[(zcrossj.eccvector)/Norm[zcrossj]/Norm[eccvector]]*180/Pi]];

Print["loaded preliminaries.m"];

(* Local Variables: *)
(* mode: mathematica *)
(* coding: utf-8-unix *)
(* End: *)
