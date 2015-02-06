Needs["ErrorBarPlots`"];

Get["preliminaries.m"];
Get["initial-cond.m"];
Get["mk-orbit.m"];
Get["chi-sq.m"];

Get["../data/obserr.symb"];

(* find the best-fitting orbit from the minimization... *)
(*  *)
(* $params = SortBy[$params, (chisq /. #) &]; *)
(* $uniqparams = Map[First, GatherBy[$params, Round[(chisq /. #), 0.01] &]]; *)

(* Print[Length[$uniqparams]]; *)

(* sol = mkorbit[$uniqparams[[1]], 2200]; *)


(* or, just type it in by hand *)
(*  *)
bestfit = {chisq     -> 195.516755567,
           alpha     ->   1.129935682889404047,
           beta      ->   3.524593420717720971,
           fkep      ->   9.087227016958496773*^-1,
           theta     ->   1.678734356242794323,
           phi       ->  -6.273240605021621619*^-1,
           gf        ->   1.096839254290853560*^1,
           sigmag2   ->   1.197338366599284811*^6,
           raShift   ->   6.930575459393565089*^-5,
           decShift  ->   3.380398930205274198*^-6,
           vraShift  ->   7.082737032032123347*^-5,
           vdecShift ->   4.856761860096908484*^-5,
           dtg1      ->  12.4869871196};

sol = mkorbit[bestfit, 2300];


(* use mathematica's plotting function to pick the points *)
(*  *)
p1 = Apply[With[{sol=#1,tbeg=#2,tend=#3},
                ParametricPlot3D[r[t]/.sol, {t,tbeg,tend},
                                 MaxRecursion->10]]&,
           sol];

data = Cases[p1, Line[data_]:>data,-4,1][[1]];

Print[Length[data]];

Export["orbit.3d.dat",
       Append[data,{}],
       "Table"];


(* also save a rotated system in which G2 stars off on the x-axis *)
(*  *)
mat = RotationTransform[{r[sol[[2]]]/.sol[[1]], {-1,0,0}}];

rotdata = Map[mat, data];

Export["orbit-rot.3d.dat",
       Append[rotdata, {}],
       "Table"];



(* and, save velocity information *)
(*  *)
ttg1 = dtg1/.bestfit;

p1 = Plot[Evaluate[v[t] /. sol[[1]]][[3]], {t, sol[[2]], sol[[3]]}];

data = Cases[p1, x_Line :> First@x, Infinity] // First;
data = Prepend[data, "# tg1 = " <> ToString[FortranForm[ttg1]]];

Export["orbit-v.dat", Append[data, {}], "Table"];
