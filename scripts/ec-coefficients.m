// ec-coefficients.m -- Dirichlet coefficients a_n of the two example curves,
// far enough out for the app's smoothed-sum evaluator to converge.
SetColumns(0);
for lab in ["11a1", "37a1"] do
    E := EllipticCurve(lab);
    L := LSeries(E : Precision := 20);
    printf "%o %o\n", lab, LGetCoefficients(L, 70);
end for;
quit;
