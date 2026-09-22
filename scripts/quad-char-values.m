// quad-char-values.m -- character tables, Dedekind coefficients, splitting
// behaviour and reference values for the Extended RH panel.
SetColumns(0);
PREC := 30;
RR := RealField(PREC); CC := ComplexField(PREC); i := CC.1;
for d in [-3, -5, -23] do
    K  := QuadraticField(d);
    OK := MaximalOrder(K);
    D  := Discriminant(OK);
    printf "### D%o\n", Abs(D);
    printf "D%o FIELD %o\n", Abs(D), d;
    printf "D%o DISC %o\n", Abs(D), D;
    printf "D%o CLASSNUMBER %o\n", Abs(D), ClassNumber(K);
    printf "D%o CLASSGROUP %o\n", Abs(D), AbelianInvariants(ClassGroup(K));
    printf "D%o ROOTSOFUNITY %o\n", Abs(D), #TorsionUnitGroup(OK);
    chi := KroneckerCharacter(D);
    printf "D%o CHI_ODD %o\n", Abs(D), IsOdd(chi);
    printf "D%o CHI_VALUES %o\n", Abs(D), [chi(n) : n in [1..Abs(D)]];
    L := LSeries(chi : Precision := PREC);
    printf "D%o CHI_L1 %o\n", Abs(D), Evaluate(L, 1);
    printf "D%o CHI_L1_CNF %o\n", Abs(D),
        2*Pi(RR)*ClassNumber(K)/(#TorsionUnitGroup(OK)*Sqrt(RR!Abs(D)));
    printf "D%o CHI_VAL_2 %o\n", Abs(D), Evaluate(L, 2);
    printf "D%o CHI_VAL_halfplus8i %o\n", Abs(D), Evaluate(L, CC!(1/2) + i*8);
    printf "D%o CHI_VAL_1.5plus2i %o\n", Abs(D), Evaluate(L, CC!(3/2) + i*2);
    LK := LSeries(K : Precision := PREC);
    printf "D%o ZETAK_COEFFS %o\n", Abs(D), LGetCoefficients(LK, 40);
    printf "D%o ZETAK_VAL_2 %o\n", Abs(D), Evaluate(LK, 2);
    printf "D%o ZETAK_VAL_halfplus8i %o\n", Abs(D), Evaluate(LK, CC!(1/2) + i*8);
    printf "D%o SPLIT %o\n", Abs(D),
        [<p, [Norm(f[1]) : f in Factorization(p*OK)]> : p in PrimesUpTo(30)];
end for;
printf "### DONE\n";
quit;
