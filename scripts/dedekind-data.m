// dedekind-data.m
//
// Data for the Extended RH panel: Dedekind zeta functions of imaginary
// quadratic fields with class number bigger than 1, and the quadratic
// Dirichlet characters they factor through.
//
// Run:
//   magma -b scripts/dedekind-data.m
//
// Zero-finding is the same sign-change method as scripts/lfunction-data.m;
// see that file for why Magma's own Zeros intrinsic cannot be used here.

SetColumns(0);
va, vb, vc := GetVersion();
printf "VERSION %o.%o-%o\n", va, vb, vc;

PREC   := 30;   // 40 exhausted memory on the conductor-20 scan; 30 completes
NCOEFF := 40;
STEP   := 1/10;

RR := RealField(PREC);
CC := ComplexField(PREC);
i  := CC.1;

function IllinoisRoot(f, a, b, fa, fb, tol, maxit)
    c := a;
    for n in [1..maxit] do
        if fb eq fa then break; end if;
        c := b - fb*(b-a)/(fb-fa);
        fc := f(c);
        if fc eq 0 or Abs(b-a) lt tol then return c; end if;
        if fc*fb lt 0 then a := b; fa := fb; else fa := fa/2; end if;
        b := c; fb := fc;
    end for;
    return c;
end function;

function CriticalLineZeros(L, height, tag)
    m      := MotivicWeight(L);
    centre := RR!(m+1)/2;
    fe     := CheckFunctionalEquation(L);
    printf "%o FECHECK %o\n", tag, fe;
    eps    := CC!Sign(L);
    error if eps eq 0, "root number unknown for " cat tag;
    printf "%o SIGN %o\n", tag, eps;
    printf "%o DEGREE %o\n", tag, Degree(L);
    rt := Sqrt(eps);
    Z  := func< t | Real(LStar(L, centre + i*RR!t)/rt) >;

    // The Z-function itself, for validating the app's line plot.
    for t in [RR | 1, 3, 6, 9, 14, 20 ] do
        printf "%o ZVAL %o %o\n", tag, t, Z(t);
    end for;

    zeros := [RR|];
    t0 := STEP; f0 := Z(t0);
    while t0 lt height do
        t1 := t0 + STEP; f1 := Z(t1);
        if f0*f1 lt 0 then
            Append(~zeros, IllinoisRoot(Z, t0, t1, f0, f1, RR!10^-30, 60));
        end if;
        t0 := t1; f0 := f1;
    end while;
    printf "%o NZEROS %o\n", tag, #zeros;
    for t in zeros do printf "%o ZERO %o\n", tag, t; end for;
    return zeros;
end function;

// Two imaginary quadratic fields where unique factorisation fails.
for spec in [ <-5, "D20">, <-23, "D23"> ] do
    d   := spec[1];
    tag := spec[2];
    K   := QuadraticField(d);
    OK  := MaximalOrder(K);
    D   := Discriminant(OK);

    printf "### %o\n", tag;
    printf "%o FIELD Q(sqrt(%o))\n", tag, d;
    printf "%o DISC %o\n", tag, D;
    printf "%o CLASSNUMBER %o\n", tag, ClassNumber(K);
    printf "%o CLASSGROUP %o\n", tag, AbelianInvariants(ClassGroup(K));
    printf "%o UNITS %o\n", tag, #TorsionUnitGroup(OK);

    // the quadratic character attached to K
    chi := KroneckerCharacter(D);
    printf "%o CHI_CONDUCTOR %o\n", tag, Conductor(chi);
    printf "%o CHI_ODD %o\n", tag, IsOdd(chi);
    printf "%o CHI_VALUES %o\n", tag, [chi(n) : n in [1..Abs(D)]];

    Lchi := LSeries(chi : Precision := PREC);
    printf "%o CHI_COEFFS %o\n", tag, LGetCoefficients(Lchi, NCOEFF);
    printf "%o CHI_L1 %o\n", tag, Evaluate(Lchi, 1);
    printf "%o CHI_L1_CLASSNUMBERFORMULA %o\n", tag,
           2*Pi(RR)*ClassNumber(K) / (#TorsionUnitGroup(OK) * Sqrt(RR!Abs(D)));
    printf "%o CHI_VAL_2 %o\n", tag, Evaluate(Lchi, 2);
    printf "%o CHI_VAL_halfplus8i %o\n", tag, Evaluate(Lchi, CC!(1/2) + i*8);
    printf "%o CHI_VAL_1.5plus2i %o\n", tag, Evaluate(Lchi, CC!(3/2) + i*2);
    _ := CriticalLineZeros(Lchi, 40, tag cat "_CHI");

    // the Dedekind zeta function itself
    LK := LSeries(K : Precision := PREC);
    printf "%o ZETAK_COEFFS %o\n", tag, LGetCoefficients(LK, NCOEFF);
    printf "%o ZETAK_VAL_2 %o\n", tag, Evaluate(LK, 2);
    printf "%o ZETAK_VAL_halfplus8i %o\n", tag, Evaluate(LK, CC!(1/2) + i*8);
    printf "%o ZETAK_VAL_1.5plus2i %o\n", tag, Evaluate(LK, CC!(3/2) + i*2);
    // No residue call here: zeta_K has a pole at s = 1 and asking Magma to
    // evaluate the completed function there aborts the run. The residue is
    // L(1, chi_D) anyway, since zeta contributes residue 1, and that is
    // printed above next to the class number formula it must equal.

    // how the rational primes split, for the Euler-product section
    printf "%o SPLIT %o\n", tag,
           [<p, [<Norm(f[1]), f[2]> : f in Factorisation(p*OK)]> : p in PrimesUpTo(30)];
end for;

// Also the Z-function of zeta itself, to validate the app's line plot there.
printf "### ZETA\n";
L1 := RiemannZeta(: Precision := PREC);
_ := CriticalLineZeros(L1, 1/2, "ZETA");

printf "### E11A1\n";
L3 := LSeries(EllipticCurve("11a1") : Precision := PREC);
_ := CriticalLineZeros(L3, 1/2, "E11A1");

printf "### E37A1\n";
L4 := LSeries(EllipticCurve("37a1") : Precision := PREC);
_ := CriticalLineZeros(L4, 1/2, "E37A1");

printf "### CHI3\n";
L2 := LSeries(KroneckerCharacter(-3) : Precision := PREC);
_ := CriticalLineZeros(L2, 1/2, "CHI3");

printf "### DONE\n";
quit;
