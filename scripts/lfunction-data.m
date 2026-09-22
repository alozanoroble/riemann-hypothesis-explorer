// lfunction-data.m
//
// Produces the data the RH exploration app needs:
//   - Dirichlet coefficients a_n
//   - zeros on the critical line
//   - reference L-values used to validate the app's JavaScript evaluators
//
// Run:
//   magma -b scripts/lfunction-data.m
//
// Magma 2.29-10 has no zero-finding intrinsic for the LSer type -- `Zeros`
// only accepts curves, function fields and number fields. So zeros are found
// here from the completed L-function instead.
//
// The method: if L is self-dual with real coefficients, motivic weight m and
// root number eps, then Lambda(s) = eps * Lambda(m+1-s), and on the critical
// line s = (m+1)/2 + it this forces
//
//     Lambda(c+it) / sqrt(eps)   to be real-valued.
//
// Its sign changes are exactly the odd-order zeros on the line, so a scan plus
// the Illinois root-finder locates them. Zeros off the line (if RH failed)
// would be invisible to this method -- that is a statement about what the app
// displays, not evidence for RH.
//
// No randomness is used, so no seed is set.

SetColumns(0);
va, vb, vc := GetVersion();
printf "VERSION %o.%o-%o\n", va, vb, vc;

PREC   := 40;   // working precision, digits
NCOEFF := 40;   // Dirichlet coefficients to tabulate
STEP   := 1/10; // critical-line scan step

RR := RealField(PREC);
CC := ComplexField(PREC);
i  := CC.1;

printf "PRECISION %o\n", PREC;

// ---------------------------------------------------------------- helpers

// Illinois variant of regula falsi: keeps the root bracketed (unlike the
// secant method) but converges far faster than bisection, which matters
// because every call to f is a completed-L-function evaluation.
function IllinoisRoot(f, a, b, fa, fb, tol, maxit)
    c := a;
    for n in [1..maxit] do
        if fb eq fa then break; end if;
        c := b - fb*(b-a)/(fb-fa);
        fc := f(c);
        if fc eq 0 or Abs(b-a) lt tol then return c; end if;
        if fc*fb lt 0 then
            a := b; fa := fb;
        else
            fa := fa/2;
        end if;
        b := c; fb := fc;
    end for;
    return c;
end function;

// Scan the critical line for sign changes and refine each one.
function CriticalLineZeros(L, height, tag)
    m      := MotivicWeight(L);
    center := RR!(m+1)/2;
    // Sign(L) is 0 until the functional equation has actually been checked,
    // which would make the rotation below divide by zero.
    fecheck := CheckFunctionalEquation(L);
    printf "%o FECHECK %o\n", tag, fecheck;
    eps    := CC!Sign(L);
    error if eps eq 0, "root number still unknown for " cat tag;
    rt     := Sqrt(eps);
    Z      := func< t | Real(LStar(L, center + i*RR!t)/rt) >;

    printf "%o CENTRE %o\n", tag, center;
    printf "%o WEIGHT %o\n", tag, m;
    printf "%o SIGN %o\n", tag, eps;
    printf "%o DEGREE %o\n", tag, Degree(L);

    zeros := [RR|];
    t0 := STEP;          // start just off the real axis
    f0 := Z(t0);
    // A zero at the central point itself (analytic rank > 0) is not found by
    // a sign change, so test it separately.
    centralval := Abs(Evaluate(L, center));
    printf "%o CENTRALVALUE %o\n", tag, centralval;

    while t0 lt height do
        t1 := t0 + STEP;
        f1 := Z(t1);
        if f0*f1 lt 0 then
            r := IllinoisRoot(Z, t0, t1, f0, f1, RR!10^-30, 60);
            Append(~zeros, r);
        end if;
        t0 := t1; f0 := f1;
    end while;
    return zeros;
end function;

procedure ShowZeros(zeros, tag)
    printf "%o NZEROS %o\n", tag, #zeros;
    for t in zeros do
        printf "%o ZERO %o\n", tag, t;
    end for;
end procedure;

// ---------------------------------------------------------------- zeta
printf "### ZETA\n";
L1 := RiemannZeta(: Precision := PREC);
printf "ZETA COEFFS %o\n", LGetCoefficients(L1, NCOEFF);
z1 := CriticalLineZeros(L1, 78, "ZETA");
ShowZeros(z1, "ZETA");
printf "ZETA VAL 2 %o\n", Evaluate(L1, 2);
printf "ZETA VAL 3 %o\n", Evaluate(L1, 3);
printf "ZETA VAL -1 %o\n", Evaluate(L1, -1);
printf "ZETA VAL 0.5+10i %o\n", Evaluate(L1, CC!(1/2) + i*10);
printf "ZETA VAL 1.5+2i %o\n", Evaluate(L1, CC!(3/2) + i*2);
printf "ZETA VAL 0.25+30i %o\n", Evaluate(L1, CC!(1/4) + i*30);

// ------------------------------------------------- Dirichlet L(s, chi_-3)
// The quadratic character attached to Q(sqrt(-3)): conductor 3, odd, order 2.
printf "### CHI3\n";
chi := KroneckerCharacter(-3);
printf "CHI3 CONDUCTOR %o\n", Conductor(chi);
printf "CHI3 ORDER %o\n", Order(chi);
printf "CHI3 ODD %o\n", IsOdd(chi);
printf "CHI3 VALUES %o\n", [<n, chi(n)> : n in [1..12]];
L2 := LSeries(chi : Precision := PREC);
printf "CHI3 COEFFS %o\n", LGetCoefficients(L2, NCOEFF);
z2 := CriticalLineZeros(L2, 50, "CHI3");
ShowZeros(z2, "CHI3");
printf "CHI3 VAL 1 %o\n", Evaluate(L2, 1);
printf "CHI3 VAL 2 %o\n", Evaluate(L2, 2);
printf "CHI3 VAL 0.5+8i %o\n", Evaluate(L2, CC!(1/2) + i*8);
printf "CHI3 VAL 1.5+2i %o\n", Evaluate(L2, CC!(3/2) + i*2);

// ------------------------------------------------- elliptic curve 11a1
printf "### E11A1\n";
E := EllipticCurve("11a1");
printf "E11A1 MODEL %o\n", aInvariants(E);
printf "E11A1 CONDUCTOR %o\n", Conductor(E);
printf "E11A1 RANK %o\n", Rank(E);
L3 := LSeries(E : Precision := PREC);
printf "E11A1 COEFFS %o\n", LGetCoefficients(L3, NCOEFF);
printf "E11A1 AP %o\n", [<p, TraceOfFrobenius(E, p)> : p in PrimesUpTo(40)
                          | Conductor(E) mod p ne 0];
z3 := CriticalLineZeros(L3, 35, "E11A1");
ShowZeros(z3, "E11A1");
printf "E11A1 VAL 1 %o\n", Evaluate(L3, 1);
printf "E11A1 VAL 2 %o\n", Evaluate(L3, 2);
printf "E11A1 VAL 1+8i %o\n", Evaluate(L3, CC!1 + i*8);
printf "E11A1 VAL 1.5+2i %o\n", Evaluate(L3, CC!(3/2) + i*2);

// ------------------------------------------------- elliptic curve 37a1
// Rank 1, so L(E,1) = 0: the contrast with 11a1 is the pedagogical point.
printf "### E37A1\n";
E2 := EllipticCurve("37a1");
printf "E37A1 MODEL %o\n", aInvariants(E2);
printf "E37A1 CONDUCTOR %o\n", Conductor(E2);
printf "E37A1 RANK %o\n", Rank(E2);
L4 := LSeries(E2 : Precision := PREC);
printf "E37A1 COEFFS %o\n", LGetCoefficients(L4, NCOEFF);
z4 := CriticalLineZeros(L4, 35, "E37A1");
ShowZeros(z4, "E37A1");
printf "E37A1 VAL 1 %o\n", Evaluate(L4, 1);
printf "E37A1 VAL 2 %o\n", Evaluate(L4, 2);

printf "### DONE\n";
quit;
