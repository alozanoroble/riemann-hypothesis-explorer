// quad-char-zeros.m -- zeros of L(s, chi_D) for the two class-number > 1
// fields, as an independent check on the LMFDB values.
SetColumns(0);
PREC := 30; STEP := 1/10;
RR := RealField(PREC); CC := ComplexField(PREC); i := CC.1;

function IllinoisRoot(f, a, b, fa, fb, tol, maxit)
    c := a;
    for n in [1..maxit] do
        if fb eq fa then break; end if;
        c := b - fb*(b-a)/(fb-fa); fc := f(c);
        if fc eq 0 or Abs(b-a) lt tol then return c; end if;
        if fc*fb lt 0 then a := b; fa := fb; else fa := fa/2; end if;
        b := c; fb := fc;
    end for;
    return c;
end function;

for D in [-20, -23] do
    tag := "D" cat IntegerToString(Abs(D));
    chi := KroneckerCharacter(D);
    L := LSeries(chi : Precision := PREC);
    _ := CheckFunctionalEquation(L);
    rt := Sqrt(CC!Sign(L));
    Z := func< t | Real(LStar(L, RR!(1/2) + i*RR!t)/rt) >;
    printf "%o SIGN %o\n", tag, Sign(L);
    zeros := [RR|]; t0 := STEP; f0 := Z(t0);
    while t0 lt 41 do
        t1 := t0 + STEP; f1 := Z(t1);
        if f0*f1 lt 0 then
            Append(~zeros, IllinoisRoot(Z, t0, t1, f0, f1, RR!10^-24, 50));
        end if;
        t0 := t1; f0 := f1;
    end while;
    printf "%o NZEROS %o\n", tag, #zeros;
    for t in zeros do printf "%o ZERO %o\n", tag, t; end for;
end for;
printf "### DONE\n";
quit;
