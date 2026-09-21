# lfunction-zeros.sage
#
# Zeros of the three example L-functions used by the RH exploration app,
# computed in Sage through Rubinstein's lcalc.
#
# Run:  sage lfunction-zeros.sage
#
# NOTE ON PROVENANCE: this script follows the documented examples in the Sage
# reference manual (sage.lfunctions.lcalc and
# sage.schemes.elliptic_curves.lseries_ell). It has NOT been executed, because
# the machine this project was built on has no Sage installed -- only Magma.
# The Magma companion, scripts/lfunction-data.m, is the one that actually
# produced the numbers in results/ and in the app. If you run this and it
# disagrees with those, trust neither until you have looked.

print("=" * 64)
print("(1) The Riemann zeta function")
print("=" * 64)

# lcalc.zeros(n) returns the first n zeros of zeta by default.
# Documented output: [14.1347251, 21.0220396, 25.0108576, 30.4248761]
for t in lcalc.zeros(20):
    print(t)

# Zeros in a window rather than the first few. The second entry of each pair
# is lcalc's own diagnostic, not part of the zero.
print("\nzeros with 10 < t < 30:")
print(lcalc.zeros_in_interval(10, 30, 0.1))


print("=" * 64)
print("(2) The Dirichlet L-function of Q(sqrt(-3))")
print("=" * 64)

from sage.libs.lcalc.lcalc_Lfunction import *

chi = DirichletGroup(3).0          # the quadratic character mod 3
print("conductor :", chi.conductor())
print("order     :", chi.order())
print("odd       :", chi.is_odd())
print("values    :", [chi(n) for n in range(1, 13)])

L = Lfunction_from_character(chi, type='int')
print("\nfirst 17 zeros:")
for t in L.find_zeros_via_N(17):
    print(t)

# L(1, chi) is a class number in disguise: for this character the analytic
# class number formula gives 2*pi*h / (w*sqrt(3)) = pi/(3*sqrt(3)).
print("\nL(1, chi)   =", RR(pi / (3 * sqrt(3))), " (class number formula)")


print("=" * 64)
print("(3) Dedekind zeta functions: fields with class number > 1")
print("=" * 64)

# zeta_K = zeta * L(s, chi_D) for a quadratic field of discriminant D, so the
# zeros of zeta_K are the zeros of zeta together with those of L(s, chi_D).
for d in [-5, -23]:
    K = QuadraticField(d)
    D = K.discriminant()
    print("\n--- Q(sqrt(%s)):  disc %s,  h = %s,  class group %s"
          % (d, D, K.class_number(), K.class_group().invariants()))
    print("zeta_K(2)     :", K.zeta_function()(2))

    chi = kronecker_character(D)
    L = Lfunction_from_character(chi, type='int')
    print("zeros of L(s, chi_D), which are the ones zeta does not supply:")
    for t in L.find_zeros_via_N(12):
        print("   ", t)

    # residue of zeta_K at s = 1, i.e. the analytic class number formula
    h = K.class_number()
    w = K.number_of_roots_of_unity()
    print("2*pi*h/(w*sqrt|D|) =", RR(2 * pi * h / (w * sqrt(abs(D)))))
    print("L(1, chi_D)        =", RR(L.value(1).real()) if hasattr(L, 'value')
          else "(compare with the line above)")


print("=" * 64)
print("(4) Elliptic curve L-functions")
print("=" * 64)

for label in ['11a1', '37a1']:
    E = EllipticCurve(label)
    print("\n--- %s: %s" % (label, E))
    print("conductor     :", E.conductor())
    print("rank          :", E.rank())
    print("root number   :", E.root_number())
    print("a_n, n <= 20  :", E.anlist(20)[1:])
    print("L(E, 1)       :", E.lseries().dokchitser()(1))
    # For a rank 1 curve the first zero reported sits at t = 0, which is the
    # central point s = 1 itself.
    print("first 10 zeros:")
    for t in E.lseries().zeros(10):
        print("   ", t)
