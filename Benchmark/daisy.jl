# ==========
# Daisy
# ==========

SUITE["Daisy"] = BenchmarkGroup()

# These examples are taken from [Project Daisy](https://github.com/malyzajko/daisy/blob/master/testcases/).

# This function measures the relative precision of the result in a more informative way than
# taking the scalar overestimation because it evaluates the precision of the lower and the
# upper range bounds separately, see Eq. (20) in [1].
function relative_precision(x, x_ref)
    x_low, x_high = inf(x), sup(x)
    x_ref_low, x_ref_high = inf(x_ref), sup(x_ref)
    rel_low = -(x_low-x_ref_low)/(x_high-x_ref_low)
    rel_high = (x_high-x_ref_high)/(x_high-x_ref_low)
    return 100 * Interval(rel_low, rel_high)
end

# holds the vector of relative precision intervals for each benchmark
RP = Vector{Interval}()

# The following intervals are used throughout the tests in Tables 3-5 in
# [1] Althoff, M., Grebenyuk, D., & Kochdumper, N. (2018). Implementation of Taylor models in CORA 2018.
#     In Proc. of the 5th International Workshop on Applied Verification for Continuous and Hybrid Systems.

a = Interval(-4.5, -0.3)
b = Interval(0.4, 0.9)
c = Interval(3.8, 7.8)
d = Interval(8.0, 10.0)
e = Interval(-10.0, 8.0)
f = Interval(1.0, 2.0)

# ==========
# sine
# ==========

dom = a
x = Taylor1(7)
p = x - (x*x*x)/6.0 + (x*x*x*x*x)/120.0 - (x*x*x*x*x*x*x)/5040.0
SUITE["Daisy"]["sin - evaluate"] = @benchmarkable evaluate($p, $dom)

approx = evaluate(p, dom)
ref = Interval(-1.0002065, 2.72505)
relprec = relative_precision(approx, ref)
push!(RP, relprec)

# ==========
# bspline
# ==========

dom = a
u = Taylor1(3)

p = (1 - u) * (1 - u) * (1 - u) / 6.0

SUITE["Daisy"]["bspline0 - evaluate"] = @benchmarkable evaluate($p, $dom)
approx = evaluate(p, dom)

relprec = relative_precision(approx, ref)
push!(RP, relprec)
ref = Interval(0.3662865, 27.7256)

p = (3 * u*u*u - 6 * u*u + 4) / 6.0

SUITE["Daisy"]["bspline1 - evaluate"] = @benchmarkable evaluate($p, $dom)
approx = evaluate(p, dom)
ref = Interval(-65.13665, 0.5630625)
relprec = relative_precision(approx, ref)
push!(RP, relprec)

p = (-3 * u*u*u  + 3*u*u + 3*u + 1) / 6.0

SUITE["Daisy"]["bspline2 - evaluate"] = @benchmarkable evaluate($p, $dom)
approx = evaluate(p, dom)
ref = Interval(0.07399765, 53.59615)
relprec = relative_precision(approx, ref)
push!(RP, relprec)
#=
p = -u*u*u / 6.0

SUITE["Daisy"]["bspline3 - evaluate"] = @benchmarkable evaluate($p, $dom)
approx = evaluate(p, dom)
ref = Interval( )
relprec = relative_precision(approx, ref)
push!(RP, relprec)
=#
# ==========
# Doppler
# ==========

dom = a×b×c

v, u, T = set_variables(Float64,["v","u","T"],order = 4)

p = (- ((331.4 + 0.6 * T)) *v) /
                            (((331.4 + 0.6 * T) + u)*((331.4 + 0.6 * T) + u))

SUITE["Daisy"]["Doppler - evaluate"] = @benchmarkable evaluate($p, $dom)
approx = evaluate(p, dom)
ref = Interval(0.000887536, 0.0134537)
relprec = relative_precision(approx, ref)
push!(RP, relprec)

# ==========
# himmilbeau
# ==========

dom = a×b
x1, x2 = set_variables(Float64,["x1","x2"],order = 5)

p = (x1*x1 + x2 - 11)*(x1 * x1 + x2 - 11) + (x1 + x2*x2 - 7)*(x1 + x2*x2 - 7)

SUITE["Daisy"]["himmilbeau - evaluate"] = @benchmarkable evaluate($p, $dom)
approx = evaluate(p, dom)
ref = Interval(-2.58927, 110.461)
relprec = relative_precision(approx, ref)
push!(RP, relprec)

# ==========
# kepler
# ==========

dom = a×b×c×d×e×f
x1, x2, x3, x4, x5, x6 = set_variables(Float64,["x1","x2","x3",
                                                "x4","x5","x6"],order = 3)
p =  x2 * x5 + x3 * x6 - x2 * x3 - x5 * x6
                                   + x1 * (-x1 + x2 + x3 - x4 + x5 + x6)

SUITE["Daisy"]["kepler0 - evaluate"] = @benchmarkable evaluate($p, $dom)
approx = evaluate(p, dom)
ref = Interval(-11.8898, 28.6359)
relprec = relative_precision(approx, ref)
push!(RP, relprec)

dom = a×b×c×d
x1, x2, x3, x4 = set_variables(Float64,["x1","x2","x3","x4"],order = 3)
p = x1 * x4 * (-x1 + x2 + x3 - x4) + x2 * (x1 - x2 + x3 + x4)
    + x3 * (x1 + x2 - x3 + x4) -x2 * x3 * x4 - x1 * x3 - x1 * x2 - x4

SUITE["Daisy"]["kepler1 - evaluate"] = @benchmarkable evaluate($p, $dom)
approx = evaluate(p, dom)
ref = Interval(-240.666, 162.141)
relprec = relative_precision(approx, ref)
push!(RP, relprec)

dom = a×b×c×d×e×f
x1, x2, x3, x4, x5, x6 = set_variables(Float64,["x1","x2","x3",
                                                "x4","x5","x6"],order = 3)
p =  x1 * x4 * (-x1 + x2 + x3 - x4 + x5 + x6)
    + x2 * x5 * (x1 - x2 + x3 + x4 - x5 + x6) +x3* x6 * (x1 + x2 - x3 + x4 + x5 - x6)
    - x2 * x3 * x4 -x1* x3* x5 - x1 * x2 * x6 - x4 * x5 * x6

SUITE["Daisy"]["kepler2 - evaluate"] = @benchmarkable evaluate($p, $dom)
approx = evaluate(p, dom)
ref = Interval(-604.201, 468.374)
relprec = relative_precision(approx, ref)
push!(RP, relprec)

# =========
# Rigidbody
# =========

dom = a×b×c
x1, x2, x3 = set_variables(Float64,["x1","x2","x3"], order = 3)

p = -x1*x2 - 2*x2*x3 - x1 - x3

SUITE["Daisy"]["Rigidbody1 - evaluate"] = @benchmarkable evaluate($p, $dom)
approx = evaluate(p, dom)
ref = Interval(-21.2701, -0.539999)
relprec = relative_precision(approx, ref)
push!(RP, relprec)

p = 2*(x1*x2*x3) + (3*x3*x3) - x2*(x1*x2*x3) + (3*x3*x3) - x2

SUITE["Daisy"]["Rigidbody2 - evaluate"] = @benchmarkable evaluate($p, $dom)
approx = evaluate(p, dom)
ref = Interval(54.9599, 362.769)
relprec = relative_precision(approx, ref)
push!(RP, relprec)
