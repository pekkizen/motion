package motion

import (
	"fmt"
	"math"
)

const (
	forceError         = false
	useWind            = true
	useFMA             = true
	countFunctionCalls = true
	minBracketLen      = 0.01
	minTolNR           = 1e-10
	maxVel             = 500.0
	useExactGrade      = false // if true 1/math.Sqrt(1+tan*tan) is used for cos
)

const (
	NewtonRaphsonMethod = 1
	NewtonHalleyMethod  = 2
	Householder3Method  = 3
	singleQuadratic     = 4
	doubleQuadratic     = 5
	doubleLinear        = 6
	singleLinear        = 7
)

// velFromPower takes the actual velocity solver function uses in func VelFromPower.
var velFromPower func(*BikeCalc, float64, float64, float64) (float64, int)

type store struct {
	tan   float64
	power float64
	wind  float64
}

type BikeCalc struct {
	store         store
	gravity       float64
	temperature   float64 // deg Celsius
	airPressure   float64 // Pascals
	rho           float64 // air density
	baseElevation float64 // meters
	wind          float64 // m/s
	power         float64 // Watts
	speed         float64 // m/s

	cd          float64 // clothing and bike slippery coefficient. This is difficult to measure
	frontalArea float64 // rider and bike and banniers, m^2
	cdA         float64 // cd * frontalArea. Primary parameter. Mostly used as such.
	crr         float64 // rolling resistance coefficient
	cbf         float64 // braking friction coefficient, G-force at level gound.
	ccf         float64 // cornering friction coefficient, for cornering speed limits.
	cDrag       float64 // 0.5 * cdA * rho
	oCDrag      float64 // 1/cDrag

	// weight dependant
	mass         float64 // total mass
	massRotating float64 // rotating mass ~ rims + tires + tubes
	mg           float64 // mass * gravity
	mgCrr        float64 // mg * crr
	mgCbf        float64 // mg * cbf
	massKin      float64 // mass + 0.9*rotatingMass
	oMassKin     float64 // 1/massKin

	// road slope dependant
	tan float64 // grade % / 100
	cos float64 // 1/math.Sqrt(1+tan*tan)
	sin float64 // tan * cos

	// road slope and weight dependant forces
	fGrav  float64 // sin * mg
	fRoll  float64 // cos * mgCrr
	fBrake float64 // cos * mgCbf
	fGR    float64 // fGrav + fRoll

	solverFunc  int
	minPower    float64
	bracketLen  float64
	tolVel      float64
	iter        int
	maxIter     int
	callsf      int
	callsSolver int

	velErr        float64
	velErrAbs     float64
	velErrMax     float64
	velErrPos     int
	callsErr      int
	calcVelErrors bool
	counter       int // event couter for testing

	prevVel  float64
	prevSin  float64
	prevWind float64
	// prevDrag  float64
	// prevForce float64
	errmsg []byte
}

// Calculator returns a new calculator. Parameters must be set
// before using the calculator. The default values are NaN.
func Calculator() *BikeCalc {
	NaN := math.NaN()

	c := &BikeCalc{
		gravity:       9.80665,  //Standard sea level value (for lat 45.5 deg)
		rho:           1.2250,   //Standard sea level value at 15 deg C
		temperature:   15,       //Standard deg Celsius.
		airPressure:   101325.0, //standard sea level atmospheric pressure (Pascals)
		baseElevation: 0,
		wind:          0,

		cbf:          NaN,
		cdA:          NaN,
		crr:          NaN,
		ccf:          NaN,
		cDrag:        NaN,
		mass:         NaN,
		massRotating: NaN,

		bracketLen:    -1,
		solverFunc:    NewtonRaphsonMethod,
		tolVel:        -1,
		calcVelErrors: false,
		minPower:      0.01,
		prevVel:       4,
	}
	c.SetWeight(0)
	c.SetGrade(0)
	c.SetVelSolver(c.solverFunc)
	return c
}

func (c *BikeCalc) Error() string {
	if c.errmsg == nil {
		return ""
	}
	return string(c.errmsg) //+ fmt.Sprintf("\n\n%#v", c)
}
func (c *BikeCalc) Dump() string {

	return "\n\nCalculator struct dump: " + fmt.Sprintf("\nn%#v", c)
}
func (c *BikeCalc) appendErr(s string) {
	c.errmsg = append(c.errmsg, s...)
}

func (c *BikeCalc) storeState() {
	c.store.tan = c.tan
	c.store.power = c.power
	c.store.wind = c.wind
}
func (c *BikeCalc) restoreState() {
	c.SetGradeExact(c.store.tan)
	c.power = c.store.power
	c.wind = c.store.wind
}

// signSq returns (x + wind) * math.Abs(x + wind).
func (c *BikeCalc) signSq(v float64) float64 {
	if !useWind {
		return v * v
	}
	v += c.wind
	return v * math.Float64frombits(math.Float64bits(v)&^(1<<63))
	// This compiles to: (go x86-64 1.22), https://godbolt.org/z/cc97eEso1
	// ADDSD   X1, X0
	// MOVQ    X0, AX
	// BTRQ    $63, AX
	// MOVQ    AX, X1
	// MULSD   X1, X0
	// RET
}

/*
Function NewtonRaphson solves v from the equation
	v x (fGR + cDrag x (v+wind) x |v+wind|) - power = 0
For speed for negative power there are normally two positive roots on downhill
slopes. The left root for heavyvbraking and low speed and the right root for
high speed and light braking. For too much negative power there is no root.
In some weird conditions: low negative power, tail wind and not very steep
downhill slope, there can be four roots around zero air speed. Asking speed
for a negative power is generally not a proper idea. On the contrary, there
is a unique power for any speed for any road slope.

For positive power there are two roots, of which the left one
is negative. Between the roots there is a derivative zero point.
The function is 3. degree polynomial, but changing air drag force
direction at air speed 0, makes is look like/behave as
2. degree polynomial with two roots. With small power, negative wind and fGR = 0
there is a double root at -wind = v, which is also derivative zero point.
This root is solvable by NewtonRaphson with up to ~50 iterations.
fRG is 0 when road slope tangent is -crr.
*/

// NewtonRaphson returns speed for power and number of iterations used.
// For power >= 0 the function is globally convergent to the rightmost root.
// Iteration is stopped when change in speed is less than tol.
// v is initial speed for iteration. NewtonRaphson can solve freewheeling speeds
// for power = 0, but for this the function VelFreewheeling is faster with
// full accuracy. For power < 0 NewtonRaphson returns a root, but this
// may not be the rightmost one.
// Max abs(error) ~ tol^2 x 0.03, mean abs(error) ~ tol^2 x 0.003.
func (c *BikeCalc) NewtonRaphson(power, tol, v float64) (vres float64, iter int) {
	if tol <= 0 {
		tol = minTolNR
	}
	// for bias adjustment, otherwise every
	// root is biased to the right.
	const bias = 0.105
	var (
		stepright = 4.0
		vAir      = v + c.wind
		cDrag     = c.cDrag
	)
	for n := 1; n < 100; n++ {
		if cDrag*vAir < 0 {
			cDrag = -cDrag
		}
		Δv, der := c.newtonIter(power, v, vAir, cDrag)

		if der > 0 && math.Abs(Δv) < tol {
			c.iter += n
			return (v + Δv) - Δv*Δv*bias, n
		}
		if der <= 0 || Δv > stepright {
			Δv = stepright
			stepright *= 1.05 // for very rare ping pong cases
		}
		v += Δv
		vAir += Δv
	}
	if len(c.errmsg) < 60 {
		c.appendErr(" Newton-Raphson did not converge: ")
	}
	c.speed = v
	return v, 0
}

func (c *BikeCalc) newtonIter(power, v, vAir, cDrag float64) (Δv, der float64) {
	if useFMA {
		der = math.FMA(cDrag*vAir, math.FMA(v, 2, vAir), c.fGR)
		Δv = -math.FMA(v, math.FMA(cDrag*vAir, vAir, c.fGR), -power) / der
		return
	}
	der = c.fGR + (cDrag*vAir)*(vAir+2*v)
	Δv = -(-power + v*c.fGR + (cDrag*vAir)*(vAir*v)) / der
	return
}

// TraubOst returns speed for power and number of iterations used.
// TraubOst is deprecated, because it is slower in any skenario.
func (c *BikeCalc) TraubOst(power, tol, v float64) (vResult float64, iter int) {
	if tol <= 0 {
		tol = minTolNR
	}
	// v = -1
	var (
		stepright = 4.0
		vAir      = v + c.wind
		cDrag     = c.cDrag
		der, Δv   float64
	)

	for n := 1; n < 100; n++ {
		if cDrag*vAir < 0 {
			cDrag = -cDrag
			// c.counter++
		}
		// if useFMA {
		// 	der = math.FMA(cDrag*vAir, math.FMA(v, 2, vAir), c.fGR)
		// 	f := math.FMA(v*cDrag, (vAir * vAir), math.FMA(v, c.fGR, -power))
		// 	Δv = -f / der
		// 	fw := -power + (v+Δv)*math.FMA(cDrag, c.signSq(v+Δv), c.fGR)
		// 	Δv = Δv * (fw - f) / math.FMA(2, fw, -f)
		// } else {
		der = c.fGR + cDrag*vAir*(vAir+2*v)
		f := -power + v*c.fGR + cDrag*vAir*(vAir*v)
		Δv = -f / der
		xAir := v + Δv + c.wind
		fw := -power + (v+Δv)*(c.fGR+cDrag*xAir*xAir)
		Δv = Δv * (fw - f) / (2*fw - f)

		// }
		// Δv, der := c.TraubOstIter(power, v, vAir, cDrag)

		if aΔv := math.Abs(Δv); der > 0 && aΔv < tol {
			c.iter += n
			v += Δv
			vAir += Δv
			if aΔv < 0.6*tol { ///&& vAir > 0 {
				// c.counter++
				return v, n
			}
			if cDrag*vAir < 0 {
				cDrag = -cDrag
				// c.counter++
			}
			// c.counter++
			Δv, _ := c.newtonIter(power, v, vAir, cDrag)
			v += Δv
			return v, n
		}
		// if der <= 0 || math.Abs(Δv) > stepright {
		if der <= 0 || Δv > stepright {
			Δv = stepright
			stepright *= 1.1
			// c.counter++
		}
		v += Δv
		vAir += Δv
	}
	if len(c.errmsg) < 60 {
		c.appendErr("Newton-Halley did not converge")
	}
	c.speed = v
	return v, 0
}

func (c *BikeCalc) TraubOstIter(power, v, vAir, cDrag float64) (Δv, d float64) {

	d = c.fGR + cDrag*vAir*(vAir+2*v)
	f := -power + v*c.fGR + v*cDrag*(vAir*vAir)
	v -= f / d
	fw := -power + v*c.fGR + v*cDrag*c.signSq(v)
	Δv = -(fw - f) * f / ((2*fw - f) * d)
	return
}

/*
NewtonHalley returns speed for power and number of iterations used.
This has 1.5 x convergence but also higher calculation cost compared to
NewtonRaphson. To same accuracy, NewtonHalley is sligthly faster (~20%) than NewtonRaphson.
NewtonHalley returns natively quite unbiased speeds. Max abs(error) < tol^3 * 0.005.
Halley's method expects continuous 2nd derivative. 2nd derivative has a discontinuity
at vAir = 0, which makes its use unsafe (propably fails) for negative air speeds .
For negative air speeds NewtonRaphson is called. Still, NewtonHalley is unsafe for
appraching the root.
*/
func (c *BikeCalc) NewtonHalley(power, tol, v float64) (vResult float64, iter int) {
	if tol <= 0 {
		tol = minTolNR
	}
	var (
		stepright = 4.0
		vAir      = v + c.wind
		cDrag     = c.cDrag
	)
	for n := 1; n < 100; n++ {
		if vAir <= 0 {
			return c.NewtonRaphson(power, 0.09*tol, v)
		}
		Δv, der := c.newtonHalleyIter(power, v, vAir, cDrag)

		if aΔv := math.Abs(Δv); der > 0 && aΔv < tol {
			c.iter += n
			v += Δv
			vAir += Δv
			if aΔv < 0.5*tol && vAir > 0 {
				return v, n
			}
			if vAir < 0 {
				cDrag = -cDrag
			}
			Δv, _ := c.newtonIter(power, v, vAir, cDrag)
			return v + Δv, n + 1
		}
		if der <= 0 || Δv > stepright {
			Δv = stepright
			stepright *= 1.05
		}
		v += Δv
		vAir += Δv
	}
	if len(c.errmsg) < 60 {
		c.appendErr("Newton-Halley did not converge :")
	}
	c.speed = v
	return v, 0
}

func (c *BikeCalc) newtonHalleyIter(power, v, vAir, cDrag float64) (Δv, der float64) {
	if useFMA {
		der = math.FMA(cDrag*vAir, math.FMA(v, 2, vAir), c.fGR)
		fun := math.FMA(v, math.FMA(cDrag*vAir, vAir, c.fGR), -power)
		Δv = -der * fun / (der*der - fun*cDrag*math.FMA(2, vAir, v))
		return
	}
	der = c.fGR + cDrag*vAir*(vAir+2*v)
	fun := -power + v*c.fGR + cDrag*vAir*(vAir*v)
	Δv = -der * fun / (der*der - fun*cDrag*(v+2*vAir))
	return
}

// Householder3 returns speed for power and number of iterations used.
// https://en.wikipedia.org/wiki/Householder%27s_method
// This has very much the same computational time to same accuracy than NewtonHalley.
// Also unsafe for negative air speeds, because of the same 2nd derivative.
// For negative air speeds NewtonRaphson is called.
// With not very much wind, Householder3 is 30% faster than NewtonRaphon.
func (c *BikeCalc) Householder3(power, tol, v float64) (vresult float64, iter int) {
	if tol <= 0 {
		tol = minTolNR * 5
	}
	if v < 1 {
		v = 1
	}
	const stepright = 4.0
	var (
		vAir  = v + c.wind
		cDrag = c.cDrag
	)
	for n := 1; n < 100; n++ {
		if vAir <= 0 {
			return c.NewtonRaphson(power, 0.07*tol, v)
		}
		Δv, der := c.householderIter(power, v, vAir, cDrag)

		if aΔv := math.Abs(Δv); der > 0 && aΔv < tol {
			c.iter += n
			v += Δv
			vAir += Δv

			if aΔv < 0.6*tol && vAir > 0 {
				return v, n
			}
			if vAir < 0 {
				cDrag = -cDrag
			}
			Δv, _ := c.newtonIter(power, v, vAir, cDrag)
			return v + Δv, n + 1
		}
		if der <= 0 || Δv > stepright {
			Δv = stepright
		}
		v += Δv
		vAir += Δv
	}
	if len(c.errmsg) < 60 {
		c.appendErr("Householder3 did not converge :")
	}
	c.speed = v
	return v, 0
}

func (c *BikeCalc) householderIter(power, v, vAir, cDrag float64) (Δv, d float64) {
	if useFMA {
		f := math.FMA(v, math.FMA(cDrag*vAir, vAir, c.fGR), -power)
		d = math.FMA(cDrag*vAir, math.FMA(v, 2, vAir), c.fGR)
		d2 := 2 * cDrag * math.FMA(2, vAir, v)
		Δv = -f * math.FMA(-0.5*f, d2, d*d) /
			math.FMA(d, math.FMA(-f, d2, d*d), f*f*cDrag)
		return
	}
	f := -power + v*c.fGR + cDrag*vAir*(vAir*v)
	d = c.fGR + cDrag*vAir*(vAir+2*v)
	d2 := 2 * cDrag * (v + 2*vAir)
	Δv = -f * (d*d - 0.5*f*d2) / (d*(d*d-f*d2) + f*f*cDrag) // d3 = 6 * cDrag
	return
}

// VelFreewheel returns speed (>= 0) for power = 0 and any wind.
// The speed v is solved from: fGR + cDrag x (v + wind)^2 = 0.
func (c *BikeCalc) VelFreewheel() (v float64) {

	v = math.Sqrt(math.Abs(c.fGR) * c.oCDrag)
	if c.fGR > 0 {
		v = -v
	}
	// v = max(1, -c.wind + math.Copysign(math.Sqrt(math.Abs(c.fGR)/c.cDrag), -c.fGR)
	if v -= c.wind; v > 0 {
		return
	}
	return 0
}

// VelError calculates  error in speed for given power.
// The error is calculated as difference between speed from
// max accuracy NewtonRaphson.
func (c *BikeCalc) VelError(power, v float64) (err, absErr float64) {
	if power < 0 {
		return
	}
	vExact, iter := c.NewtonRaphson(power, minTolNR, v)
	c.iter -= iter
	if iter <= 0 {
		return
	}
	c.callsErr++
	err = (v - vExact) // vExact
	absErr = math.Abs(err)
	if err >= 0 {
		c.velErrPos++
	}
	c.velErr += err
	c.velErrAbs += absErr
	if absErr > c.velErrMax {
		c.velErrMax = absErr
	}
	return
}

// velGuess returns speed estimate by adjusting previous speed by road gradient sin change.
// This works well if we are calculating consecutive virtual road segments, eg. GPX route data.
func (c *BikeCalc) velGuess() (v float64) {

	v = c.prevVel * (1 + 10*(c.prevSin-c.sin) + 0.07*(c.prevWind-c.wind))
	if v > 1 {
		return
	}
	return 1
}

func (c *BikeCalc) storePrev(v float64) {
	c.prevVel = v
	c.prevSin = c.sin
	c.prevWind = c.wind
	// c.prevDrag = c.Fdrag(v)
	// c.prevForce = c.fGR + c.Fdrag(v)
}

// VelFromPower returns speed for power using initial speed velguess.
// For power < 0 VelFromPower returns 0, false. For power < minPower (c.SetMinPower)
// VelFreewheeling is returned. As a velocity solver function is used the one
// set by SetVelSover function. NewtonRaphson by default.
func (c *BikeCalc) VelFromPower(power, velguess float64) (v float64, ok bool) {
	if power < 0 {
		c.appendErr("VelFromPower called with power < 0")
		return 0, false
	}
	if power <= c.minPower {
		return c.VelFreewheel(), true
	}
	if forceError {
		c.cDrag = 0
	}
	c.power = power
	if velguess < 0 {
		velguess = c.velGuess()
	}
	c.callsSolver++

	v, iter := velFromPower(c, power, c.tolVel, velguess)

	if iter > c.maxIter {
		c.maxIter = iter
	}
	if c.errmsg != nil { //|| v < 0 || v > maxVel {
		return v, false
	}
	if c.calcVelErrors {
		c.VelError(power, v)
	}
	c.storePrev(v)
	return v, true
}

func (c *BikeCalc) DerivativeRoot() (root float64) {
	w := c.wind
	f := c.fGR * c.oCDrag
	if f <= 0 {
		// for c.fGR + c.cDrag*vAir*(vAir+2*v) = 0, vAir > 0
		return (-2*w + math.Sqrt(w*w-3*f)) * (1.0 / 3)
	}
	// for c.fGR - c.cDrag*vAir*(vAir+2*v) = 0, vAir < 0
	return (-2*w - math.Sqrt(w*w+3*f)) * (1.0 / 3)
}

// PowerFromVel returns power for speed v with wind.
func (c *BikeCalc) PowerFromVel(v float64) (power float64) {
	if !useWind {
		return v * (c.fGR + c.cDrag*v*v)
	}
	return v * (c.fGR + c.cDrag*c.signSq(v))
}

// FlatPower returns power for speed v on flat (grade = 0) road
// for wind = 0.
func (c *BikeCalc) FlatPower(v float64) (power float64) {
	if v <= 0 {
		return 0
	}
	return v * (c.mgCrr + c.cDrag*v*v)
}

// FlatSpeed returns speed for power on flat (grade = 0) road for wind = 0.
func (c *BikeCalc) FlatSpeed(power float64) (v float64) {
	if power <= 0 {
		return 0
	}
	c.storeState()
	c.SetGradeExact(0)
	c.SetWind(0)
	v, _ = c.NewtonRaphson(power, minTolNR, 6)
	c.restoreState()
	return
}

// GradeFromVelAndPower returns slope tangent for speed and power
// and wind = 0. The tangent tan is calculated from sin solved from:
// sin*mg + sqrt(1-sin^2)*mg*crr + cDrag*v^2 - power/v = 0
func (c *BikeCalc) GradeFromVelAndPower(v, power float64) (tan float64) {

	sin := (power/v - c.cDrag*v*v) / c.mg // sin for crr = 0
	r := c.crr * c.crr
	d := r * (1 + r - sin*sin)
	if d < 0 || v <= 0 {
		return math.NaN()
	}
	sin = (sin - math.Sqrt(d)) / (1 + r)
	tan = sin / math.Sqrt(1-sin*sin)
	return
}

// PowerFromVerticalUp returns power from vertical up speed (m/h).
// The power is calculated by riding uphill slope with tangent grade.
func (c *BikeCalc) PowerFromVerticalUp(v, grade float64) (power float64) {
	if v <= 0 || grade <= 0 {
		return math.NaN()
	}
	c.storeState()
	c.SetGradeExact(grade)
	v /= (c.sin * 3600) //speed on road m/s
	power = v * (c.fGR + c.cDrag*v*v)
	c.restoreState()
	return
}

// VerticalUpFromPower returns vertical speed up (m/h) for power.
// Speed is calculated by riding uphill slope with tangent grade.
func (c *BikeCalc) VerticalUpFromPower(power, grade float64) (vUp float64) {
	if power <= 0 || grade <= 0 {
		return math.NaN()
	}
	c.storeState()
	c.SetGradeExact(grade)
	c.SetWind(0)
	vUp, _ = c.NewtonRaphson(power, minTolNR, 4)
	vUp *= c.sin * 3600 //speed up m/h
	c.restoreState()
	return
}

// CdAfromVelAndPower returns CdA from speed and power (and Crr)
// on flat ground road and no wind.
func (c *BikeCalc) CdAfromVelAndPower(v, power float64) (CdA float64) {
	if power <= 0 || v <= 0 {
		return math.NaN()
	}
	CdA = (power - v*c.mgCrr) / (0.5 * c.rho * v * v * v) // can be < 0
	return
}

// CrrFromVelAndPower returns Crr from speed and power (and CdA)
// on flat ground road.
func (c *BikeCalc) CrrFromVelAndPower(v, power float64) (Crr float64) {
	if power <= 0 || v <= 0 {
		return math.NaN()
	}
	Crr = (c.Fdrag(v)*v - power) / (-v * c.mg)
	return
}

// VelFromTurnRadius returns approx. max cornering speed for turn radius and
// cornering friction coefficient ccf. The cornering speed v solved from
// m x v^2 / radius = m x gravity x ccf
// This is very simplified version of the "true" function, which should
// model downhill turn front wheel slip at least ect. In steeper downhill
// slopes also downhill gravity force should be taken into account.
func (c *BikeCalc) VelFromTurnRadius(radius float64) float64 {

	return math.Sqrt(radius * c.gravity * c.ccf)
}

/*
Latitude equation for gravity
https://en.wikipedia.org/wiki/Gravity_of_Earth

Elevation correction for gravity
https://www.geology.cwu.edu/facstaff/tim/TEACHING/Geophysics/gravity_geoid.pdf
The free air correction amounts to -3.1 x 10-6 m /s2 per meter of elevation.
"assuming a crustal density of 2.7 g/cm3, the Bouguer correction
is 1.1 x 10-6 m/s2 per meter of elevation."
*/

// LocalGravity returns local gravitational acceleration for latitude and elevation.
// Formula for gravity as a function of latitude is the WGS (World Geodetic System) 84
// Ellipsoidal Gravity Formula
func (c *BikeCalc) LocalGravity(degLatitude, metersEle float64) (gravity float64) {
	const eleCorrection = (-3.1 + 1.1) * 1e-6 // free air correction + Bouguer correction

	sin := math.Sin(degLatitude * (math.Pi / 180))
	sin *= sin
	gravity = 9.78032534 * (1 + 0.001931853*sin) / math.Sqrt(1-0.00669438*sin)
	gravity += eleCorrection * metersEle
	return
}

/*
	https://en.wikipedia.org/wiki/Density_of_air
	https://en.wikipedia.org/wiki/Barometric_formula#Density_equations

RhoFromEle calculates air density at a given elevation based on the base
elevation and its corresponding temperature and air pressure values.
Use functions

	(c *BikeCalc).SetBaseElevation(x float64)
	(c *BikeCalc).SetTemperature(x float64)
	(c *BikeCalc).SetAirPressure(x float64)
	(c *BikeCalc) LocalGravity(degLatitude, baseElevation) // minor effect on rho

to set parameters. At least the proper temperature should be given.
Default values for parameters

	gravity:       9.80665  //Standard sea level value (for lat 45.5 deg)
	rho:           1.2250   //Standard sea level value at 15 deg C
	temperature:   15       //Standard deg Celsius.
	airPressure:   101325.0 //standard sea level atmospheric pressure, Pascals
	baseElevation: 0

The function RhoFromEle assumes that the air is dry. The presence
of water vapor in the air reduces its density slightly, as water vapor
is lighter than dry air.
*/
func (c *BikeCalc) RhoFromEle(elevation float64) (rho, temperature float64) {
	const (
		M = 0.0289644 //molar mass of Earth's air (kg/mol)
		R = 8.3144598 //universal gas constant (N·m/(mol·K))
		L = -0.0065   //temperature lapse rate (per 1 m up)
	)
	var (
		P    = c.airPressure          // at baseElevation
		G    = c.gravity              // at baseElevation
		T    = c.temperature + 273.15 // at baseElevation
		dEle = elevation - c.baseElevation
	)
	temperature = c.temperature + dEle*L

	// Base elevation rho at pressure P and temperature T
	rho = P * M / (R * T)
	if dEle == 0 {
		return
	}
	// Elevation adjustment from Base elevation air density
	// rho to dEle meters above or under base elevation
	x := dEle * L / T
	y := -(1 + G*M/(R*L))
	rho *= math.Pow(1+x, y)
	return
}
