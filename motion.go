package motion

import (
	"math"
)

const (
	forceError    = false
	useWind       = true
	minBracketLen = 0.05
	minTolNR      = 1e-12
	maxVel        = 500.0
	useExactGrade = false // use 1/math.Sqrt(1+tan*tan) for cos
)
const (
	newtonRaphson   = 1
	newtonHalley    = 2
	singleQuadratic = 3
	doubleQuadratic = 4
	doubleLinear    = 5
)

// velFromPower keeps the actual velocity solver function.
var velFromPower func(*BikeCalc, float64, float64, float64) (float64, int)

type store struct {
	tan   float64
	power float64
	wind  float64
}

type BikeCalc struct {
	errmsg []byte
	store  store

	gravity       float64
	temperature   float64
	airPressure   float64
	rho           float64
	baseElevation float64
	wind          float64
	power         float64

	cd          float64
	frontalArea float64
	cdA         float64 // cd * frontalArea
	crr         float64 // rolling resistance coefficient
	cbf         float64 // braking friction
	ccf         float64 // cornering friction coefficient, for cornering speeds
	cDrag       float64 // 0.5 * cdA * rho

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
	iterNR      int
	maxIter     int
	callsf      int
	callsSolver int

	velErr        float64
	velErrAbs     float64
	velErrMax     float64
	velErrPos     int
	callsErr      int
	calcVelErrors bool

	prevVel float64
	prevTan float64
}

// Calculator returns a new calculator. Riding parameters must be set
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
		massRotating:  0,

		cbf:   NaN,
		cdA:   NaN,
		crr:   NaN,
		ccf:   NaN,
		cDrag: NaN,
		mass:  NaN,

		bracketLen:    0.25,
		solverFunc:    newtonRaphson,
		tolVel:        minTolNR,
		calcVelErrors: false,
		minPower:      0.01,
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
	return string(c.errmsg)
}
func (c *BikeCalc) appendErr(s string) {
	s = ": " + s
	c.errmsg = append(c.errmsg, s...)
}

func (c *BikeCalc) storeState() {
	c.store.tan = c.tan
	c.store.power = c.power
	c.store.wind = c.wind
}
func (c *BikeCalc) restoreState() {
	c.SetGrade(c.store.tan)
	c.power = c.store.power
	c.wind = c.store.wind
}

// signSq returns (x + c.wind) * math.Abs(x + c.wind).
func (c *BikeCalc) signSq(x float64) float64 {
	if !useWind {
		return x * x
	}
	x += c.wind
	return x * math.Float64frombits(math.Float64bits(x)&^(1<<63)) // inline cost 20
	// return x * math.Abs(x) // same but inline cost 24
}

/*
Function NewtonRaphson solves v from the equation
	v (fGR + cDrag (v+wind) |v+wind|) - power = 0
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
2. degree polynomial with two roots. With zero power, negative wind and fGR = 0
there is a double root at -wind = v, which is also derivative zero point.
This root is solvable by NewtonRaphson with up to ~50 iterations.
fRG is 0 when road slope tangent = -crr.
*/

// NewtonRaphson returns speed for power and number of iterations used.
// For power >= 0 the function is globally convergent to the rightmost root.
// Iteration is stopped when change in speed is less than tol.
// v is initial speed for iteration. NewtonRaphson can solve freewheeling speed
// for power = 0, but for this the function VelFreewheeling is faster with
// full accuracy. For power < 0 NewtonRaphson returns -1, 0.
// Max abs(error) < tol^2 x 0.03, mean abs(error) < tol^2 x 0.003.
func (c *BikeCalc) NewtonRaphson(power, tol, v float64) (vel float64, iter int) {
	if power < 0 {
		return -1, 0
	}
	if tol < minTolNR {
		tol = minTolNR
	}
	tol *= tol
	const (
		stepright = 4
		bias      = 0.065
	)
	vAir := v + c.wind
	for n := range 100 {
		dva := c.cDrag * vAir
		if dva < 0 {
			dva = -dva
		}
		var (
			fSum  = c.fGR + dva*vAir
			deriv = fSum + dva*v*2
			Δv    = (power - v*fSum) / deriv
		)
		if Δv*Δv < tol && deriv > 0 {
			c.iterNR += n + 1
			Δv *= 1 - Δv*bias //bias adjustment
			return v + Δv, n
		}
		if deriv <= 0 || Δv > stepright {
			Δv = stepright
		}
		v += Δv
		vAir += Δv
	}
	// Never here for power >= 0 and tol >= 1e-12 and cDrag > 0.
	// Approaching double root may take near 50 iterations.
	if len(c.errmsg) < 60 {
		c.appendErr("Newton-Raphson did not converge")
	}
	return v, 0
}

// NewtonHalley returns speed for power and number of iterations used.
// This has 1.5 x convergence but also 1.5 x calculation cost compared to
// NewtonRaphson. You can use 5 x larger tol than with NewtonRaphson, but
// execution times are about the same.  NewtonHalley returns natively
// quite unbiased speeds. Max abs(error) < tol^3 * 0.005.
func (c *BikeCalc) NewtonHalley(power, tol, v float64) (vel float64, iter int) {
	if power < 0 {
		return -1, 0
	}
	if tol < minTolNR {
		tol = minTolNR
	}
	const stepright = 4.0
	vAir := v + c.wind
	for n := range 100 {
		cDrag := c.cDrag
		if vAir < 0 {
			cDrag = -cDrag
		}
		var (
			fSum  = c.fGR + (cDrag*vAir)*vAir
			deriv = fSum + (cDrag*vAir)*v*2
			Δv    = (power - v*fSum) / deriv
		)
		switch {
		case deriv <= 0 || Δv > stepright:
			Δv = stepright
		case vAir*(vAir+Δv) <= 0:
			// Halley's method expects continuous 2nd derivative.
			// 2nd derivative has discontinuity at vAir = 0.
			// NewtonRaphson is used near it.
			c.iterNR += n + 1
			return c.NewtonRaphson(power, tol*0.2, v+Δv)
		default:
			Δv *= deriv / (deriv + Δv*cDrag*(v+2*vAir)) //Halley's booster for Newton's method
			if math.Abs(Δv) < tol {
				c.iterNR += n + 1
				return v + Δv, n
			}
		}
		v += Δv
		vAir += Δv
	}
	if len(c.errmsg) < 60 {
		c.appendErr("Newton-Halley did not converge")
	}
	return v, 0
}

// VelError calculates  error in speed for given power.
// The error is calculated as difference between speed from
// max accuracy NewtonRaphson.
func (c *BikeCalc) VelError(power, v float64) {
	if power <= c.minPower { //why not 0?
		return
	}
	vExact, iter := c.NewtonRaphson(power, minTolNR, v)
	c.iterNR -= iter
	if iter == 0 {
		return
	}
	c.callsErr++
	err := (v - vExact) / vExact
	absErr := math.Abs(err)
	if err >= 0 {
		c.velErrPos++
	}
	c.velErr += err
	c.velErrAbs += absErr
	if absErr > c.velErrMax {
		c.velErrMax = absErr
	}
}

// velGuess returns speed estimate by adjusting previous speed by road gradient change.
// For a case we are calculating consecutive virtual road segments, eg. GPX route data.
func (c *BikeCalc) velGuess() (v float64) {
	v = c.prevVel * (1 + 11*(c.prevTan-c.tan))
	if v < 1 {
		v = 1
	}
	return
}

func (c *BikeCalc) storePrevVel(v float64) {
	c.prevVel = v
	c.prevTan = c.tan
}

// VelFromPower returns speed for power using initial speed velguess.
// For power < 0 VelFromPower returns 0, false. For power < minPower (c.SetMinPower)
// VelFreewheeling is returned. As a velocity solver function is used the one
// set by SetVelSover function. NewtonRaphson by default.
func (c *BikeCalc) VelFromPower(power, velguess float64) (vel float64, ok bool) {
	if power < 0 {
		c.appendErr("VelFromPower called with power < 0")
		return 0, false
	}
	if power <= c.minPower {
		vel, ok = c.VelFreewheeling(), true
		return
	}
	if forceError {
		c.cDrag = 0
	}
	if velguess < 0 {
		velguess = c.velGuess()
	}
	c.callsSolver++

	vel, iter := velFromPower(c, power, c.tolVel, velguess)

	if iter > c.maxIter {
		c.maxIter = iter
	}
	if c.errmsg != nil {
		return vel, false
	}
	if c.calcVelErrors {
		c.VelError(power, vel)
	}
	c.storePrevVel(vel)
	return vel, true
}

// PowerFromVel returns power for speed v with wind.
func (c *BikeCalc) PowerFromVel(v float64) (power float64) {
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
func (c *BikeCalc) FlatSpeed(power float64) (vel float64) {
	if power <= 0 {
		return 0
	}
	c.storeState()
	c.SetGrade(0)
	c.wind = 0
	vel, _ = c.NewtonRaphson(power, minTolNR, 6)
	c.restoreState()
	return
}

// GradeFromVelAndPower returns slope tangent for speed and power
// and wind = 0. The tangent tan is calculated from sin solved from:
// sin*mg + sqrt(1-sin^2)*mg*crr + cDrag*v^2 - power/v = 0
func (c *BikeCalc) GradeFromVelAndPower(v, power float64) (tan float64) {
	if v <= 0 {
		return math.NaN()
	}
	sin := (power/v - c.cDrag*v*v) / c.mg //sin for crr = 0
	if math.Abs(sin) > 1 {
		return math.NaN()
	}
	r := c.crr * c.crr
	d := r * (1 + r - sin*sin)
	sin = (sin - math.Sqrt(d)) / (1 + r)
	tan = sin / math.Sqrt(1-sin*sin)
	return
}

// VelFreewheeling returns speed (>= 0) for power = 0 and any wind.
// The speed is v solved from: fGR + cDrag x (v + wind)^2 = 0.
func (c *BikeCalc) VelFreewheeling() (vel float64) {

	vel = math.Sqrt(math.Abs(c.fGR / c.cDrag))
	if c.fGR > 0 {
		vel = -vel
	}
	vel -= c.wind
	if vel < 0 {
		return 0
	}
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
	c.wind = 0
	vUp, _ = c.NewtonRaphson(power, minTolNR, 4)
	vUp *= c.sin * 3600 //speed up m/h
	c.restoreState()
	return
}

// CdAfromVelAndPower returns CdA from speed and power on flat ground road
// and no wind.
func (c *BikeCalc) CdAfromVelAndPower(v, power float64) (CdA float64) {
	if power <= 0 || v <= 0 {
		return math.NaN()
	}
	CdA = (power - v*c.mgCrr) / (0.5 * c.rho * v * v * v) // can be < 0
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
	const eleCorrection = (-3.1 + 1.1) * 1e-6 // -free air correction + Bouguer correction

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
		P    = c.airPressure
		G    = c.gravity
		T    = c.temperature + 273.15
		dEle = elevation - c.baseElevation
	)
	temperature = c.temperature + dEle*L

	// Base elevation rho at pressure P and temperature T
	rho = P * M / (R * T)
	if dEle == 0 {
		return
	}
	// Elevation adjustment from Base elevation air density
	// rho at dEle meters above or under base elevation
	x := dEle * L / T
	y := -(1 + G*M/(R*L))
	rho *= math.Pow(1+x, y)
	return
}
