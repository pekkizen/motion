/*
Package motion documentation goes here.
*/
package motion

import (
	"math"
)

const (
	forceERR        = false
	minAcceleration = 0.001
	minBracketLEN   = 1e-8
	minTOLNR        = 1e-12
	maxVEL          = 500.0
	maxVELS         = "500 m/s"
	noWIND          = false
)
const (
	newtonRaphSON   = 1
	singleQuadRATIC = 2
	doubleQuadRATIC = 3
	doubleLINEAR    = 4
)

// type solver func(*BikeCalc, float64, float64, float64) (float64, int)
// var velSolver solver

var velSolver func(*BikeCalc, float64, float64, float64) (float64, int)

type temporaryStore struct {
	tan   float64
	power float64
	wind  float64
}

type BikeCalc struct {
	errmsg []byte
	store  temporaryStore

	gravity       float64
	temperature   float64
	airPressure   float64
	rho           float64
	baseElevation float64
	wind          float64
	power         float64
	mass          float64
	massRotating  float64

	cd          float64
	frontalArea float64
	cdA         float64 // cd * frontalArea
	crr         float64
	cbf         float64
	cDrag       float64 // 0.5 * cdA * rho

	// weight dependant
	mg       float64 // mass * gravity
	mgCrr    float64 // mg * crr
	mgCbf    float64 // mg * cbf
	massKin  float64 // mass + 0.9*rotatingMass
	omassKin float64 // 1/massKin

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
	tolNR       float64
	minPower    float64
	bracketLen  float64
	iterNR      int
	maxIter     int
	callsF      int
	callsSolver int
	callsFW     int
	callsPFV    int

	velErr    float64
	velErrAbs float64
	velErrMax float64
	velErrPos int
	callsErr  int
	velErrors bool

	prevVel float64
	prevTan float64
}

// Calculator returns a new calculator.
func Calculator() *BikeCalc {
	c := &BikeCalc{
		gravity:       9.80665,  //Standard sea level value (for lat 45.5 deg)
		rho:           1.2250,   //Standard sea level value at 15 deg C
		temperature:   15,       //Standard deg Celsius.
		airPressure:   101325.0, //standard sea level atmospheric pressure (Pascals)
		baseElevation: 0,
		cbf:           0.25,
		cdA:           0.6,
		crr:           0.007,

		bracketLen: 0.5,
		solverFunc: newtonRaphSON,
		velErrors:  false,
		minPower:   1,
	}
	c.cDrag = 0.5 * c.cdA * c.rho
	c.mass = math.NaN()
	c.SetWeight(0)
	c.SetVelSolver(newtonRaphSON)
	return c
}

type errstr struct {
	s string
}



func (e *errstr) Error() string {
	return e.s
}
func (c *BikeCalc) Error() error {
	if len(c.errmsg) == 0 {
		return nil
	}
	return &errstr{string(c.errmsg)}
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
	c.SetGrade(c.store.tan)
	c.power = c.store.power
	c.wind = c.store.wind
}

func (c *BikeCalc) signSq(x float64) float64 {
	if noWIND {
		return x * x
	}
	x += c.wind
	return x * math.Float64frombits(math.Float64bits(x)&^(1<<63)) // inline cost 20
	// return  x * math.Abs(x) // same but inline cost 24
}

/*
Function NewtonRaphson
For speed for negative power there are normally two positive roots on downhill
slopes. The left root for heavyvbraking and low speed and the right root for
high speed and light braking. In some weird conditions: low negative power,
tail wind and not very steep downhill slope, there can be four roots around
zero air speed. Asking speed for a negative power is generally not a proper
idea. On the contrary, there is a unique power for any speed for any road slope.

For positive power there are two roots, of which the left one
is negative. Between the roots there is as a derivative zero point.
The function is 3. degree polynomial, but changing air drag force
direction at air speed 0, makes is look like/behave as
2. degree polynomial with two roots.

const bias is used to shift NR root to the left by dis
Otherwise ~all NR roots are on the right side of the real root.
*/

// NewtonRaphson returns speed for power and number of iterations used.
// For power > 0 the function is globally convergent to the rightmost root.
// Iteration is stopped when change in speed is less than tol. This gives
// max error < ~tol^2 x 0.1. v is initial speed for iteration.
// NewtonRaphson can solve freewheeling speed for power = 0, but for this
// the function VelFreewheeling is faster with full accuracy.
// For negative power NewtonRaphson returns 0, 0.
func (c *BikeCalc) NewtonRaphson(power, tol, v float64) (vel float64, iter int) {
	const (
		move = 4
		bias = 0.1
	)
	if power < 0 {
		return 0, 0
	}
	if tol < 0 {
		tol = c.tolNR
	}
	if tol < minTOLNR {
		tol = minTOLNR
	}
	vAir := v + c.wind
	for i := 1; i < 100; i++ {
		var (
			x     = math.Abs(c.cDrag * vAir)
			fSum  = c.fGR + vAir*x
			deriv = fSum + 2*v*x
			Δv    = (power - v*fSum) / deriv
		)
		if deriv <= 0 || Δv > move {
			v += move
			vAir += move
			continue
		}
		if math.Abs(Δv) < tol {
			c.iterNR += i
			return v + Δv*(1-Δv*bias), i
		}
		v += Δv
		vAir += Δv
	}
	if c.errmsg == nil {
		c.appendErr(": Newton-Raphson did not converge") //never happened for power > 0
	}
	return v, 0
}

func (c *BikeCalc) velError(power, vel float64) {
	if power <= c.minPower {
		return
	}
	vExact, iter := c.NewtonRaphson(power, minTOLNR, vel)
	c.iterNR -= iter
	if iter == 0 {
		return
	}
	c.callsErr++
	err := vel - vExact
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
func (c *BikeCalc) velGuess() float64 {
	v := c.prevVel * (1 + 12*(c.prevTan-c.tan))
	if v > 1 {
		return v
	}
	return 1
}

// VelFromPower returns speed for power using initial speed velguess.
// For power < 0 VelFromPower returns 0, false. For powers < minPower
// VelFreewheeling is returned. As velocity solver function is used the one
// set by SetVelSover function. NewtonRaphson by default.
func (c *BikeCalc) VelFromPower(power, velguess float64) (vel float64, ok bool) {
	var iter int
	if power < 0 {
		return 0, false
	}
	if power < c.minPower {
		vel, ok = c.VelFreewheeling(), true
		return
	}
	c.callsSolver++
	if forceERR && c.callsSolver == 100 {
		c.cDrag = 0
	}
	if vel = velguess; vel == 0 {
		vel = c.velGuess()
	}
	vel, iter = velSolver(c, power, -1, vel)

	if c.errmsg != nil {
		return vel, false
	}
	ok = true
	if iter > c.maxIter {
		c.maxIter = iter
	}
	if c.velErrors && ok {
		c.velError(power, vel)
	}
	c.prevVel = vel
	c.prevTan = c.tan
	return
}

// PowerFromVel returns power for speed v.
func (c *BikeCalc) PowerFromVel(v float64) (power float64) {
	c.callsPFV++
	return v * (c.fGR + c.cDrag*c.signSq(v))
}

// FlatPower returns power for speed v on flat (grade = 0) road.
func (c *BikeCalc) FlatPower(v float64) (power float64) {
	if v <= 0 {
		return 0
	}
	return v * (c.mgCrr + c.cDrag*v*v)
}

// FlatSpeed returns speed for power on flat (grade = 0) road.
func (c *BikeCalc) FlatSpeed(power float64) (vel float64) {
	if power <= 0 {
		return 0
	}
	c.storeState()
	c.SetGrade(0)
	c.wind = 0
	vel, _ = c.NewtonRaphson(power, minTOLNR, 8)
	c.restoreState()
	return
}

// GradeFromVelAndPower returns slope tangent for speed and power.
// The tangent tan is calculated from sin from
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

// VelFreewheeling return speed (>= 0) for power = 0.
// The speed is v from
// fGR + cDrag x (v + wind)^2 = 0
func (c *BikeCalc) VelFreewheeling() float64 {
	c.callsFW++

	v := math.Sqrt(math.Abs(c.fGR / c.cDrag))
	if c.fGR > 0 {
		v = -v
	}
	v -= c.wind
	if v < 0 {
		return 0
	}
	return v
}

// PowerFromVerticalUp returns power from vertical up speed (m/h).
// The power is calculated by riding uphill slope with tangent grade.
func (c *BikeCalc) PowerFromVerticalUp(v, grade float64) (power float64) {
	if v <= 0 || grade <= 0 {
		return math.NaN()
	}
	c.storeState()
	c.SetGrade(grade)

	v /= (c.sin * 3600) //speed on road m/s
	power = v * (c.fGR + c.cDrag*v*v)
	c.restoreState()
	return
}

// VerticalUpFromPower returns vertical speed up (m/h) for power.
// Speed is calculated by riding uphill slope with tangent grade.
func (c *BikeCalc) VerticalUpFromPower(power, grade float64) (velUp float64) {
	if power <= 0 || grade <= 0 {
		return math.NaN()
	}
	c.storeState()
	c.SetGrade(grade)
	c.wind = 0
	velUp, _ = c.NewtonRaphson(power, minTOLNR, 4)
	velUp *= c.sin * 3600 //speed up m/h
	c.restoreState()
	return
}

// CdAfromVelAndPower returns CdA from speed and power on flat ground road.
func (c *BikeCalc) CdAfromVelAndPower(v, power float64) (CdA float64) {
	if power <= 0 || v <= 0 {
		return math.NaN()
	}
	CdA = (power - v*c.mgCrr) / (0.5 * c.rho * v * v * v) // can be < 0
	return
}

// VelFromTurnRadius returns approx. max cornering speed for turn radius and
// cornering friction coefficient.
// The cornering speed is v from
// m x v^2 / radius = m x cos x gravity x ccf
func (c *BikeCalc) VelFromTurnRadius(radius, ccf float64) float64 {

	return math.Sqrt(radius * c.cos * c.gravity * ccf)
}

// LocalGravity returns local gravitational acceleration for latitude and elevation.
// https://en.wikipedia.org/wiki/Gravity_of_Earth
func (c *BikeCalc) LocalGravity(degLat, metersEle float64) (gravity float64) {
	const eleCorrection = (-3.086 + 1.1) * 1e-6

	cos := math.Cos(2 * degLat * (math.Pi / 180.0))
	gravity = 9.780327 * (1.0026454 - cos*(0.0026512-cos*0.0000058))
	gravity += eleCorrection * metersEle
	return
}

/*
 https://en.wikipedia.org/wiki/Density_of_air
 https://en.wikipedia.org/wiki/Barometric_formula#Density_equations

 Sea level rho at pressure P and temperature T
 		rho = P*M / (R*T)
 Elevation adjustment for sea level air density
 		x = ele * L / T
 		y = -(1.0 + g*M/(R*L))
 		rho = rho * (1+x)^y

 (1+x)^y = exp(z), z = y * ln(1+x)
 For small x and z we can use Taylor series:
 		ln(1+x) = x - x^2/2 + x^3/3 - x^4/4 ...
 		exp(z) = 1 + z + z^2/2 + z^3/6 + z^4/24 ...
 			==>
 		z = y * (x - x^2 /2)
 		rho = rho * (1 + z + z^2 /2 + z^3 /6)
 Error under 3000 m  < 0.005%
 Error at 5000 m 0.2%
*/

// RhoFromEle calculates air density and temperature for given elevation from
// base elevation and base elevation's temperature and air pressure.
// If air pressure is not known, calculation is based on
// the Standard Atmosphere sea level air pressure.
func (c *BikeCalc) RhoFromEle(ele float64) (rho, temperature float64) {

	const (
		M  = 0.0289644 //molar mass of Earth's air (kg/mol)
		R  = 8.3144598 //universal gas constant (N·m/(mol·K))
		P0 = 101325.0  //standard sea level atmospheric pressure (Pascals)
		g0 = 9.80665   //Standard sea level gravity. local gravity for lat 45.5 deg
		L  = -0.0065   //temperature lapse rate (per 1 m up)
	)
	P := c.airPressure
	if P <= 0 {
		P = P0
		c.airPressure = P0
		c.baseElevation = 0
	}
	g := c.gravity
	if g < 9.7 {
		g = g0
		c.gravity = g0
	}
	T := c.temperature + 273.15

	//rho at base elevation at temperature T and air pressure P
	rho = P * M / (R * T)

	// Calculate rho at dEle meters above or under base elevation
	dEle := ele - c.baseElevation
	x := dEle * L / T
	y := -(1.0 + g*M/(R*L))
	z := y * x * (1 - x*0.5)
	rho *= 1 + z*(1+z*(0.5+z*(1.0/6))) // rho *= math.Pow(1+x, y)
	temperature = c.temperature + dEle*L
	return
}
