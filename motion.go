package motion

import "math"

const (
	forceError		= false
	minBracketLen  	= 1e-10
	minTolNR       	= 1e-12
	maxSpeed       	= 1000.0
	maxSpeedS		= "1000 m/s"
	stdAirPressure 	= 101325.0 //standard sea level atmospheric pressure (Pascals)
)
const (
	newtonRaphson       = 1
	newtonHalley		= 2
	singleQuadratic		= 3
	doubleQuadratic     = 4
	singleLinear		= 5
	doubleLinear        = 6
	threeLinear         = 7
	bisect              = 8
)



type store struct {
	grade float64
	power float64
	wind  float64
}

// BikeCalc --
type BikeCalc struct {
	errmsg			[]byte
	store         	store
	
	solverNo		int
	gravity       	float64 
	temperature   	float64
	airPressure   	float64
	rho           	float64
	baseElevation 	float64

	cd           	float64
	frontalArea  	float64
	cdA          	float64
	crr          	float64
	cbf          	float64
	ccf          	float64
	weight       	float64
	weightWheels 	float64
	
	wind  			float64
	power 			float64

	cDrag			float64
	sin      		float64
	cos      		float64
	tan   			float64
	fGrav    		float64
	fRoll    		float64
	fGR      		float64
	fBrake   		float64
	mg       		float64
	mgCrr    		float64
	mgCbf    		float64
	cMassKin 		float64
	massKin  		float64

	tolNR      		float64
	minPower   		float64
	bracketLen 		float64

	prevVel    		float64
	prevTan    		float64

	roundsNR    	int
	maxIter			int
	callsFunc    	int
	callsSolver  	int
	callsFW      	int
	callsPFV     	int
	velErr       	float64
	velErrAbs    	float64
	velErrMax    	float64
	velErrSS     	float64
	velErrPos    	int
	velErrOverCL 	int
	callsErr     	int
	velErrors    	bool

}

var solverFunc func (*BikeCalc, float64, float64, float64) (float64, int)

// Calculator returns new semi-initialized calculator
func Calculator() *BikeCalc {
	c := &BikeCalc{
		gravity:       	9.80665, //Standard sea level value (for lat 45.5 deg)
		rho:           	1.2250,  //Standard sea level value at 15 deg C
		temperature:   	15,      //Standars deg Celsius.
		airPressure:   	stdAirPressure, //standard sea level
		baseElevation: 	0,
		cbf:           	0.25,
		ccf:           	0.175,
		cdA:         	0.7,
		crr:          	0.007,
		cos:			1,
		bracketLen:    	1.25,
		prevVel:       	5,
		tolNR:         	0.05,
		velErrors:     	false,
		minPower:      	1,
	}
	c.SetVelSolver(newtonRaphson)
	return c
}

type errstr struct {
	s string
}
func (e *errstr) Error() string {
	return e.s
}
// Error returns cumulative error and resets it to nil.
func (c *BikeCalc) Error() error {
	if len(c.errmsg) == 0 {
		return nil
	}
	e := &errstr{string(c.errmsg)}
	c.errmsg = nil
	return  e
}
func (c *BikeCalc) appendErr(s string) {
	c.errmsg = append(c.errmsg, s...)
}

// this makes functions using it more inlineable.
func abs(x float64) float64 {
	if x > 0 {
		return x
	}
	return -x
}

func max(x, y float64) float64 {
	if x > y {
		return x
	}
	return y
}

func (c *BikeCalc) storeState() {
	c.store.grade = c.tan
	c.store.power = c.power
	c.store.wind = c.wind
}

func (c *BikeCalc) restoreState() {
	c.SetGradeExact(c.store.grade)
	c.power = c.store.power
	c.wind = c.store.wind
}

// NewtonRaphson condition "deriv <= 0 || Δv > 3" handles +/-Inf, NaN and too long right step.
// deriv > 0 --> too long left step is not possible or very rare and is recoverable.

// NewtonRaphson returns speed for power and number of iterations used. 
// Function is globally convergent to the rightmost root.
// iter = 0 --> calculation failed, which until now has not been seen for power >= 1 w, 
// air drag coeffient cDrag clearly > 0, reasonable initial speed and tol > 0.01.
// Iteration is stopped when change in speed is less than tol. This gives error <= tol^2 x 0.1.
// v is initial speed for iteration.
//
func (c *BikeCalc) NewtonRaphson(power, tol, v float64) (vel float64, iter int) {
	const bias = 0.1 //bias correction coefficient
	const iterLim = 100
	if tol < minTolNR {
		tol = minTolNR
	}
	// fGR, cDrag := c.fGR, c.cDrag
	vAir := v + c.wind
	for i := 1; i < iterLim; i++ {
		c.roundsNR++

		s := 		abs(c.cDrag * vAir)
		forces := 	c.fGR + s*vAir
		deriv := 	forces + 2*s*v
		Δv := 		(power - v*forces) / deriv

		if deriv <= 0 || Δv > 3 {
			v += 3; vAir += 3
			continue
		} 
		if abs(Δv) < tol {
			return v + Δv * (1 - Δv*bias), i
		}
		vAir += Δv; v += Δv; 
	}
	if c.errmsg == nil {
		c.appendErr(": Newton-Raphson did not converge")
	} 
	return v, 0
}

// NewtonRaphsonHalley is NewtonRaphson with Halley acceleration for the final round.
// If Newton condition "abs(Δv) < tol" fails, Halley iteration can give the same accuracy for
// condition "abs(Δv) < 4.0*tol". Halley acceleration can save 0.3 - 0.5 iterations but 
// this doesn't not make it measurable faster than plain NewtonRaphson.
// 2. derivative = 2*c.cDrag*sign(vAir)*(v + 2*vAir) is discontinuous at vAir = 0.
// So we use Halley only when vAir is clearly > 0. Also at root, approximately.
//
func (c *BikeCalc) NewtonRaphsonHalley(power, tol, v float64) (vel float64, iter int) {
	const  bias = 0.1
	const iterLim = 100
	if tol < minTolNR {
		tol = minTolNR
	}
	tolHal := 4.0*tol
	vAir := v + c.wind

	for i := 1; i < iterLim; i++ {
		c.roundsNR++

		s := 		abs(c.cDrag * vAir)
		forces := 	c.fGR + s*vAir
		deriv := 	forces + 2*s*v
		Δv := 		(power - v*forces) / deriv

		if deriv <= 0 || Δv > 3 {
			v += 3; vAir += 3
			continue
		} 
		if abs(Δv) < tol {
			return v + Δv * (1 - Δv*bias), i
		}
		if vAir > tolHal && abs(Δv) < tolHal {
			return v + Δv * deriv/(deriv + Δv*c.cDrag*(v + 2*vAir)), i
		}
		v += Δv; vAir += Δv
	}
	if c.errmsg == nil {
		c.appendErr("Newton-Raphson-Halley did not converge")
	}  
	return v, 0
}

//VelError --
func (c *BikeCalc) VelError(power, vel float64) {
	
	if power <= c.minPower {
		return 
	}
	vExact, iter := c.NewtonRaphson(power, 1e-13, vel)
	c.roundsNR -= iter
	if iter == 0 {
		return 
	}
	err := vel - vExact
	if err > -2.2e-16*vel {
		c.velErrPos++
	}
	c.velErr += err
	c.velErrAbs += abs(err)
	if abs(err) > c.velErrMax {
		c.velErrMax = abs(err)
	}
	c.velErrSS += err * err
	c.callsErr++
	if err*err > 3.8416*c.velErrSS/float64(c.callsErr) {
		c.velErrOverCL++
	}
	return 
}

// VelGuess returns speed estimate by adjusting previous speed by road grade change.
func (c *BikeCalc) VelGuess() float64 {
	v := c.prevVel * (1 + 14*(c.prevTan-c.tan))
	if v > 1 {
		return v
	}
	return 1
}

//SetPrevVel ...
func (c *BikeCalc) SetPrevVel(v float64) {
	c.prevVel = v
	c.prevTan = c.tan
}

// VelFromPower returns speed for power using initial speed velguess. 
// Default or set by SetVelSolver velocity solver is used.
func (c *BikeCalc) VelFromPower(power, velguess float64) (vel float64, ok bool) {

	if power < 0 {
		return 0, false
	}
	if power < c.minPower {
		vel, ok = c.VelFreewheeling(), true
		return 
	}
	c.callsSolver++
	if forceError && c.callsSolver >= 100 {
		c.cDrag = 0
	}
	if vel = velguess; vel == 0 {
		vel = c.VelGuess()
	}
	tol := c.tolNR
	if c.solverNo > newtonHalley && c.solverNo < bisect {
		tol = c.bracketLen
	}
	vel, iter := solverFunc(c, power, tol, vel)

	if ok = iter > 0; !ok {
		return c.tryBisect(power)
	}
	if iter > c.maxIter {
		c.maxIter = iter
	}
	c.SetPrevVel(vel)
	if c.velErrors {
		c.VelError(power, vel)
	}
	return 
}

func (c *BikeCalc) tryBisect(power float64) (float64, bool) {
	if c.solverNo == bisect {
		return 0, false
	}
	if len(c.errmsg) < 35 { 
		c.appendErr(": Bisect used")
	} 
	v, i := c.Bisect(power, 0, 32) // tol = 0 --> NR accuracy
	return v, i > 0
}

// PowerFromVel returns power for speed v.
func (c *BikeCalc) PowerFromVel(v float64) (power float64) {
	c.callsPFV++
	return v*(c.fGR + c.cDrag*abs(v+c.wind)*(v+c.wind))
}

// FlatPower returns power for speed v on flat (slope=0) road and no wind.
func (c *BikeCalc) FlatPower(v float64) (power float64) {
	if v <= 0 {
		return 0
	}
	return v*(c.mgCrr + c.cDrag*v*v)
}

// FlatSpeed returns speed for power on flat (slope=0) road and no wind.
func (c *BikeCalc) FlatSpeed(power float64) (vel float64, ok bool) {

	if power < c.minPower {
		return 0, false
	}
	c.storeState()
	defer c.restoreState()
	c.wind = 0
	c.SetGrade(0)
	vel, iter := c.NewtonRaphson(power, minTolNR, 6)
	ok = iter > 0
	return 
}

// GradeFromVelAndPower returns slope tangent for speed, power and wind = 0.
// Equation: sin*mg + sqrt(1-sin^2)*mg*Crr + cDrag*v^2 - power/v = 0
func (c *BikeCalc) GradeFromVelAndPower(v, power float64) (tan float64) {
	if v <= 0 {
		return math.NaN()
	}
	sin := (power/v - c.cDrag*v*v) / c.mg 	//sin for crr = 0
	if abs(sin) > 1 {
		return math.NaN()
	}
	r := c.crr * c.crr
	d := r * (1 + r - sin*sin)
	sin = (sin - math.Sqrt(d)) / (1 + r)
	tan = sin / math.Sqrt(1 - sin*sin)
	return 
}

// VelFreewheeling return speed (>= 0) for power = 0.
// Equation: fGR + cDrag x (v + wind)^2 = 0
func (c *BikeCalc) VelFreewheeling() (v float64) {
	c.callsFW++ 

	v = math.Sqrt(abs(c.fGR / c.cDrag))
	if c.fGR > 0 { v = -v }
	v -= c.wind
	if v < 0 { v = 0 }
	return 
}

// PowerFromVerticalUp returns power from vertical up speed (m/h). 
// Power is calculated by riding uphill slope with tangent grade.
func (c *BikeCalc) PowerFromVerticalUp(v, grade float64) (power float64) {
	if v < 0 || grade < 0.01 {
		return 0
	}
	c.storeState()
	defer c.restoreState()
	c.SetGrade(grade)

	v /= (c.sin * 3600) //speed on road m/s
	return v * (c.fGR + c.cDrag*v*v)
}

// VerticalUpFromPower returns vertical speed up (m/h) for power. 
// Speed is calculated by riding uphill slope with tangent grade.
func (c *BikeCalc) VerticalUpFromPower(power, grade float64) (float64, bool) {
	if power <= 0 || grade < 0.01 {
		return 0, false
	}
	c.storeState()
	defer c.restoreState()
	c.SetGrade(grade)
	c.wind = 0
	v, iter := c.NewtonRaphson(power, minTolNR, 4)
	return  v * c.sin * 3600, iter > 0
}

// CdAfromVelAndPower returns CdA from speed and power on flat road.
func (c *BikeCalc) CdAfromVelAndPower(vel, power float64) (CdA float64) {
	if power <= 0 || vel <= 0 {
		return 0
	}
	CdA = (power - vel*c.mgCrr) / (0.5 * c.rho * vel * vel * vel)
	if CdA < 0 {
		CdA = 0
	}
	return 
}

// VelFromTurnRadius returns turn speed for turn radius and preset 
// cornering friction coefficient.
func (c *BikeCalc) VelFromTurnRadius(radius float64) (vel float64) {
	const minCcf = 0.05
	const slopeFactor = 0.5
	ccf := c.ccf
	if c.tan < 0 {
		ccf += c.tan * slopeFactor
		if ccf < minCcf {
			ccf = minCcf
		}
	}
	return math.Sqrt(radius * c.gravity * ccf)
}

// LocalGravity calculates local gravitational acceleration for latitude and elevation.
// Elevation correction for mean density rock underneath.
// https://en.wikipedia.org/wiki/Gravity_of_Earth
func (c *BikeCalc) LocalGravity(degLat, metersEle float64) (gravity float64) {
	const eleCorrection = (-3.086 + 1.1) * 1e-6
	
	cos := math.Cos(2 * degLat * math.Pi / 180)
	gravity = 9.780327 * (1.0026454 - cos*(0.0026512 - cos*0.0000058))
	gravity += eleCorrection * metersEle
	return
}

// https://en.wikipedia.org/wiki/Density_of_air
// https://en.wikipedia.org/wiki/Barometric_formula#Density_equations
//
// Sea level rho at pressure P and temperature T
// 		rho = P*M / (R*T)
// Elevation adjustment for sea level air density
// 		x = ele * L / T
// 		y = -(1.0 + g*M/(R*L))
// 		rho = rho * (1+x)^y
//
// (1+x)^y = exp(z), z = y * ln(1+x)
// For small x and y we can use Taylor series:
// 		ln(1+x) = x - x^2/2 + x^3/3 - x^4/4 ...
// 		exp(z) = 1 + z + z^2/2 + z^3/6 + z^4/24 ...
// 			==>
// 		z = y * (x - x^2 /2)
// 		rho = rho * (1 + z + z^2 /2 + z^3 /6)
// Error under 3000 m  < 0.005%
// Error at 5000 m 0.2%

// RhoFromEle calculates air density and temperature for given elevation from
// preset base elevation and base elevation's temperature and air density.
func (c *BikeCalc) RhoFromEle(ele float64) (rho, temperature float64) {

	const (
		M  = 0.0289644 //molar mass of Earth's air (kg/mol)
		R  = 8.3144598 //universal gas constant (N·m/(mol·K))
		P0 = 101325.0  //standard sea level atmospheric pressure (Pascals)
		g0  = 9.80665  //Standard sea level gravity. local gravity for lat 45.5 deg
		L  = -0.0065   //temperature lapse rate (per 1 m up)
	)
	baseEle := c.baseElevation
	P := c.airPressure
	if P <= 0 {
		P = P0
		baseEle = 0
	}
	g := c.gravity
	if g < 9.7 {
		g = g0
	}
	T := c.temperature + 273.15

	//rho at base elevation at temperature T and air pressure P
	rho = P * M / (R * T)

	//calculate rho at dEle meters above or under baseElevation
	dEle := ele - baseEle
	x := dEle * L / T
	y  := -(1.0 + g*M/(R*L))
	// "exact" rho
	//rho *= math.Pow(1+x, y)
	z := y * x*(1 - x*0.5)
	rho *= 1 + z*(1 + z*(0.5 + z* (1.0 / 6.0)))
	temperature = c.temperature + dEle*L
	return 
}
