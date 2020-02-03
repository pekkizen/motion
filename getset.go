package motion

import "math"

// SetBracket ...
func (c *BikeCalc) SetBracket(x float64) {
	if x < minBracketLen {
		x = minBracketLen
	}
	c.bracketLen = x
}

// SetBaseElevation ...
func (c *BikeCalc) SetBaseElevation(x float64) {
	if x >= 0 {
		c.baseElevation = x
	}
}

// SetTemperature ...
func (c *BikeCalc) SetTemperature(x float64) {
	c.temperature = x
}

// SetAirPressure ...
func (c *BikeCalc) SetAirPressure(x float64) {
	if x > 0 {
		c.airPressure = x * 100 //hp to pascals
	}
}

// SetCdA ...
func (c *BikeCalc) SetCdA(x float64) {
	if x > 0 {
		c.cdA = x
	}
	c.cDrag = 0.5 * c.cdA * c.rho
}

// SetCd ...
func (c *BikeCalc) SetCd(x float64) {
	if x > 0 {
		c.cd = x
	}
	c.SetCdA(c.frontalArea * c.cd)
}

// SetFrontalArea ...
func (c *BikeCalc) SetFrontalArea(x float64) {
	if x > 0 {
		c.frontalArea = x
	}
	c.SetCdA(c.frontalArea * c.cd)
}

// SetRho sets air density and calculates compound air drag coefficient cDrag.
func (c *BikeCalc) SetRho(x float64) {
	if x > 0 {
		c.rho = x
	}
	c.cDrag = 0.5 * c.cdA * c.rho
}

// SetCbf sets coefficient of braking friction.
func (c *BikeCalc) SetCbf(x float64) {
	if x > 0 {
		c.cbf = x
	}
	c.mgCbf = c.mg * c.cbf
	c.fBrake = c.cos * c.mgCbf
}

// SetCcf sets coefficient of cornering friction.
func (c *BikeCalc) SetCcf(x float64) {
	if x > 0 {
		c.ccf = x
	}
}

// SetCrr sets coefficient of rolling resistance and dependant forces fRoll and fGR.
func (c *BikeCalc) SetCrr(x float64) {
	if x > 0 {
		c.crr = x
	}
	c.mgCrr = c.mg * c.crr
	c.fRoll = c.cos * c.mgCrr
	c.fGR = c.fGrav + c.fRoll
}

// Calculate cos and sin as function of tan
// tan = grade% / 100
// Trigonometry: sin = sin(arctan(tan))
//				 cos = cos(arctan(tan))
// Pythagoras:	 cos = 1/sqrt(1 + tan^2).
// 				 sin = tan * cos
// Taylor serie:
// 		1/sqrt(1+s) = 1 - 1/2*s + 3/8*s^2 - 5/16s^3 + O(s^4)  =>
// 		cos = 1 - 1/2*tan^2 + 3/8*tan^4 + O(tan^6)
// 
// Taylor approximation with slightly modified coefficients for tan < 0.2
// 		cos = 1 - tan^2*(0.4998 - tan^2*0.36)		
//		sin = tan * cos								// Error in sin < 1e-6 for grade < 20%

// SetGrade calculates sin and cos and dependant forces from tan = grade% / 100
func (c *BikeCalc) SetGrade(tan float64) {
	c.tan = tan
	// Without hardware sqrt use these two lines for c.cos.
	// tan  *= tan
	// c.cos = 1 - tan*(0.4998 - tan*0.36)
	c.cos = 1 / math.Sqrt(1 + tan*tan)
	c.sin = c.tan * c.cos
	c.setForces()
}

// SetGradeExact uses exact formula, always.
func (c *BikeCalc) SetGradeExact(tan float64) {
	c.tan = tan
	c.cos = 1 / math.Sqrt(1 + tan*tan)
	c.sin = tan * c.cos
	c.setForces()
}

func (c *BikeCalc) setForces() {
	c.fBrake = c.cos * c.mgCbf
	c.fRoll  = c.cos * c.mgCrr
	c.fGrav  = c.sin * c.mg
	c.fGR    = c.fGrav + c.fRoll
}

// SetGravity sets gravity and calculates all weight dependant forces and intermediates.
func (c *BikeCalc) SetGravity(x float64) {
	if x > 0 {
		c.gravity = x
	}
	c.SetWeight(0)
}

// SetMinPower sets power limit between freewheeling speeds and powered speeds.
func (c *BikeCalc) SetMinPower(x float64) {
	c.minPower = x
}

// SetVelTol sets iteration stop tolerance for Newton-Raphson and Bisect methods.
func (c *BikeCalc) SetVelTol(x float64) {
	if x < minTolNR {
		x = minTolNR
	}
	c.tolNR = x
}

// SetPower ...
func (c *BikeCalc) SetPower(x float64) {
	c.power = x
}

// SetVelErrors ...
func (c *BikeCalc) SetVelErrors(b bool) {
	c.velErrors = b
}

// SetVelSolver sets function for solving speed from power.
// 
func (c *BikeCalc) SetVelSolver(i int) {
	if i < newtonRaphson || i > bisect {
		i = newtonRaphson
	}
	c.solverNo = i
	switch c.solverNo {
	case newtonRaphson:
		solverFunc = (*BikeCalc).NewtonRaphson
	case newtonHalley:
		solverFunc = (*BikeCalc).NewtonRaphsonHalley
	case singleQuadratic:
		// solverFunc = (*BikeCalc).SingleQuadratic
		solverFunc =  (*BikeCalc).SingleQuadraticTriple
	case doubleQuadratic:
		solverFunc = (*BikeCalc).DoubleQuadratic
	case singleLinear:
		solverFunc = (*BikeCalc).SingleLinear
	case doubleLinear:
		solverFunc = (*BikeCalc).DoubleLinear
	case threeLinear:
		solverFunc = (*BikeCalc).ThreeLinear
	case bisect:
		solverFunc = (*BikeCalc).Bisect
	}
}

// SetWeight ...
func (c *BikeCalc) SetWeight(x float64) {
	if x > 0 {
		c.weight = x
	}
	c.mg = c.weight * c.gravity
	c.mgCrr = c.mg * c.crr
	c.mgCbf = c.mg * c.cbf
	c.SetWeightWheels(0)
	c.setForces()
}

// SetWeightWheels ...
func (c *BikeCalc) SetWeightWheels(x float64) {
	if x > 0 {
		c.weightWheels = x
	}
	c.massKin = c.weight + c.weightWheels
	c.cMassKin = 0.5 * c.massKin
}

// SetWind ...
func (c *BikeCalc) SetWind(x float64) {
	c.wind = x
}

// AirPressure ...
func (c *BikeCalc) AirPressure() float64 {
	return c.airPressure
}

// BaseElevation ...
func (c *BikeCalc) BaseElevation() float64 {
	return c.baseElevation
}

// Cbf ...
func (c *BikeCalc) Cbf() float64 {
	return c.cbf
}

// CdA --
func (c *BikeCalc) CdA() float64 {
	return c.cdA
}

// Cdrag ...
func (c *BikeCalc) Cdrag() float64 {
	return c.cDrag
}

// CmassKin ...
func (c *BikeCalc) CmassKin() float64 {
	return c.cMassKin
}

// Fbrake ...
func (c *BikeCalc) Fbrake() float64 {
	return c.fBrake
}

// Fdrag ...
func (c *BikeCalc) Fdrag(v float64) float64 {
	return c.cDrag * abs(v + c.wind) * (v + c.wind)
}

// Fgrav ...
func (c *BikeCalc) Fgrav() float64 {
	return c.fGrav
}

// Froll ...
func (c *BikeCalc) Froll() float64 {
	return c.fRoll
}

// Fgr ...
func (c *BikeCalc) Fgr() float64 {
	return c.fGR
}

// Gravity ...
func (c *BikeCalc) Gravity() float64 {
	return c.gravity
}
// Grade ...
func (c *BikeCalc) Grade() float64 {
	return c.tan
}

// MaxIter ...
func (c *BikeCalc) MaxIter() int {
	return c.maxIter
}

// Rho ...
func (c *BikeCalc) Rho() float64 {
	return c.rho
}

// SolverCalls ...
func (c *BikeCalc) SolverCalls() int {
	return c.callsSolver
}

// SolverRounds ...
func (c *BikeCalc) SolverRounds() int {
	if c.solverNo <= newtonHalley {
		return c.roundsNR
	}
	return c.callsFunc
}

// FreewheelCalls ...
func (c *BikeCalc) FreewheelCalls() int {
	return c.callsFW
}
// PowerFromVelCalls ...
func (c *BikeCalc) PowerFromVelCalls() int {
	return c.callsPFV
}

// VelErrorSD ...
func (c *BikeCalc) VelErrorSD() float64 {
	if c.callsErr == 0 {
		return 0
	}
	return math.Sqrt(c.velErrSS / float64(c.callsErr))
}

// VelErrorOverCL ...
func (c *BikeCalc) VelErrorOverCL() float64 {
	if c.callsErr == 0 {
		return 0
	}
	return float64(c.velErrOverCL) / float64(c.callsErr) 
}

// VelErrorMax ...
func (c *BikeCalc) VelErrorMax() float64 {
	return c.velErrMax
}

// VelErrorMean ...
func (c *BikeCalc) VelErrorMean() float64 {
	if c.callsErr == 0 {
		return 0
	}
	return c.velErr / float64(c.callsErr)
}

// VelErrorAbsMean ...
func (c *BikeCalc) VelErrorAbsMean() float64 {
	if c.callsErr == 0 {
		return 0
	}
	return c.velErrAbs / float64(c.callsErr)
}

// VelErrorPos ...
func (c *BikeCalc) VelErrorPos() float64 {
	if c.callsErr == 0 {
		return 0
	}
	return float64(c.velErrPos) / float64(c.callsErr)
}