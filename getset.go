package motion

import "math"

func (c *BikeCalc) SetBracket(ms float64) {
	if ms < minBracketLen {
		ms = minBracketLen
	}
	c.bracketLen = ms
}

func (c *BikeCalc) SetBaseElevation(m float64) { c.baseElevation = m }

func (c *BikeCalc) SetTemperature(C float64) { c.temperature = C }

func (c *BikeCalc) SetAirPressure(hPa float64) {
	if hPa > 0 {
		c.airPressure = hPa * 100 // hPa to Pascal
	}
}

func (c *BikeCalc) SetCdA(x float64) {
	if x > 0 {
		c.cdA = x
	}
	c.cDrag = 0.5 * c.cdA * c.rho
}

func (c *BikeCalc) SetCd(x float64) {
	if x > 0 {
		c.cd = x
	}
	c.SetCdA(c.frontalArea * c.cd)
}

func (c *BikeCalc) SetFrontalArea(x float64) {
	if x > 0 {
		c.frontalArea = x
	}
	c.SetCdA(c.frontalArea * c.cd)
}

// SetRho sets air density and updates the compound air drag coefficient cDrag.
func (c *BikeCalc) SetRho(x float64) {
	if x > 0 {
		c.rho = x
	}
	c.cDrag = 0.5 * c.cdA * c.rho
}

// SetCbf sets the braking friction coefficient.
func (c *BikeCalc) SetCbf(x float64) {
	if x > 0 {
		c.cbf = x
	}
	c.mgCbf = c.mg * c.cbf
	c.fBrake = c.cos * c.mgCbf
}

// SetCrr sets the rolling resistance coefficient and dependant forces fRoll and fGR.
func (c *BikeCalc) SetCrr(x float64) {
	if x > 0 {
		c.crr = x
	}
	c.mgCrr = c.mg * c.crr
	c.fRoll = c.cos * c.mgCrr
	c.fGR = c.fGrav + c.fRoll
}

func (c *BikeCalc) SetCcf(x float64) { c.ccf = x }

/*
 Calculate cos and sin as function of tan
 tan = grade% / 100
 Trigonometry: 	sin = sin(arctan(tan))
				cos = cos(arctan(tan))
 Pythagoras:	cos = 1/sqrt(1 + tan^2).
 				sin = tan * cos
 Taylor serie:
 		1/sqrt(1+s) = 1 - 1/2*s + 3/8*s^2 - 5/16s^3 + O(s^4)  =>
 		cos = 1 - 1/2*tan^2 + 3/8*tan^4 + O(tan^6)

 Taylor approximation with slightly modified coefficients for tan < 0.2
 		cos = 1 - tan^2*(0.4998 - tan^2*0.36)
		sin = tan * cos		// Error  ~< 1e-6 for tan < 0.2
 Hardware Sqrt is not slower than Taylor approximation.
*/

// SetGrade calculates sin and cos and dependant forces from tan = grade% /100.
func (c *BikeCalc) SetGrade(tan float64) {
	if c.tan == tan {
		return
	}
	c.tan = tan
	c.cos = 1 / math.Sqrt(1+tan*tan)
	c.sin = tan * c.cos
	c.setForces()
}

func (c *BikeCalc) setForces() {
	c.fBrake = c.cos * c.mgCbf
	c.fRoll = c.cos * c.mgCrr
	c.fGrav = c.sin * c.mg
	c.fGR = c.fGrav + c.fRoll
}

// SetGravity sets gravity and calculates all weight dependant forces and intermediates.
func (c *BikeCalc) SetGravity(g float64) {
	if g > 0 {
		c.gravity = g
	}
	c.SetWeight(0)
}

// SetTolNR sets iteration stop tolerance for Newton-Raphson method.
func (c *BikeCalc) SetTolNR(ms float64) {
	if ms < minTolNR {
		ms = minTolNR
	}
	c.tolNR = ms
}

func (c *BikeCalc) SetMinPower(w float64) {
	if w > 0 {
		c.minPower = w
	}
}

// func (c *BikeCalc) SetPower(x float64) { c.power = x }

func (c *BikeCalc) SetVelErrors(b bool) { c.velErrors = b }

func (c *BikeCalc) SetVelSolver(i int) {
	c.solverFunc = i
	switch i {
	case newtonRaphson:
		velFromPower = (*BikeCalc).NewtonRaphson

	case singleQuadratic:
		velFromPower = (*BikeCalc).Quadratic

	case doubleQuadratic:
		velFromPower = (*BikeCalc).DoubleQuadratic

	case doubleLinear:
		velFromPower = (*BikeCalc).DoubleLinear

	default:
		velFromPower = (*BikeCalc).NewtonRaphson
		c.solverFunc = newtonRaphson
	}
	c.tolVel = c.tolNR
	if c.solverFunc != newtonRaphson {
		c.tolVel = c.bracketLen
	}
}

// SetWeight sets total weight and dependant forces.
func (c *BikeCalc) SetWeight(kg float64) {
	if kg > 0 {
		c.mass = kg
	}
	c.mg = c.mass * c.gravity
	c.mgCrr = c.mg * c.crr
	c.mgCbf = c.mg * c.cbf
	c.setForces()
	c.SetWeightRotating(0)
}

// SetWeightRotating sets rotating mass.
// rotatingMass = weight (tyre + tube + rim)
func (c *BikeCalc) SetWeightRotating(kg float64) {
	if kg > 0 {
		c.massRotating = kg
	}
	// aproximating 0.9 reduction of the wheel radius, because
	// the mass is not rotating at the outer edge of the wheel.
	c.massKin = c.mass + c.massRotating*0.9
	c.oMassKin = 1 / c.massKin
}

// SetWind sets +head/-tail wind speed in m/s.
func (c *BikeCalc) SetWind(ms float64) { c.wind = ms }

func (c *BikeCalc) AirPressure() float64 { return c.airPressure / 100 }

func (c *BikeCalc) BaseElevation() float64 { return c.baseElevation }

func (c *BikeCalc) Cbf() float64 { return c.cbf }

func (c *BikeCalc) Crr() float64 { return c.crr }

func (c *BikeCalc) CdA() float64 { return c.cdA }

func (c *BikeCalc) Cdrag() float64 { return c.cDrag }

func (c *BikeCalc) MassKin() float64 { return c.massKin }

func (c *BikeCalc) Fbrake() float64 { return c.fBrake }

func (c *BikeCalc) Fdrag(v float64) float64 { return c.cDrag * c.signSq(v) }

func (c *BikeCalc) Fgrav() float64 { return c.fGrav }

func (c *BikeCalc) Froll() float64 { return c.fRoll }

func (c *BikeCalc) Fgr() float64 { return c.fGR }

func (c *BikeCalc) Gravity() float64 { return c.gravity }

func (c *BikeCalc) Grade() float64 { return c.tan }

func (c *BikeCalc) Rho() float64 { return c.rho }

func (c *BikeCalc) SolverCalls() int { return c.callsSolver }

func (c *BikeCalc) MaxIter() int { return c.maxIter }

func (c *BikeCalc) SolverRounds() int {
	if c.solverFunc == newtonRaphson {
		return c.iterNR
	}
	return c.callsf
}

func (c *BikeCalc) SolverFunc() int { return c.solverFunc }

// func (c *BikeCalc) FreewheelCalls() int { return c.callsfw }

// func (c *BikeCalc) PowerFromVelCalls() int { return c.callspfv }

func (c *BikeCalc) VelErrorMax() float64 { return c.velErrMax }

func (c *BikeCalc) VelErrorMean() float64 {
	if c.callsErr == 0 {
		return 0
	}
	return c.velErr / float64(c.callsErr)
}

func (c *BikeCalc) VelErrorAbsMean() float64 {
	if c.callsErr == 0 {
		return 0
	}
	return c.velErrAbs / float64(c.callsErr)
}

func (c *BikeCalc) VelErrorPos() float64 {
	if c.callsErr == 0 {
		return 0
	}
	return float64(c.velErrPos) / float64(c.callsErr)
}
