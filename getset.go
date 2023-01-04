package motion

import "math"

func (c *BikeCalc) SetBracket(x float64) {
	if x < minBracketLEN {
		x = minBracketLEN
	}
	c.bracketLen = x
}

func (c *BikeCalc) SetBaseElevation(x float64) {
	if x >= 0 {
		c.baseElevation = x
	}
}

func (c *BikeCalc) SetTemperature(x float64) {
	c.temperature = x
}

func (c *BikeCalc) SetAirPressure(x float64) {
	if x < 0 || x > 900 {
		c.airPressure = x * 100 // hPa to pascals
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

// SetCbf sets coefficient of braking friction.
func (c *BikeCalc) SetCbf(x float64) {
	if x > 0 {
		c.cbf = x
	}
	c.mgCbf = c.mg * c.cbf
	c.fBrake = c.cos * c.mgCbf
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
//		sin = tan * cos							// Error in sin < 1e-6 for grade < 20%

// SetGrade calculates sin and cos and dependant forces from tan = grade % / 100
func (c *BikeCalc) SetGrade(tan float64) {
	c.tan = tan
	// c.cos = 1 - tan*tan*(0.4998 - tan*tan*0.36)
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
func (c *BikeCalc) SetGravity(x float64) {
	if x > 0 {
		c.gravity = x
	}
	c.SetWeight(0)
}

// SetTolNR sets iteration stop tolerance for Newton-Raphson methods.
func (c *BikeCalc) SetTolNR(x float64) {
	if x < minTOLNR {
		x = minTOLNR
	}
	c.tolNR = x
}

func (c *BikeCalc) SetMinPower(x float64) {
	if x > 0 {
		c.minPower = x
	}

}

func (c *BikeCalc) SetPower(x float64) {
	c.power = x
}

func (c *BikeCalc) SetVelErrors(b bool) {
	c.velErrors = b
}

// SetVelSolver selects function for solving speed from power.
func (c *BikeCalc) SetVelSolver(i int) {
	if i < newtonRaphson || i > singleQuadratic {
		i = newtonRaphson
	}
	c.velSolver = i
	// switch c.solverNo {
	// case newtonRaphson:
	// 	solverFunc = (*BikeCalc).NewtonRaphson
	// case singleQuadratic:
	// 	solverFunc = (*BikeCalc).Quadratic
	// case singleLinear:
	// 	solverFunc = (*BikeCalc).SingleLinear
	// case doubleLinear:
	// 	solverFunc = (*BikeCalc).DoubleLinear
	// case doubleQuadratic:
	// 	solverFunc = (*BikeCalc).DoubleQuadratic
	// case bisect:
	// 	solverFunc = (*BikeCalc).Bisect
	// }
}

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

func (c *BikeCalc) SetWeightWheels(x float64) {
	if x > 0 {
		c.weightWheels = x
	}
	c.massKin = c.weight + c.weightWheels*0.9
	c.oMassKin = 1 / c.massKin
}

func (c *BikeCalc) SetWind(x float64) {
	c.wind = x
}

func (c *BikeCalc) AirPressure() float64 {
	return c.airPressure / 100 // hPa
}

func (c *BikeCalc) BaseElevation() float64 {
	return c.baseElevation
}

func (c *BikeCalc) Cbf() float64 {
	return c.cbf
}

func (c *BikeCalc) CdA() float64 {
	return c.cdA
}

func (c *BikeCalc) Cdrag() float64 {
	return c.cDrag
}

func (c *BikeCalc) MassKin() float64 {
	return c.massKin
}

func (c *BikeCalc) Fbrake() float64 {
	return c.fBrake
}

func (c *BikeCalc) Fdrag(v float64) float64 {
	v += c.wind
	return c.cDrag * math.Abs(v) * v
}

func (c *BikeCalc) Fgrav() float64 {
	return c.fGrav
}

func (c *BikeCalc) Froll() float64 {
	return c.fRoll
}

func (c *BikeCalc) Fgr() float64 {
	return c.fGR
}

func (c *BikeCalc) Gravity() float64 {
	return c.gravity
}

func (c *BikeCalc) Grade() float64 {
	return c.tan
}

func (c *BikeCalc) Rho() float64 {
	return c.rho
}

func (c *BikeCalc) SolverCalls() int {
	return c.callsSolver
}

func (c *BikeCalc) MaxIter() int {
	return c.maxIter
}

func (c *BikeCalc) SolverRounds() int {
	if c.velSolver <= newtonRaphson {
		return c.iterNR
	}
	return c.callsFun
}

func (c *BikeCalc) FreewheelCalls() int {
	return c.callsFW
}

func (c *BikeCalc) PowerFromVelCalls() int {
	return c.callsPFV
}

func (c *BikeCalc) VelErrorMax() float64 {
	return c.velErrMax
}

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
