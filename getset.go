package motion

import "math"

func (c *BikeCalc) SetBracket(x float64) {
	if x < minBracketLen {
		x = minBracketLen
	}
	c.bracketLen = x
	if c.solverFunc != newtonRaphson {
		c.tolVel = c.bracketLen
	}
}

// SetTolNR sets iteration stop tolerance for Newton-Raphson method.
func (c *BikeCalc) SetTolNR(x float64) {
	if x < minTolNR {
		x = minTolNR
	}
	c.tolNR = x
	if c.solverFunc == newtonRaphson {
		c.tolVel = c.tolNR
	}
}

func (c *BikeCalc) SetBaseElevation(m float64) { c.baseElevation = m }

func (c *BikeCalc) SetTemperature(t float64) { c.temperature = t }

func (c *BikeCalc) SetAirPressure(hPa float64) {
	c.airPressure = hPa * 100 // hPa to Pascal
}

func (c *BikeCalc) SetCdA(x float64) {
	if x > 0 {
		c.cdA = x
	}
	if c.rho > 0 && c.cdA > 0 {
		c.cDrag = 0.5 * c.cdA * c.rho
	}
}

func (c *BikeCalc) SetCd(x float64) {
	if x > 0 {
		c.cd = x
	}
	if c.frontalArea > 0 && c.cd > 0 {
		c.SetCdA(c.frontalArea * c.cd)
	}
}

func (c *BikeCalc) SetFrontalArea(x float64) {
	if x > 0 {
		c.frontalArea = x
	}
	if c.cd > 0 && c.frontalArea > 0 {
		c.SetCdA(c.frontalArea * c.cd)
	}
}

// SetRho sets air density and updates the compound air drag coefficient cDrag.
func (c *BikeCalc) SetRho(x float64) {
	if x > 0 {
		c.rho = x
	}
	if c.cdA > 0 && c.rho > 0 {
		c.cDrag = 0.5 * c.cdA * c.rho
	}
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
	 ////////////////////////////////////////
	 tan = grade% / 100
	 Trigonometry: 	sin = sin(arctan(tan))
					cos = cos(arctan(tan))
	 Pythagoras:	cos = 1/sqrt(1 + tan^2).
	 				sin = tan * cos
	 Taylor serie:
	 		1/sqrt(1+s) = 1 - 1/2*s + 3/8*s^2 - 5/16*s^3 + O(s^4)  =>
	 		cos = 1 - 1/2*tan^2 + 3/8*tan^4 - 5/16*tan^6 + O(tan^8)

	Below is rational polynomial approximation func cosFromTanP2 for
	abs(tan) < 0.3 by R minimaxApprox package.	Polynomial approximation is
	efficient for 1/sqrt(1+x) and especially when x = tan * tan and tan's
	in this application	are small. Ratio of two polynomials seems more
	efficient than a single	polynomial with same number of terms.
	R:	r$> f <- function(x) 1/sqrt(1+x)
		r$> minimaxApprox(f, 0, 0.09, degree=c(2,2))
*/

// SetGrade calculates sin and cos and dependant forces from tan = grade% /100.
func (c *BikeCalc) SetGrade(tan float64) {
	if useSystemSqrt {
		c.cos = cosFromTanSqrt(tan)
	} else {
		c.cos = cosFromTanP22(tan)
	}
	c.tan = tan
	c.sin = tan * c.cos
	c.setForces()
}

// SetGradeP2 calculates sin and cos and dependant forces from tan = grade% /100.
func (c *BikeCalc) SetGradeP2(tan float64) {
	c.cos = cosFromTanP22(tan)
	if false {
		c.cos = cosFromTanP2NR(tan)
	}
	c.tan = tan
	c.sin = tan * c.cos
	c.setForces()
}

func cosFromTanSqrt(tan float64) (cos float64) {
	cos = 1 / math.Sqrt(1+tan*tan)
	return
}

// cosFromTanP22 returns 1/math.Sqrt(1+tan*tan) by a ratio of two 2. degree
// polynomials. Max error < 6.5e-10 for abs(tan) < 0.3. In benchmark loop
// hardware 1/math.Sqrt(1+tan*tan) is ~75% slower than cosFromTanP22.
// 3 ns vs. 1.7 ns. In cpu profiling the difference is not so clear.
func cosFromTanP22(tan float64) (cos float64) {
	const (
		a1 = 0.73656502
		a2 = 0.05920391
		b1 = 1.2365650
		b2 = 0.3024874
	)
	tan *= tan
	cos = (1 + tan*(a1+tan*a2)) / (1 + tan*(b1+tan*b2))
	return
}

// Max error < 2.4e-10 for abs(tan) < 0.3
func cosFromTanP2NR(tan float64) (cos float64) {
	const (
		a1 = -0.4987452
		a2 = 0.3364923
	)
	tan *= tan
	z := 1 + tan*(a1+tan*a2)
	x := 0.5 * (1 + tan)
	cos = z * (1.5 - x*z*z) // Newton-Raphson iteration
	return
}

// setForces updates all road slope (and mass) dependant forces.
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

func (c *BikeCalc) SetMinPower(w float64) {
	if w < 0 {
		w = 0
	}
	c.minPower = w
}

// func (c *BikeCalc) SetPower(x float64) { c.power = x }

func (c *BikeCalc) SetVelErrors(b bool) { c.calcVelErrors = b }

func (c *BikeCalc) SetVelSolver(i int) {
	c.solverFunc = i
	switch i {
	case newtonRaphson:
		velFromPower = (*BikeCalc).NewtonRaphson

	case newtonHalley:
		velFromPower = (*BikeCalc).NewtonHalley

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
	if c.solverFunc != newtonRaphson && c.solverFunc != newtonHalley {
		c.tolVel = c.bracketLen
	}
}

// SetWeight sets total weight and updates all dependant forces
// and intermediates.
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
// rotatingMass = weight of tyre + tube + rim.
func (c *BikeCalc) SetWeightRotating(kg float64) {
	const rotatingMassReducingFactor = 0.9
	// aproximating 0.9 reduction of the wheel radius, because
	// the mass is not rotating at the outer edge of the wheel.
	if kg > 0 {
		c.massRotating = kg
	}
	c.massKin = c.mass + c.massRotating*rotatingMassReducingFactor
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

func (c *BikeCalc) Sin() float64 { return c.sin }

func (c *BikeCalc) Cos() float64 { return c.cos }

func (c *BikeCalc) Rho() float64 { return c.rho }

func (c *BikeCalc) SolverCalls() int { return c.callsSolver }

func (c *BikeCalc) MaxIter() int { return c.maxIter }

func (c *BikeCalc) SolverRounds() int {
	if c.solverFunc == newtonRaphson || c.solverFunc == newtonHalley {
		return c.iterNR
	}
	return c.callsf
}

func (c *BikeCalc) SolverFunc() int { return c.solverFunc }

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
