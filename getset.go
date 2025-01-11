package motion

import "math"

func (c *BikeCalc) SetBracket(x float64) {
	if x >= 0 && x < minBracketLen {
		x = minBracketLen
	}
	c.bracketLen = x
	if c.solverFunc != NewtonRaphsonMethod &&
		c.solverFunc != NewtonHalleyMethod &&
		c.solverFunc != Householder3Method {
		c.tolVel = c.bracketLen
	}
}

func (c *BikeCalc) SetVelTol(x float64) {
	if x >= 0 && x < minTolNR {
		x = minTolNR
	}
	c.tolVel = x
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
	c.cDrag = 0.5 * c.cdA * c.rho
	if c.cDrag <= 0 {
		c.cDrag = math.NaN()
	}
	c.oCDrag = 1 / c.cDrag
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
	c.SetCdA(0)
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
	c.setForces()
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

	Below is rational polynomial approximation functions cosFromTan for
	abs(tan) < 0.3 by R minimaxApprox package.	Polynomial approximation is
	efficient for 1/sqrt(1+x) and especially when x = tan^2 and tan's
	in this application	are small < 0.3, tan^2 < = 0.09. Ratio of two polynomials
	seems more efficient than a single	polynomial with the same number of terms.

	R:	r$> f <- function(x) 1/sqrt(1+x)
		r$> minimaxApprox(f, 0, 0.09, degree=c(2,2)) // ratio of two 2. degree polynomials
		cos = (1 + tan*(a2+tan*a4)) / (1 + tan*(b2+tan*b4))
		Accurracy is ~3e-10 for abs(tan) < 0.3.
*/

// SetGrade sets road slope and calculates all slope (ans wight) dependant forces and intermediates.
// SetGrade approximates cos by ta by rationof two 2. degree polynomials.s.
func (c *BikeCalc) SetGrade(tan float64) {
	if useExactGrade {
		c.SetGradeExact(tan)
		return
	}
	cos := cosFromTanP22(tan)
	c.tan = tan
	sin := tan * cos
	c.cos = cos
	c.sin = sin
	c.fRoll = cos * c.mgCrr
	c.fGrav = sin * c.mg
	c.fBrake = cos * c.mgCbf
	c.fGR = c.fGrav + c.fRoll
}

// SetGradeExact sets road slope and calculates all slope (ans wight) dependant forces and intermediates.
// SetGradeExact uses full acccuratio, but is skigthly slower than SetGrade.
func (c *BikeCalc) SetGradeExact(tan float64) {
	c.cos = 1 / math.Sqrt(1+tan*tan)
	c.tan = tan
	c.sin = tan * c.cos
	c.setForces()
}

// setForces updates all road slope and mass dependant forces.
func (c *BikeCalc) setForces() {
	c.fRoll = c.cos * c.mgCrr
	c.fGrav = c.sin * c.mg
	c.fBrake = c.cos * c.mgCbf
	c.fGR = c.fGrav + c.fRoll
}

/*
cosFromTanP22 returns inverse square root 1/math.Sqrt(1+tan^2) by a ratio of two
2. degree polynomials of tan^2. Max error ~ 6e-10 for abs(tan) < 0.3.
*/
func cosFromTanP22(tan float64) (cos float64) {
	const (
		a2 = 0.73656502
		a4 = 0.05920391
		b2 = 1.2365650
		b4 = 0.3024874
	)
	tan *= tan
	if useFMA {
		cos = (math.FMA(tan, math.FMA(tan, a4, a2), 1)) /
			(math.FMA(tan, math.FMA(tan, b4, b2), 1))
		return
	}
	cos = (1 + tan*(a2+tan*a4)) / (1 + tan*(b2+tan*b4))
	return
}

/*
cosFromTanP2NR returns 1/math.Sqrt(1+tan^2) by a single 2. degree
polynomial and one Newton-Raphson iteration. This is a kind of "QUAKE" method,
where "magic" number manipulation is replaced by single 2 degree polynomial
approximation. Max error < 2.4e-10 for abs(tan) < 0.3
*/
// func cosFromTanP2NR(tan float64) (cos float64) {
// 	const (
// 		a1 = -0.4987452
// 		a2 = 0.3364923
// 	)
// 	tan *= tan
// 	z := 1 + tan*(a1+tan*a2)
// 	cos = z * (1.5 - 0.5*(1+tan)*z*z) // Newton-Raphson iteration
// 	return
// }

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

func (c *BikeCalc) SetVelErrors(b bool) { c.calcVelErrors = b }

func (c *BikeCalc) SetVelSolver(i int) {
	c.solverFunc = i
	if i < 0 || i > 7 {
		c.solverFunc = Householder3Method
	}
	switch c.solverFunc {
	case Householder3Method: // 3
		velFromPower = (*BikeCalc).Householder3
		if c.tolVel < 0 {
			c.tolVel = 0.8
		}
	case NewtonRaphsonMethod: // 1
		velFromPower = (*BikeCalc).NewtonRaphson
		if c.tolVel < 0 {
			c.tolVel = 0.06
		}
	case NewtonHalleyMethod: // 2
		velFromPower = (*BikeCalc).NewtonHalley
		if c.tolVel < 0 {
			c.tolVel = 0.5
		}

	case singleQuadratic: // 4
		velFromPower = (*BikeCalc).Quadratic
		if c.bracketLen < 0 {
			c.bracketLen = 0.45
		}
	case doubleQuadratic: // 5
		velFromPower = (*BikeCalc).DoubleQuadratic
		if c.bracketLen < 0 {
			c.bracketLen = 3
		}
	case doubleLinear: // 6
		velFromPower = (*BikeCalc).DoubleLinear
		if c.bracketLen < 0 {
			c.bracketLen = 2
		}
	case singleLinear: // 7
		velFromPower = (*BikeCalc).SingleLinear
		if c.bracketLen < 0 {
			c.bracketLen = 0.15
		}
	}
	if c.solverFunc != NewtonRaphsonMethod && c.solverFunc != NewtonHalleyMethod &&
		c.solverFunc != Householder3Method {
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
// rotatingMass = weight of tyre + tube + rim + some of spokes
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
func (c *BikeCalc) SetWind(ms float64) {
	c.wind = ms
	if !useWind {
		c.wind = 0
	}
}

func (c *BikeCalc) AirPressure() float64 { return c.airPressure / 100 }

func (c *BikeCalc) BaseElevation() float64 { return c.baseElevation }

func (c *BikeCalc) Cbf() float64 { return c.cbf }
func (c *BikeCalc) Counter() int { return c.counter }

func (c *BikeCalc) Crr() float64 { return c.crr }

func (c *BikeCalc) CdA() float64 { return c.cdA }

func (c *BikeCalc) Cdrag() float64 { return c.cDrag }

func (c *BikeCalc) Weight() float64 { return c.mass }

func (c *BikeCalc) WeightKin() float64 { return c.massKin }

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
	if c.solverFunc == NewtonRaphsonMethod || c.solverFunc == NewtonHalleyMethod ||
		c.solverFunc == Householder3Method {
		return c.iter
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

func (c *BikeCalc) VelTol() float64 { return c.tolVel }
