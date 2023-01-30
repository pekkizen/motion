package motion

import "math"

// This is the function we are solving: find v, so that function(v) = 0.
func (c *BikeCalc) function(v float64) float64 {
	c.callsFun++
	if c.wind == 0 {
		return v*(c.fGR+c.cDrag*v*v) - c.power
	}
	return v*(c.fGR+c.cDrag*math.Abs(v+c.wind)*(v+c.wind)) - c.power
}

// Bisect --
func (c *BikeCalc) Bisect(power, tol, uplim float64) (float64, int) {
	if uplim < 2 {
		uplim = 2
	}
	if tol == 0 {
		tol = c.tolNR * c.tolNR * 0.2
	}
	if tol < 0 {
		tol = 0
	}
	c.power = power
	before := c.callsFun
	vL, vR := 0.0, uplim
	for c.function(vR) < 0 {
		if vR > maxVEL {
			c.appendErr(": Bisect stopped at speed > " + maxVELS)
			return vR, 0
		}
		vL = vR
		vR *= 2
	}
	for vR-vL > tol {
		v := (vL + vR) / 2
		if v == vR || v == vL { //adjacent numbers, small tol
			return vL, c.callsFun - before
		}
		if c.function(v) > 0 {
			vR = v
		} else {
			vL = v
		}
	}
	if c.function(vL) > 0 {
		c.appendErr(" Bisect failed: function(vLeft) > 0.")
		return vL, 0
	}
	return (vL + vR) / 2, c.callsFun - before
}

// bracket return bracket vL < vR, vR - vL <= len and fR <= 0 && 0 <= fL.
// Returning fR < 0 and message bracket not found --> function is negative to very high speeds.
// c.power must be > 0, so that function(0) < 0. Otherwise zero or near zero root may be bracketed.
func (c *BikeCalc) bracket(len, vel float64) (vL, fL, vR, fR float64) {

	vL = vel
	if vL < len {
		vL = len
	}
	fL = c.function(vL)
	for fL > 0 {
		vR, fR = vL, fL
		if vL -= len; vL <= 0 {
			return 0, -c.power, vR, fR
		}
		if fL = c.function(vL); fL <= 0 {
			return
		}
	}
	vR = vL + len
	fR = c.function(vR)
	for fR < 0 {
		vL, fL = vR, fR
		if vR += len; vR > maxVEL {
			c.appendErr(": Bracket not found before " + maxVELS)
			return
		}
		fR = c.function(vR)
	}
	return
}

// BDQRF Bisected Direct Quadratic Regula Falsi
// Applied Mathematical Sciences, Vol. 4, 2010, no. 15, 709 - 718
// Bisected Direct Quadratic Regula Falsi
// Robert G. Gottlieb and Blair F. Thompson
// Odyssey Space Research
// 1120 NASA Parkway, Houston, Texas

// Quadratic interpolation of the root by three velocity points.
// Parameter conditions: v0 < v2, f0 <= 0 && 0 <= f2 --> there is a root between v0 and v2.
// Condition => the discriminant of the square root is always >= 0 (in the paper above).
// v1 = (v0+v2)/2, which makes the simple interpolation formula below possible.
func (c *BikeCalc) quadratic(v0, f0, v2, f2 float64) float64 {

	v1 := (v0 + v2) * 0.5
	f1 := c.function(v1)

	d := f2 - f0
	d += math.Sqrt(d*d - 8*f1*(f2+f0-2*f1))
	return v1 - 2*f1*(v2-v0)/d
}

// Quadratic --
// This performs even with NR, in practical accuray
func (c *BikeCalc) Quadratic(power, len, velGuess float64) (vel float64, callsFun int) {
	c.power = power
	before := c.callsFun
	vL, fL, vR, fR := c.bracket(len, velGuess)
	if fL > 0 || fR < 0 {
		return vR, 0
	}
	vel = c.quadratic(vL, fL, vR, fR)
	callsFun = c.callsFun - before
	return
}

func linear(v0, f0, v1, f1 float64) float64 {
	return (v0*f1 - v1*f0) / (f1 - f0)
}

// SingleLinear --
func (c *BikeCalc) SingleLinear(power, len, vel float64) (float64, int) {
	c.power = power
	before := c.callsFun
	v0, f0, v1, f1 := c.bracket(len, vel)
	if f0 > 0 || f1 < 0 {
		return v1, 0
	}
	return linear(v0, f0, v1, f1) + len*len*0.014, c.callsFun - before
}

// DoubleLinear --
func (c *BikeCalc) DoubleLinear(power, len, vel float64) (float64, int) {
	c.power = power
	before := c.callsFun
	v0, f0, v1, f1 := c.bracket(len, vel)
	if f0 > 0 || f1 < 0 {
		return v1, 0
	}
	vel = linear(v0, f0, v1, f1)
	len *= len * 0.04
	vel += len * 0.5
	vel = linear(c.bracket(len, vel)) + len*len*0.02
	return vel, c.callsFun - before
}

// DoubleQuadratic --
func (c *BikeCalc) DoubleQuadratic(power, len, vel float64) (float64, int) {
	c.power = power
	before := c.callsFun
	v0, f0, v2, f2 := c.bracket(len, vel)
	vel = c.quadratic(v0, f0, v2, f2)

	len *= len * len * 0.001
	// vel is in (root-len, root+len) with probability > 95%
	v0, f0, v2, f2 = c.bracket(len, vel)

	vel = c.quadratic(v0, f0, v2, f2)
	return vel, c.callsFun - before
}
