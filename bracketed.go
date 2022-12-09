package motion

import "math"

// This is the function we are solving: find v, so that function(v) = 0.
func (c *BikeCalc) function(v float64) float64 {
	c.callsFunc++
	w := v + c.wind
	fD := c.cDrag * w * w
	if w < 0 {
		fD = -fD
	}
	return v*(c.fGR+fD) - c.power
	// return v*(c.fGR+c.cDrag*math.Abs(v+c.wind)*(v+c.wind)) - c.power
}

func (c *BikeCalc) bracketNotFound(s string) {
	c.appendErr(": ")
	c.appendErr(s)
	c.appendErr(": bracket not found before " + maxSpeedS)
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
	before := c.callsFunc
	vL, vR := 0.0, uplim
	for c.function(vR) < 0 {
		if vR > maxSpeed {
			c.appendErr(": Bisect stopped at speed > " + maxSpeedS)
			return vR, 0
		}
		vL = vR
		vR *= 2
	}
	for vR-vL > tol {
		v := (vL + vR) / 2
		if v == vR || v == vL { //adjacent numbers, small tol
			return vL, 1
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
	return (vL + vR) / 2, c.callsFunc - before
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
		if vR += len; vR > maxSpeed {
			c.bracketNotFound("Bracket")
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

// quadratic interpolation of the root.
// Parameter conditions: v0 < v2, f0 <= 0 && 0 <= f2 --> there is a root between v0 and v2.
// Condition -> the discriminant of the square root is always >= 0 (in the paper above).
func (c *BikeCalc) quadratic(v0, f0, v1, f1, v2, f2 float64) float64 {

	d := f2 - f0
	d += math.Sqrt(d*d - 8*f1*(f2+f0-2*f1))
	return v1 - 2*f1*(v2-v0)/d
}

func linear(v0, f0, v1, f1 float64) float64 {
	return (v0*f1 - v1*f0) / (f1 - f0)
}

//SingleLinear --
func (c *BikeCalc) SingleLinear(power, len, vel float64) (float64, int) {
	c.power = power
	before := c.callsFunc
	v0, f0, v1, f1 := c.bracket(len, vel)
	if f0 > 0 || f1 < 0 {
		return v1, 0
	}
	return linear(v0, f0, v1, f1) + len*len*0.014, c.callsFunc - before
}

//DoubleLinear --
func (c *BikeCalc) DoubleLinear(power, len, vel float64) (float64, int) {
	c.power = power
	before := c.callsFunc
	v0, f0, v1, f1 := c.bracket(len, vel)
	if f0 > 0 || f1 < 0 {
		return v1, 0
	}
	vel = linear(v0, f0, v1, f1)
	len *= len * 0.04
	vel += len * 0.5
	vel = linear(c.bracket(len, vel)) + len*len*0.02
	return vel, c.callsFunc - before
}

//SingleQuadratic --
func (c *BikeCalc) SingleQuadratic(power, len, vel float64) (float64, int) {
	c.power = power
	before := c.callsFunc
	v0, f0, v2, f2 := c.bracket(len, vel)
	if f0 > 0 || f2 < 0 {
		return v2, 0
	}
	v1 := (v0 + v2) / 2
	f1 := c.function(v1)
	return c.quadratic(v0, f0, v1, f1, v2, f2), c.callsFunc - before
}

// DoubleQuadratic --
func (c *BikeCalc) DoubleQuadratic(power, len, vel float64) (float64, int) {
	c.power = power
	before := c.callsFunc
	v0, f0, v2, f2 := c.bracket(len, vel)
	v1 := (v0 + v2) / 2
	f1 := c.function(v1)
	vel = c.quadratic(v0, f0, v1, f1, v2, f2)

	len *= len * len * 0.001
	// vel is in (root-len, root+len) with probability >= 95%
	v0, f0, v2, f2 = c.bracket(len, vel)
	v1 = (v0 + v2) / 2
	f1 = c.function(v1)

	vel = c.quadratic(v0, f0, v1, f1, v2, f2)
	return vel, c.callsFunc - before
}
