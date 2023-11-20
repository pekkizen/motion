package motion

import "math"

// This is the function we are solving: find v, so that function(v) = 0.
func (c *BikeCalc) function(v float64) float64 {
	// c.callsf++ // this is expensive and for statistics only
	return v*(c.fGR+c.cDrag*c.signSq(v)) - c.power
}

// Bisect returns the best floating point approximation of the root, with tol = 0.
func (c *BikeCalc) Bisect(power, tol float64) float64 {

	c.power = power
	vL, vR := 0.0, 10.0
	for c.function(vR) < 0 {
		if vR > 1000 {
			c.appendErr("Bisect: no root before 1000 m/s")
			return math.NaN()
		}
		vL = vR
		vR += 10
	}
	for vR-vL > tol {
		v := (vL + vR) / 2
		if v == vR || v == vL { //adjacent numbers
			if math.Abs(c.function(vL)) > math.Abs(c.function(vR)) {
				vL = vR
			}
			return vL
		}
		if c.function(v) > 0 {
			vR = v
		} else {
			vL = v
		}
	}
	return (vL + vR) / 2
}

// bracket returns bracket vL < vR, vR - vL <= len and fL <= 0 && 0 <= fR.
// c.power should be > 0, so that function(0) < 0. Otherwise zero or near zero
// root may be bracketed.
func (c *BikeCalc) bracket(len, vel float64) (vL, fL, vR, fR float64) {
	vL = vel
	if vL < len {
		vL = len
	}
	fL = c.function(vL)
	for fL > 0 {
		vR, fR = vL, fL
		vL -= len
		if vL <= 0 {
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
		vR += len
		if vR > 1000 {
			c.appendErr("bracket: no root before 1000 m/s")
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

// Quadratic interpolation of the root by three equidistance velocity points.
// Parameter conditions: v0 < v2, f0 <= 0 && 0 <= f2 -->
// there is a root between v0 and v2. Condition -> the discriminant
// of the square root is always >= 0 (in the paper above).
func (c *BikeCalc) quadratic(v0, f0, v2, f2 float64) float64 {

	v1 := (v0 + v2) * 0.5
	f1 := c.function(v1)

	d := f2 - f0
	d += math.Sqrt(d*d - 8*f1*(f2+f0-2*f1))
	return v1 - 2*f1*(v2-v0)/d
}

// Linear interpolation of the root by two velocity points.
// Regula falsi interpolation.
func linear(v0, f0, v1, f1 float64) float64 {
	return (v0*f1 - v1*f0) / (f1 - f0)
}

// Quadratic returns speed for power by a single quadratic interpolation.
// First a bracket (v, v+len) which have the root, is searched by function
// bracket. With a good initial value guess, this performs even with NR.
// The returned speeds are highly unbiased, mean signed error near zero.
// error < len^3 * 0.005, ie. cubic convergence.
func (c *BikeCalc) Quadratic(power, len, v float64) (vel float64, callsfun int) {
	if power < 0 {
		return -1, 0
	}
	c.power = power
	if len < minBracketLen {
		len = minBracketLen
	}
	before := c.callsf
	v = c.quadratic(c.bracket(len, v))
	if c.errmsg != nil {
		return v, 0
	}
	return v, c.callsf - before
}

// DoubleQuadratic returns speed for power by two consecutive quadratic
// interpolations.
func (c *BikeCalc) DoubleQuadratic(power, len, v float64) (vel float64, callsfun int) {
	if power < 0 {
		return -1, 0
	}
	if len < minBracketLen {
		len = minBracketLen
	}
	before := c.callsf
	c.power = power
	v = c.quadratic(c.bracket(len, v))
	if c.errmsg != nil {
		return v, 0
	}
	len *= len * len * 0.001 // root is in (v-len, v+len) with probability > 95%
	v = c.quadratic(c.bracket(len, v))
	return v, c.callsf - before
}

// SingleLinear returns speed for power by a single linear interpolation.
func (c *BikeCalc) SingleLinear(power, len, v float64) (vel float64, callsfun int) {
	const bias = 0.014
	if power < 0 {
		return -1, 0
	}
	if len < minBracketLen {
		len = minBracketLen
	}
	before := c.callsf
	c.power = power
	v = linear(c.bracket(len, v)) + len*len*bias
	if c.errmsg != nil {
		return v, 0
	}
	return v, c.callsf - before
}

// DoubleLinear returns speed for power by two consecutive linear
// interpolations.
func (c *BikeCalc) DoubleLinear(power, len, v float64) (vel float64, callsfun int) {
	const bias = 0.014
	if power < 0 {
		return -1, 0
	}
	if len < minBracketLen {
		len = minBracketLen
	}
	before := c.callsf
	c.power = power
	v = linear(c.bracket(len, v)) + len*len*bias
	if c.errmsg != nil {
		return v, 0
	}
	len *= len * 0.025
	v = linear(c.bracket(len, v)) + len*len*bias
	return v, c.callsf - before
}
