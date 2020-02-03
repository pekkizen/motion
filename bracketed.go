package motion
import "math"

// This is the function we are solving: find v, so that function(v) = 0.
func (c *BikeCalc) function(v float64) float64 {
	c.callsFunc++
	return v*(c.fGR + c.cDrag * abs(v+c.wind)*(v+c.wind)) - c.power
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
		tol =  c.tolNR*c.tolNR*0.2
	}
	if tol < 0 {
		tol = 0
	}
	c.power = power
	before := c.callsFunc
	vL, vR := 0.0, uplim
	for c.function(vR) < 0 {
		if vR > maxSpeed {
			c.appendErr(": Bisect stopped at speed > "+ maxSpeedS)
			return vR, 0
		}
		vL = vR
		vR *= 2
	}
	for vR-vL > tol {
		v := (vL + vR) / 2
		if v == vR || v == vL { //adjacent numbers, small tol
			if abs(c.function(vL)) < abs(c.function(vR)) {
				return vL, 1
			}
			return vR, 1
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

// Bracket return bracket v0 < v1, v1 - v0 <= len and f0 <= 0 && 0 <= f1.
// Returning f0 > 0 || f1 < 0 --> bracket not found --> function is negative to very high speeds.
// c.power must be > 0, so that function(0) < 0. Otherwise zero or near zero root may be bracketed.
func (c *BikeCalc) Bracket(len, vel float64) (v0, f0, v1, f1 float64) {

	v0 = vel
	if v0 < len {
		v0 = len
	}
	f0 = c.function(v0)
	for f0 > 0 {
		v1, f1 =  v0, f0
		if v0 -= len; v0 <= 0 {
			return 0, -c.power, v1, f1
		} 
		if f0 = c.function(v0); f0 <= 0 {
			return
		}
	}
	v1 = v0 + len
	f1 = c.function(v1)
	for f1 < 0 {
		v0, f0 =  v1, f1
		if v1 += len; v1 > maxSpeed {
			c.bracketNotFound("Bracket")
			return
		}
		f1 = c.function(v1)
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
// If f0 <= 0 && 0 <= f2, the discriminant of the square root is always positive (in the article above).
func (c *BikeCalc) quadratic(v0, f0, v2, f2 float64) float64 {
	v1 := (v0 + v2) / 2 
	f1 := c.function(v1)
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
	v0, f0, v1, f1 := c.Bracket(len, vel)
	if f0 > 0 || f1 < 0 {
		return v1, 0
	}
	return linear(v0, f0, v1, f1) + len*len*0.015, c.callsFunc - before 
}

//DoubleLinear --
func (c *BikeCalc) DoubleLinear(power, len, vel float64) (float64, int) {
	c.power = power
	before := c.callsFunc
	v0, f0, v1, f1 := c.Bracket(len, vel)
	if f0 > 0 || f1 < 0 {
		return v1, 0
	}
	vel = linear(v0, f0, v1, f1)
	len *= len * 0.04
	vel += len * 0.5
	vel = linear(c.Bracket(len, vel)) + len*len*0.02
	return vel, c.callsFunc - before 
}

//ThreeLinear --
func (c *BikeCalc) ThreeLinear(power, len, vel float64) (float64, int) {
	c.power = power
	before := c.callsFunc
	v0, f0, v1, f1 := c.Bracket(len, vel)
	if f0 > 0 || f1 < 0 {
		return v1, 0
	}
	vel = linear(v0, f0, v1, f1)
	len *= len * 0.04 
	vel += len * 0.5
	vel = linear(c.Bracket(len, vel))
	len *= len * 0.04 
	vel += len * 0.5
	vel = linear(c.Bracket(len, vel)) + len*len*0.02
	return vel, c.callsFunc - before 
}

//SingleQuadratic --
func (c *BikeCalc) SingleQuadratic(power, len, vel float64) (float64, int) {
	c.power = power
	before := c.callsFunc
	v0, f0, v1, f1 := c.Bracket(len, vel)
	if f0 > 0 || f1 < 0 {
		return v1, 0
	}
	v := c.quadratic(v0, f0, v1, f1)
	return v, c.callsFunc - before 
}

// SingleQuadraticTriple is faster SingleQuadratic, if vel is mostly guessed right.
// Can challenge Newton-Raphson.
// Propably still faster if function evaluations (3 or more) by SIMD/vector instructions.
func (c *BikeCalc) SingleQuadraticTriple(power, len, vel float64) (float64, int) {
	c.power = power
	before := c.callsFunc
	v0 := vel - 0.25*len
	if v0 < len {
		v0 = len
	}
	f0 := c.function(v0)
	if f0 > 0 {
		v0 -= len
		f0 = c.function(v0)
	}
	v1, v2 := v0 + 0.5*len, v0 + len
	f1, f2 := c.function(v1), c.function(v2)

	if f0 <= 0 && 0 <= f2 {
		d := f2 - f0
		d  += math.Sqrt(d*d - 8*f1*(f2+f0-2*f1))
		return v1 - 2*f1*len/d, c.callsFunc - before 
	}
	vel = v2 + len
	if f0 > 0 {
		vel = v0 - len
	}
	vel, i := c.SingleQuadratic(power, len, vel)
	if i == 0 {
		return vel, 0
	}
	return vel, c.callsFunc - before 
}

//DoubleQuadratic --
func (c *BikeCalc) DoubleQuadratic(power, len, vel float64) (float64, int) {
	c.power = power
	before := c.callsFunc
	v0, f0, v1, f1 := c.Bracket(len, vel)
	if f0 > 0 || f1 < 0 {
		return v1, 0
	}
	vel = c.quadratic(v0, f0, v1, f1)
	len *= len * len * 0.001
	// vel on välissä (root-len, root+len) todennäköisyydellä >= 95%
	vel = c.quadratic(c.Bracket(len, vel))
	return vel, c.callsFunc - before 
}
