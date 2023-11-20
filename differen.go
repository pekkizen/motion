package motion

import "math"

// Acceleration returns acceleration at speed v and power power.
// 1. derivative of speed v with respect to time t: dv/dt = acceleration a = F/m.
// where m is mass and F = power/v - fGR - cDrag*|v+wind|*(v+wind) is the sum of
// resisting/assisting forces
func (c *BikeCalc) Acceleration(v, power float64) (acce float64) {
	return (power/v - c.fGR - c.cDrag*c.signSq(v)) * c.oMassKin
}

// accelerationVelDerivatives returns acceleration and its 1. and 2. derivative
// with respect to speed. 1. derivative of speed v with respect to
// time t: dv/dt = acceleration a = F/m. At v = 0 and power > 0 acceleration is +inf.
func (c *BikeCalc) accelerationVelDerivatives(v, power float64) (a, da_dv, d2a_dv2 float64) {
	ov := 1 / v
	cDrag := c.cDrag
	if !noWind {
		v += c.wind
		if v < 0 {
			cDrag = -cDrag
		}
	}
	power *= ov
	a = c.oMassKin * (power - c.fGR - cDrag*v*v)
	da_dv = c.oMassKin * (-power*ov - 2*cDrag*v)
	d2a_dv2 = c.oMassKin * 2 * (power*ov*ov - cDrag)
	return
	// can inline (*BikeCalc).accelerationVelDerivatives with cost 73
}

// DeltaTime returns change in speed Δvel and distance Δdist for driving
// Δtime seconds with power power and starting with an initial speed v0.
func (c *BikeCalc) DeltaTime(Δtime, v0, power float64) (Δvel, Δdist float64) {

	a, da_dv, d2a_dv2 := c.accelerationVelDerivatives(v0, power)

	// 2. and 3. derivative of speed with respect to time
	var (
		da_dt   = a * da_dv                 // d2v/dt2
		d2a_dt2 = a*a*d2a_dv2 + da_dt*da_dv // d3v/dt3
	)
	Δvel = Δtime * (a + Δtime*(da_dt*0.5+Δtime*d2a_dt2*(1.0/6)))
	Δdist = Δtime * (v0 + 0.5*Δvel)
	return
}

// DeltaDist returns change in speed Δvel and time Δtime for driving
// Δdist meters with power power and starting with initial speed v0.
func (c *BikeCalc) DeltaDist(Δdist, v0, power float64) (Δvel, Δtime float64) {

	a, da_dv, d2a_dv2 := c.accelerationVelDerivatives(v0, power)

	// 1. derivative of speed v with respect to distance x: dv/dx = a/v
	// 2. and 3. derivative of speed with respect to distance below
	var (
		ov      = 1 / v0
		dv_dx   = a * ov                       // dv/dx = f
		f1      = (da_dv - dv_dx) * ov         // df/dv
		f2      = (d2a_dv2 - 2*f1) * ov        // d2f/dv2
		d2v_dx2 = dv_dx * f1                   // df/dx		= d2v/dx2
		d3v_dx3 = dv_dx * dv_dx * (f2 + f1*f1) // d2f/dx2 	= d3v/dx3
	)
	Δvel = Δdist * (dv_dx + Δdist*(d2v_dx2*0.5+Δdist*d3v_dx3*(1.0/6)))
	Δtime = Δdist / (v0 + 0.5*Δvel)
	return
}

/*
https://physics.stackexchange.com/questions/15587/how-to-get-distance-when-acceleration-is-not-constant

When acceleration is a function of velocity a(v). Then the time as
a function of velocity is
    Δtime  = ∫ 1/a(v) dv from v0 to v1
and the position as a function of velocity is
    Δdist = ∫ v/a(v) dv from v0 to v1

If we have air drag force as a function of velocity fd(v), then the air
drag energy in interval (v0, v1) is
	joule drag = ∫ v/a(v) * fd(v) dv from v0 to v1

Single step solving  Δtime and Δdist for accelerating from v0 to v1 with power power.

vM     = (v0 + v1)/2               - mean speed
fGR    = fGrav + fRoll             - gravity + rolling resistance force
P      =                           - rider power (w)
fRider = (P/v0 + P/v1) / 2         - mean rider pedaling force. Mean of end points
vA     = vM + wind                 - mean air speed
fDrag  = cDrag * abs(vA) * vA      - mean air drag force. Midpoint value
ΔKE    = 0.5*mass*(v1^2 - v0^2)    - change in kinetic energy
forces = fRider - fGR - fDrag      - resisting/assisting forces
acceleration = forces / mass       - "mean" acceleration in (v0, v1)

Calc Δdist and Δtime by
Δdist = ΔKE / forces or
Δdist = 0.5 * (v1^2 - v0^2) / acceleration
Δtime        = Δdist / vM
or
Δtime = Δvel / acceleration
Δdist = Δtime * vM
*/

// DeltaVel returns distance, time and air drag energy for a single step
// ac/decelerating from initial speed v0 to speed v0+Δvel using power power.
// When v0 = 0 or v0+Δvel = 0 acceleration is +inf and DeltaVel
// returns zero distance and time for any Δvel. If returned Δtime is non
// positive ac/decelerating is not possibe by the power given.
func (c *BikeCalc) DeltaVel(Δvel, v0, power float64) (Δdist, Δtime, jDrag float64) {
	var (
		vm     = v0 + 0.5*Δvel
		frider = power * vm / (v0 * (v0 + Δvel)) // rider force, mean of end points
		fdrag  = c.cDrag * c.signSq(vm)          // air drag force, midpoint value
	)
	Δtime = Δvel * c.massKin / (frider - c.fGR - fdrag) // Δvel/acceleration
	Δdist = Δtime * vm
	jDrag = Δdist * fdrag
	return
	// can inline (*BikeCalc).DeltaVel with cost 77 (budget 80)
}

// DeltaVelBrake returns distance, time and air drag energy for braking from
// initial speed v0 to speed v0+Δvel. If returned Δtime <= 0, braking friction
// coefficient Cbf is too small to slow down.
func (c *BikeCalc) DeltaVelBrake(Δvel, v0 float64) (Δdist, Δtime, jDrag float64) {
	var (
		vm    = v0 + 0.5*Δvel
		fdrag = c.cDrag * c.signSq(vm)
	)
	Δtime = Δvel * c.massKin / (-c.fBrake - c.fGR - fdrag)
	Δdist = Δtime * vm
	jDrag = Δdist * fdrag
	return
}

// DeltaVelSimpson returns distance, time and air drag energy for a single step
// ac/decelerating from initial speed v0 to speed v0+Δvel using power power.
// The time, distance and air drag energy integrals are solved numerically by
// Simpson's 1/3 rule.
func (c *BikeCalc) DeltaVelSimpson(Δvel, v0, power float64) (Δdist, Δtime, jDrag float64) {
	v1 := v0 + 0.5*Δvel
	v2 := v0 + Δvel

	d0 := c.cDrag * c.signSq(v0)
	d1 := c.cDrag * c.signSq(v1)
	d2 := c.cDrag * c.signSq(v2)

	a0 := v0 / (power - v0*(c.fGR+d0))
	a1 := 4 * v1 / (power - v1*(c.fGR+d1))
	a2 := v2 / (power - v2*(c.fGR+d2))

	Δvel *= c.massKin * (1.0 / 6)
	Δtime = Δvel * (a0 + a1 + a2)
	Δdist = Δvel * (a0*v0 + a1*v1 + a2*v2)
	jDrag = Δvel * (a0*v0*d0 + a1*v1*d1 + a2*v2*d2)
	return
}

// DeltaVelBoole returns distance, time and air drag energy for a single step
// ac/deceleration from initial speed v0 to speed v0+Δvel using power power.
// The time, distance and air drag energy integrals are solved numerically by Boole's rule.
func (c *BikeCalc) DeltaVelBoole(Δvel, v0, power float64) (Δdist, Δtime, jDrag float64) {
	v1 := v0 + 0.25*Δvel
	v2 := v0 + 0.50*Δvel
	v3 := v0 + 0.75*Δvel
	v4 := v0 + Δvel

	d0 := c.cDrag * c.signSq(v0)
	d1 := c.cDrag * c.signSq(v1)
	d2 := c.cDrag * c.signSq(v2)
	d3 := c.cDrag * c.signSq(v3)
	d4 := c.cDrag * c.signSq(v4)

	a0 := 07 * v0 / (power - v0*(c.fGR+d0))
	a1 := 32 * v1 / (power - v1*(c.fGR+d1))
	a2 := 12 * v2 / (power - v2*(c.fGR+d2))
	a3 := 32 * v3 / (power - v3*(c.fGR+d3))
	a4 := 07 * v4 / (power - v4*(c.fGR+d4))

	Δvel *= c.massKin * (1.0 / 90)
	Δtime = Δvel * (a0 + a1 + a2 + a3 + a4)
	Δdist = Δvel * (a0*v0 + a1*v1 + a2*v2 + a3*v3 + a4*v4)
	jDrag = Δvel * (a0*v0*d0 + a1*v1*d1 + a2*v2*d2 + a3*v3*d3 + a4*v4*d4)
	return
}

func velSteps(v0, v1, Δvel float64) (steps int, h float64) {
	steps = int(math.Abs((v1-v0)/Δvel)) + 1
	h = (v1 - v0) / float64(steps)
	return
}

// Brake returns distance, time and wind resistance energy (J) for
// braking from speed v0 to v1. Braking is calculated by stepsize Δvel.
// If the braking friction coefficient Cbf is too small, -1 is returned
// for all parameters.
func (c *BikeCalc) Brake(v0, v1, Δvel float64) (dist, time, jDrag float64) {
	if v0 <= v1 {
		return
	}
	steps, Δvel := velSteps(v0, v1, Δvel)

	for ; steps > 0; steps-- {
		Δdist, Δtime, drag := c.DeltaVelBrake(Δvel, v0)
		if Δtime <= 0 {
			return -1, -1, -1
		}
		dist += Δdist
		time += Δtime
		jDrag += drag
		v0 += Δvel
	}
	return
}

// AccelerateVelBoole returns distance, time and air resistance energy (J) for
// ac/decelerating from speed v0 to v1 with power power. Acceleration is
// calculated by stepsize Δvel. If power is not enough for acceleration
// or too much for deceleration -1 is returned for all parameters.
// AccelerateVelBoole uses Boole's rule.
func (c *BikeCalc) AccelerateVelBoole(v0, v1, power, Δvel float64) (dist, time, jDrag float64) {
	if v0 == v1 {
		return
	}
	if power == 0 {
		if v1 == 0 {
			v1 = v0 * 1e-14
		}
		if v0 == 0 {
			v0 = 1e-100
		}
	}
	steps, Δvel := velSteps(v0, v1, Δvel)

	if Δvel*c.Acceleration(v1, power) <= 0 {
		return -1, -1, -1
	}
	for ; steps > 0; steps-- {
		Δdist, Δtime, drag := c.DeltaVelBoole(Δvel, v0, power)
		dist += Δdist
		time += Δtime
		jDrag += drag
		v0 += Δvel
	}
	return
}

// AccelerateVel returns distance, time and air resistance energy (J) for
// ac/decelerating from speed v0 to v1 with power power. Acceleration is
// calculated by stepsize Δvel. If power is not enough for acceleration
// or too much for deceleration -1 is returned for all parameters.
// AccelerateVel uses Simpson's 1/3 rule.
func (c *BikeCalc) AccelerateVel(v0, v1, power, Δvel float64) (dist, time, jDrag float64) {
	if v0 == v1 {
		return
	}
	if power == 0 {
		if v1 == 0 {
			v1 = v0 * 1e-14
		}
		if v0 == 0 {
			v0 = 1e-100
		}
	}
	var (
		steps, h = velSteps(v0, v1, Δvel)
		d1       = c.cDrag * c.signSq(v0)
		a1       = 1 / (power/v0 - c.fGR - d1)
		fGR      = c.fGR
		cDrag    = c.cDrag
		a2, d2   float64
	)
	if h*c.Acceleration(v1, power) <= 0 {
		return -1, -1, -1
	}
	time = a1
	dist = a1 * v0
	jDrag = a1 * v0 * d1

	v1 = v0 - 0.5*h
	v2 := v0
	for ; steps > 0; steps-- {
		v1 += h
		v2 += h
		d1 = cDrag * c.signSq(v1)
		d2 = cDrag * c.signSq(v2)

		a1 = 4 * v1 / (power - v1*(fGR+d1))
		a2 = 2 * v2 / (power - v2*(fGR+d2))
		// a1 = 4 / (power/v1 - v1*(fGR+d1))
		// a2 = 2 / (power/v2 - v2*(fGR+d2))

		time += a1 + a2              // integral 1/a(v) dv
		dist += a1*v1 + a2*v2        // integral v/a(v) dv
		jDrag += a1*v1*d1 + a2*v2*d2 // integral v/a(v) * fdrag(v) dv
	}
	a2 *= 0.5
	time -= a2
	dist -= a2 * v2
	jDrag -= a2 * v2 * d2

	h *= c.massKin * (1.0 / 6)
	time *= h
	dist *= h
	jDrag *= h
	return
}

// AccelerateVelMiddlePoint returns distance, time and wind resistance energy (J) for
// ac/decelerating from speed v0 to v1 with power power. Acceleration is
// calculated by stepsize Δvel. If power is not enough for acceleration
// or too much for deceleration -1 is returned for all parameters.
// AccelerateVelMiddlePoint uses Middle Point rule.
func (c *BikeCalc) AccelerateVelMiddlePoint(v0, v1, power, Δvel float64) (dist, time, jDrag float64) {
	if v0 == v1 {
		return
	}
	var (
		steps, h = velSteps(v0, v1, Δvel)
		v        = v0 - 0.5*h
	)
	if h*c.Acceleration(v1, power) <= 0 {
		return -1, -1, -1
	}
	for ; steps > 0; steps-- {
		v += h

		d := c.cDrag * c.signSq(v)
		a := v / (power - v*(c.fGR+d))

		time += a          // integral 1/a(v) dv
		dist += a * v      // integral v/a(v) dv
		jDrag += a * v * d // integral v/a(v) * fdrag(v) dv
	}
	h *= c.massKin
	time *= h
	dist *= h
	jDrag *= h
	return
}

// MaxEntryVel returns approx. maximum entry speed for braking to speed
// vExit within distance dist. For wind = 0.
// Solve vEntry from equation:
// 0.5*mass*(vEntry^2-vExit^2) = dist*(fGR+fBrake+cDrag*0.5*(vEntry^2+vExit^2))
func (c *BikeCalc) MaxEntryVel(dist, vExit float64) float64 {
	var (
		vv   = vExit * vExit
		fSum = c.fGR + c.fBrake + 0.5*c.cDrag*vv
		s    = 0.5 * (c.massKin - dist*c.cDrag)
	)
	if fSum <= 0 || dist <= 0 {
		return vExit
	}
	if s <= 0 {
		return maxVEL
	}
	return math.Sqrt((dist*fSum + 0.5*c.massKin*vv) / s)
}

// MaxBrakeStopVel returns approx. maximum speed for braking to stop
// within distance dist. For wind = 0. Up to 50 meters distances,
// MaxBrakeStopVel gives speeds, for which the exact braking distance
// differs less than 1% from original parameter dist.
// MaxBrakeStopVel(dist) = MaxEntryVel(dist, 0).
// Solves vEntry from equation:
// 0.5 * mass * vEntry^2 = dist*(fGR + fBrake + 0.5*cDrag*vEntry^2)
func (c *BikeCalc) MaxBrakeStopVel(dist float64) float64 {
	var (
		fSum = c.fGR + c.fBrake
		s    = 0.5 * (c.massKin - dist*c.cDrag)
	)
	if fSum <= 0 || dist <= 0 {
		return 0
	}
	if s <= 0 {
		return maxVEL
	}
	return math.Sqrt(dist * fSum / s)
}
