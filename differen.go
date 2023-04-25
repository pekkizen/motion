package motion

import "math"

// Acceleration returns acceleration at speed v and power power.
// 1. derivative of speed v with respect to time t: dv/dt = acceleration a = F/m.
// where m is mass and F is the sum of resisting/assisting forces.
func (c *BikeCalc) Acceleration(v, power float64) (acce float64) {
	return (power/v - c.fGR - c.cDrag*c.signSq(v)) * c.omassKin
}

// func (c *BikeCalc) noAcceleration(v, Δvel, power float64) bool {
// 	return Δvel*c.Acceleration(v, power) <= 0
// 	// return Δvel/math.Abs(Δvel)*c.Acceleration(v, power) <= 0
// }

// accelerationVelDerivatives returns acceleration and its 1. and 2. derivative
// with respect to speed.
func (c *BikeCalc) accelerationVelDerivatives(v, power float64) (a, da_dv, d2a_dv2 float64) {
	ov := 1 / v
	v += c.wind
	cDrag := c.cDrag
	if v < 0 {
		cDrag = -cDrag
	}
	power *= ov
	a = c.omassKin * (power - c.fGR - cDrag*v*v)
	da_dv = c.omassKin * (-power*ov - 2*cDrag*v)
	d2a_dv2 = c.omassKin * 2 * (power*ov*ov - cDrag)
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
	// 2. and 3. derivative of speed with respect to distance 
	var (
		ov      = 1.0 / v0
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
    t1 - t0  = ∫ 1/a(v) dv from v0 to v1    (5)
and the position as a function of velocity is
    x1 - x0 = ∫ v/a(v) dv from v0 to v1     (6)

If we have air drag force as a function of v fd(v), then the air drag energy in
interval (v0, v1) is
	joule drag = ∫ v/a(v) * fd(v) dv from v0 to v1
I have no proper derivation of this, but in the numerical integrals this formula gives right numbers.

Single step solving  Δtime and Δdist for driving from v0 to v1 with power power.

vM     = (v0 + v1)/2               - mean speed
fGR    = fGrav + fRoll             - gravity + rolling resistance force
P      =                           - rider power (w)
fRider = (P/v0 + P/v1) / 2         - mean rider pedaling force. Mean of end points
vA     = vM + wind                 - mean air speed
fDrag  = cDrag * abs(vA) * vA      - mean air drag force. Midpoint value
ΔKE    = 0.5*mass*(v1^2 - v0^2)    - change in kinetic energy
forces = fRider - fGR - fDrag      - resisting/assisting forces
acceleration = forces / mass       - "mean" acceleration in (v0, v1)

Solve Δdist and Δtime from
ΔKE          = Δdist * forces or
v1^2 - v0^2  = 2 * Δdist * acceleration
Δtime        = Δdist / vM

 Or solve Δtime and Δdist from simpler equivalent formulas
 Δvel  = Δtime * acceleration
 Δdist = Δtime * vM
*/

// DeltaVel returns distance, time and air drag energy for  a single step
// ac/decelerating from initial speed v0 to speed v0+Δvel using power power.
// When v0 = 0 or v0+Δvel = 0 acceleration goes to infinity and DeltaVel
// returns zero distance and time for any Δvel. If returned Δtime is non
// positive ac/decelerating is not possibe by the power given.
func (c *BikeCalc) DeltaVel(Δvel, v0, power float64) (Δdist, Δtime, jDrag float64) {
	var (
		vm     = v0 + 0.5*Δvel
		frider = power * vm / (v0 * (v0 + Δvel)) // rider force, mean of end points
		fdrag  = c.cDrag * c.signSq(vm)          // drag force, midpoint value
	)
	Δtime = Δvel * c.massKin / (frider - c.fGR - fdrag)
	Δdist = Δtime * vm
	jDrag = Δdist * fdrag
	return
	// can inline (*BikeCalc).DeltaVel with cost 77 (budget 80)
}

// DeltaVelSimpson returns distance, time and air drag energy for a single step
// ac/decelerating from initial speed v0 to speed v0+Δvel using power power.
// The time, distance and air drag integrals are solved  numerically by
// Simpson's 1/3 rule.
func (c *BikeCalc) DeltaVelSimpson(Δvel, v0, power float64) (Δdist, Δtime, jDrag float64) {
	v1 := v0 + 0.5*Δvel
	v2 := v0 + Δvel

	d0 := c.cDrag * c.signSq(v0)
	d1 := c.cDrag * c.signSq(v1)
	d2 := c.cDrag * c.signSq(v2)

	a0 := 1 / (power/v0 - c.fGR - d0)
	a1 := 4 / (power/v1 - c.fGR - d1)
	a2 := 1 / (power/v2 - c.fGR - d2)

	Δvel *= c.massKin * (1.0 / 6)
	Δtime = Δvel * (a0 + a1 + a2)
	Δdist = Δvel * (a0*v0 + a1*v1 + a2*v2)
	jDrag = Δvel * (a0*v0*d0 + a1*v1*d1 + a2*v2*d2)
	return
}

// DeltaVelBoole returns distance, time and air drag energy for a single step
// ac/deceleration from initial speed v0 to speed v0+Δvel using power power.
// The time and distance integrals are solved  numerically by by Boole's rule.
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

	// a0 := 7 / (power/v0 - c.fGR - d0)
	// a1 := 32 / (power/v1 - c.fGR - d1)
	// a2 := 12 / (power/v2 - c.fGR - d2)
	// a3 := 32 / (power/v3 - c.fGR - d3)
	// a4 := 7 / (power/v4 - c.fGR - d4)
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

// DeltaVelBrake returns distance, time and air drag energy for braking from
// initial speed v0 to speed v0+Δvel. If returned Δtime <= 0, braking friction
// coefficient Cbf is too small to brake/slow down.
func (c *BikeCalc) DeltaVelBrake(Δvel, v0 float64) (Δdist, Δtime, jDrag float64) {

	vm := v0 + 0.5*Δvel
	fdrag := c.cDrag * c.signSq(vm)
	Δtime = Δvel * c.massKin / (-c.fBrake - c.fGR - fdrag)
	Δdist = Δtime * vm
	jDrag = Δdist * fdrag
	return
}

func velSteps(v0, v1, Δvel float64) (steps int, ΔvelOut float64) {
	steps = int(math.Abs((v1-v0)/Δvel)) + 1
	ΔvelOut = (v1 - v0) / float64(steps)
	return
}

// Brake returns distance, time and wind resistance energy (J) for
// braking from speed v0 to v1. Braking is calculated by stepsize Δvel.
// If the braking friction coefficient Cbf is too small, -1 is returned
// for all parameters.
func (c *BikeCalc) Brake(v0, v1, Δvel float64) (dist, time, drag float64) {
	if v0 <= v1 {
		return
	}
	steps, Δvel := velSteps(v0, v1, Δvel)

	for ; steps > 0; steps-- {
		Δdist, Δtime, jDrag := c.DeltaVelBrake(Δvel, v0)
		if Δtime <= 0 {
			dist, time, drag = -1, -1, -1
			return
		}
		dist += Δdist
		time += Δtime
		drag += jDrag
		v0 += Δvel
	}
	return
}

// AccelerateBoole returns distance, time and wind resistance energy (J) for
// ac/decelerating from speed v0 to v1 with power power. Acceleration is
// calculated by stepsize Δvel. If power is not enough for acceleration
// or too much for deceleration -1 is returned for all parameters.
// AccelerateBoole uses Boole's rule.
func (c *BikeCalc) AccelerateBoole(v0, v1, power, Δvel float64) (dist, time, drag float64) {
	if v0 == v1 {
		return
	}
	if power == 0 {
		if v1 == 0 {
			v1 = v0 * 1e-14
		} else if v0 == 0 {
			v0 = 1e-100
		}
	}
	steps, Δvel := velSteps(v0, v1, Δvel)

	if Δvel*c.Acceleration(v1, power) <= 0 {
		dist, time, drag = -1, -1, -1
		return
	}
	for ; steps > 0; steps-- {
		Δdist, Δtime, jdrag := c.DeltaVelBoole(Δvel, v0, power)
		// Δdist, Δtime, jdrag := c.DeltaVelSimpson(Δvel, v0, power)
		dist += Δdist
		time += Δtime
		drag += jdrag
		v0 += Δvel
	}
	return
}

// Accelerate returns distance, time and wind resistance energy (J) for
// ac/decelerating from speed v0 to v1 with power power. Acceleration is
// calculated by stepsize Δvel. If power is not enough for acceleration
// or too much for deceleration -1 is returned for all parameters.
// Accelerate uses Simpson's 1/3 rule.
func (c *BikeCalc) Accelerate(v0, v1, power, Δvel float64) (dist, time, drag float64) {
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
		a2, d2   float64
	)
	if h*c.Acceleration(v1, power) <= 0 {
		dist, time, drag = -1, -1, -1
		return
	}
	time = a1
	dist = a1 * v0
	drag = a1 * v0 * d1

	v1 = v0 - 0.5*h
	v2 := v0
	for ; steps > 0; steps-- {
		v1 += h
		v2 += h
		d1 = c.cDrag * c.signSq(v1)
		d2 = c.cDrag * c.signSq(v2)

		// a1 = 4 / (power/v1 - c.fGR - d1)
		// a2 = 2 / (power/v2 - c.fGR - d2)
		a1 = 4 * v1 / (power - v1*(c.fGR+d1))
		a2 = 2 * v2 / (power - v2*(c.fGR+d2))

		time += a1 + a2             // integral 1/a(v) dv
		dist += a1*v1 + a2*v2       // integral v/a(v) dv
		drag += a1*v1*d1 + a2*v2*d2 // integral v/a(v) * fdrag(v) dv
	}
	a2 /= 2
	time -= a2
	dist -= a2 * v2
	drag -= a2 * v2 * d2
	h *= c.massKin * (1.0 / 6)
	time *= h
	dist *= h
	drag *= h
	return
}

// AccelerateMiddlePoint returns distance, time and wind resistance energy (J) for
// ac/decelerating from speed v0 to v1 with power power. Acceleration is
// calculated by stepsize Δvel. If power is not enough for acceleration
// or too much for deceleration -1 is returned for all parameters.
// AccelerateMiddlePoint uses Middle Point rule.
func (c *BikeCalc) AccelerateMiddlePoint(v0, v1, power, Δvel float64) (dist, time, drag float64) {
	if v0 == v1 {
		return
	}
	var (
		steps, h = velSteps(v0, v1, Δvel)
		v        = v0 - 0.5*h
	)
	if h*c.Acceleration(v1, power) <= 0 {
		dist, time, drag = -1, -1, -1
		return
	}
	for ; steps > 0; steps-- {
		v += h

		d := c.cDrag * c.signSq(v)
		a := v / (power - v*(c.fGR+d))
		// a := 1 / (power/v - c.fGR-d)

		time += a         // integral 1/a(v) dv
		dist += a * v     // integral v/a(v) dv
		drag += a * v * d // integral v/a(v) * fdrag(v) dv
	}
	h *= c.massKin
	time *= h
	dist *= h
	drag *= h
	return
}

// MaxEntryVel returns maximum entry speed for braking to speed
// vExit within distance dist. For wind = 0.
// Solves vEntry from equation:
// mass*(vEntry^2-vExit^2)/2 = dist*(fGR+fBrake+cDrag*(vEntry^2+vExit^2))/2
func (c *BikeCalc) MaxEntryVel(dist, vExit float64) float64 {
	vv := vExit * vExit
	fSum := c.fGR + c.fBrake + 0.5*c.cDrag*vv
	if fSum <= 0 || dist <= 0 {
		return vExit
	}
	s := 0.5 * (c.massKin - dist*c.cDrag)
	if s <= 0 {
		return maxVEL
	}
	return math.Sqrt((dist*fSum + 0.5*c.massKin*vv) / s)
}

// MaxBrakeStopVel returns maximum speed for braking to stop
// within distance dist. For wind = 0.
// MaxBrakeStopVel(dist) == MaxEntryVel(dist, 0)
// Solves vEntry from equation:
// 0.5 * mass * vEntry^2 = dist*(fGR + fBrake + 0.5*cDrag*vEntry^2)
func (c *BikeCalc) MaxBrakeStopVel(dist float64) float64 {
	fSum := c.fGR + c.fBrake
	if fSum <= 0 || dist <= 0 {
		return 0
	}
	s := 0.5 * (c.massKin - dist*c.cDrag)
	if s <= 0 {
		return maxVEL
	}
	return math.Sqrt(dist * fSum / s)
}
