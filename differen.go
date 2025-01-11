package motion

import "math"

// Acceleration returns acceleration at speed v and power power.
// 1. derivative of speed v with respect to time t: dv/dt = acceleration a = F/m.
// where m is mass and F = power/v - fGR - cDrag*|v+wind|*(v+wind) is the sum of
// resisting/assisting forces
func (c *BikeCalc) Acceleration(v, power float64) (acce float64) {
	if v == 0 && power == 0 { // ?
		return (-c.fGR - c.cDrag*c.signSq(v)) * c.oMassKin
	}
	return (power/v - (c.fGR + c.cDrag*c.signSq(v))) * c.oMassKin
}

// accelerationVelDeriv returns acceleration and its 1. and 2. derivative
// with respect to speed. 1. derivative of speed v with respect to
// time t: dv/dt = acceleration a = F/m. At v = 0 and power > 0 acceleration is +inf.
func (c *BikeCalc) accelerationVelDeriv(v, power float64) (a, da_dv, d2a_dv2 float64) {
	ov := 1 / v
	cDrag := c.cDrag
	if useWind {
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
func (c *BikeCalc) DeltaTime(Δtime, v0, power float64) (Δvel, Δdist, jdrag float64) {

	a, da_dv, d2a_dv2 := c.accelerationVelDeriv(v0, power)

	// 2. and 3. derivative of speed with respect to time
	var (
		da_dt   = a * da_dv                     // d2v/dt2
		d2a_dt2 = a * (a*d2a_dv2 + da_dv*da_dv) // d3v/dt3
	)

	Δvel = Δtime * (a + Δtime*(da_dt*0.5+Δtime*d2a_dt2*(1.0/6)))
	vm := v0 + 0.5*Δvel
	Δdist = Δtime * vm
	jdrag = Δtime * vm * (c.cDrag * c.signSq(vm))
	return
}

// DeltaDist returns change in speed Δvel and time Δtime for driving
// Δdist meters with power power and starting with initial speed v0.
func (c *BikeCalc) DeltaDist(Δdist, v0, power float64) (Δvel, Δtime, jdrag float64) {

	a, da_dv, d2a_dv2 := c.accelerationVelDeriv(v0, power)

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
	vm := v0 + 0.5*Δvel
	Δtime = Δdist / vm
	jdrag = Δdist * c.cDrag * c.signSq(vm)
	return
}

/*
When acceleration is a function of velocity a(v). Then the time as
a function of velocity is
    Δtime  = ∫ 1/a(v) dv from v0 to v1
and the distance as a function of velocity is
    Δdist = ∫ v/a(v) dv from v0 to v1

If we have air drag force as a function of velocity fd(v), then the air
drag energy in the interval (v0, v1) is
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
Δvel  = v1 - v0
Δtime = Δvel / acceleration
Δdist = Δtime * vM

Alternatives for rider force calculation:
fRider = power * ((math.Log(v0+Δvel) - math.Log(v0)) / Δvel) // exact mean by exact integral
fRider = power * (1.0/6) * (1/v0 + 4/vm + 1/(v0+Δvel))
fRider = power * ((math.Log(v0+Δvel) - math.Log(v0)) / Δvel)    // integral by Simpson's rule
fRider = power / vm // by midpoint value, most error
fDrag = (1.0 / 6) * c.cDrag * (c.signSq(v0) + 4*c.signSq(vm) + c.signSq(v0+Δvel)) // air drag force by Simpson's rule
*/

// DeltaVel returns distance, time and air drag energy for a single step
// ac/decelerating from initial speed v0 to speed v0+Δvel using power power.
// When v0 = 0 or v0+Δvel = 0 acceleration is +inf and DeltaVel
// returns zero distance and time for any Δvel. If returned Δtime is non
// positive or very large positive the ac/decelerating is not possibe by the power given.
func (c *BikeCalc) DeltaVel(Δvel, v0, power float64) (Δdist, Δtime, Δjdrag float64) {
	var (
		vm     = v0 + 0.5*Δvel
		fRider = power * vm / (v0 * (v0 + Δvel)) // rider force, mean of end points
		fDrag  = c.cDrag * c.signSq(vm)          // air drag force, midpoint value
	)
	Δtime = Δvel * c.massKin / (-c.fGR + fRider - fDrag) // Δvel / ~mean acceleration
	Δdist = Δtime * vm
	Δjdrag = Δtime * vm * fDrag
	return
	// can inline, just, one single operation more exceeds limit
}

// DeltaVelBrake returns distance, time and air drag energy for braking from
// initial speed v0 to speed v0+Δvel. If returned Δtime <= 0, braking friction
// coefficient Cbf is too small to slow down.
func (c *BikeCalc) DeltaVelBrake(Δvel, v0 float64) (Δdist, Δtime, jdrag float64) {
	var (
		vm    = v0 + 0.5*Δvel
		fDrag = c.cDrag * c.signSq(vm)
	)
	Δtime = Δvel * c.massKin / (-c.fBrake - c.fGR - fDrag)
	Δdist = Δtime * vm
	jdrag = Δtime * vm * fDrag
	return
	// can inline
}

// DeltaVelDist returns force (braking or riding) and time, force energy and air drag energy
// for decelerating from speed v0 to speed v1 using exactly distance dist.
func (c *BikeCalc) DeltaVelDist(v0, v1, dist float64) (force, time, jForce, jdrag float64) {
	var (
		vm    = (v0 + v1) * 0.5
		fDrag = c.cDrag * c.signSq(vm)
	)
	force = (v0-v1)*c.massKin*vm/dist - c.fGR - fDrag
	time = dist / vm
	jdrag = dist * fDrag
	jForce = dist * force
	return
}

// DeltaVelSimpson returns distance, time and air drag energy for a single step
// ac/decelerating from initial speed v0 to speed v0+Δvel using power power.
// The time, distance and air drag energy integrals are solved numerically by
// Simpson's 1/3 rule. This is a single Simpson step calculation.
func (c *BikeCalc) DeltaVelSimpson(Δvel, v0, power float64) (Δdist, Δtime, jdrag float64) {
	v1 := v0 + 0.5*Δvel
	v2 := v0 + Δvel*0.999 // avoid zero acceleration at end point

	d0 := c.cDrag * c.signSq(v0)
	d1 := c.cDrag * c.signSq(v1)
	d2 := c.cDrag * c.signSq(v2)

	a0 := v0 / (power - v0*(c.fGR+d0))
	a1 := 4 * v1 / (power - v1*(c.fGR+d1))
	a2 := v2 / (power - v2*(c.fGR+d2))

	Δvel *= c.massKin * (1.0 / 6)
	Δtime = Δvel * (a0 + a1 + a2)                   // integral 1/a(v) dv
	Δdist = Δvel * (a0*v0 + a1*v1 + a2*v2)          // integral v/a(v) dv
	jdrag = Δvel * (a0*v0*d0 + a1*v1*d1 + a2*v2*d2) // integral v/a(v) * fdrag(v) dv
	return
}

// DeltaVelBoole returns distance, time and air drag energy for a single (short) step
// ac/deceleration from initial speed v0 to speed v0+Δvel using power power.
// The time, distance and air drag energy integrals are solved numerically by Boole's rule.
// This is a single Boole step calculation.
func (c *BikeCalc) DeltaVelBoole(Δvel, v0, power float64) (Δdist, Δtime, jdrag float64) {
	v1 := v0 + 0.25*Δvel
	v2 := v0 + 0.50*Δvel
	v3 := v0 + 0.75*Δvel
	v4 := v0 + Δvel*0.999 // !!!!!!!!!

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
	Δtime = Δvel * (a0 + a1 + a2 + a3 + a4)                               // integral 1/a(v) dv
	Δdist = Δvel * (a0*v0 + a1*v1 + a2*v2 + a3*v3 + a4*v4)                // integral v/a(v) dv
	jdrag = Δvel * (a0*v0*d0 + a1*v1*d1 + a2*v2*d2 + a3*v3*d3 + a4*v4*d4) // integral v/a(v) * fdrag(v) dv
	return
}

func velSteps(v0, v1, Δvel float64) (steps int, ΔvelExact float64) {
	steps = int(math.Abs((v1-v0)/Δvel)) + 1
	ΔvelExact = (v1 - v0) / float64(steps)
	return
}

// Brake returns distance, time and wind resistance energy (J) for
// braking from speed v0 to v1. Braking is calculated by stepsize Δvel.
// If the braking friction coefficient Cbf is too small, -1 is returned
// for all parameters.
func (c *BikeCalc) Brake(v0, v1, Δvel float64) (dist, time, jdrag float64) {
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
		jdrag += drag
		v0 += Δvel
	}
	return
}

// AccelerateVelBoole returns distance, time and air resistance energy (J) for
// ac/decelerating from speed v0 to v1 with power power. Acceleration is
// calculated by stepsize Δvel. If power is not enough for acceleration
// or too much for deceleration -1 is returned for all parameters.
// AccelerateVelBoole uses Boole's rule.
func (c *BikeCalc) AccelerateVelBoole(v0, v1, power, Δvel float64) (dist, time, jdrag float64) {
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

	if Δvel*c.Acceleration(v1, power) < 0 {
		return -1, -1, -1
	}
	for ; steps > 0; steps-- {
		Δdist, Δtime, drag := c.DeltaVelBoole(Δvel, v0, power)
		dist += Δdist
		time += Δtime
		jdrag += drag
		v0 += Δvel
	}
	return
}

// AccelerateVel returns distance, time and air resistance energy (J) for
// ac/decelerating from speed v0 to v1 with power power. Acceleration is
// calculated by stepsize Δvel. If power is not enough for acceleration
// or too much for deceleration -1 is returned for all parameters.
// AccelerateVel uses Simpson's 1/3 rule. AccelerateVel is the preferred/fastest
// function for calculating accelerations by velocity steps.
func (c *BikeCalc) AccelerateVel(v0, v1, power, Δvel float64) (dist, time, jdrag float64) {
	if v0 == v1 {
		return // 0, 0, 0
	}
	if power == 0 { // avoid NaNs of division by zero
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
	if h*c.Acceleration(v1, power) < 0 {
		return -1, -1, -1
	}
	time = a1
	dist = a1 * v0
	jdrag = a1 * v0 * d1

	v1 = v0 - 0.5*h
	v2 := v0
	for ; steps > 0; steps-- {
		v1 += h
		v2 += h
		d1 = cDrag * c.signSq(v1)
		d2 = cDrag * c.signSq(v2)

		a1 = 4 * v1 / (power - v1*(fGR+d1))
		a2 = 2 * v2 / (power - v2*(fGR+d2))

		time += a1 + a2              // integral 1/a(v) dv
		dist += a1*v1 + a2*v2        // integral v/a(v) dv
		jdrag += a1*v1*d1 + a2*v2*d2 // integral v/a(v) * fdrag(v) dv
	}
	a2 *= 0.5
	time -= a2
	dist -= a2 * v2
	jdrag -= a2 * v2 * d2

	h *= c.massKin * (1.0 / 6)
	time *= h
	dist *= h
	jdrag *= h
	return
}

// AccelerateVelMiddlePoint returns distance, time and wind resistance energy (J) for
// ac/decelerating from speed v0 to v1 with power power. Acceleration is
// calculated by stepsize Δvel. If power is not enough for acceleration
// or too much for deceleration -1 is returned for all parameters.
// AccelerateVelMiddlePoint uses Middle Point rule. For small Δvel, AccelerateVelMiddlePoint is a
// proof/ckeck method for other functions.
func (c *BikeCalc) AccelerateVelMiddlePoint(v0, v1, power, Δvel float64) (dist, time, jdrag float64) {
	if v0 == v1 {
		return
	}
	var (
		steps, h = velSteps(v0, v1, Δvel)
		v        = v0 - 0.5*h
	)
	if h*c.Acceleration(v1-0.5*h, power) < 0 {
		return -1, -1, -1
	}
	for ; steps > 0; steps-- {
		v += h

		d := c.cDrag * c.signSq(v)
		a := v / (power - v*(c.fGR+d))

		time += a
		dist += a * v      // integral v/a(v) dv
		jdrag += a * v * d // integral v/a(v) * fdrag(v) dv
	}
	h *= c.massKin
	time *= h
	dist *= h
	jdrag *= h
	return
}

// MaxEntryVelNoWind returns approx. maximum entry speed for braking to speed
// vExit within distance dist. For wind = 0.
// Solve vEntry from equation:
// 0.5*mass*(vEntry^2-vExit^2) = dist*(fGR+fBrake+cDrag*0.5*(vEntry^2+vExit^2))
func (c *BikeCalc) MaxEntryVelNoWind(dist, vExit float64) (vel float64) {
	var (
		s    = 2 / (c.massKin - dist*c.cDrag)
		fSum = c.fGR + c.fBrake + 0.5*c.cDrag*(vExit*vExit)
	)
	if fSum <= 0 || dist <= 0 {
		return vExit
	}
	if s <= 0 {
		return maxVel
	}
	return math.Sqrt((dist*fSum + 0.5*c.massKin*(vExit*vExit)) * s)
}

// MaxEntryVelWind solves v0 from equation:
// 0.5 * mass * (v0^2 - v1^2)= dist*(fGR + fBrake + cDrag*0.5*(v0+wind)^2+(v1 + wind)^2)
// Problem is that we cannot change cDrag's sign for negative air speed. Still this
// generally gives much better approximation than MaxEntryVelNoWind.
func (c *BikeCalc) MaxEntryVelWind(dist, vExit float64) (vEntry float64) {
	var (
		md = c.massKin / dist
		z  = -c.cDrag*c.wind*(vExit+c.wind) - c.fGR - c.fBrake
		b  = -c.cDrag * c.wind
		a  = md - c.cDrag
	)
	z = b*b - 2*a*(z+0.5*(vExit*vExit)*(-md-c.cDrag))

	if z <= 0 || dist <= 0 {
		return vExit
	}
	if a <= 0 {
		return maxVel
	}
	return (-b + math.Sqrt(z)) / a
}

func (c *BikeCalc) MaxEntryVel(dist, vExit float64) (vEntry float64) {
	if !useWind {
		return c.MaxEntryVelNoWind(dist, vExit)
	}
	return c.MaxEntryVelWind(dist, vExit)
}

// MaxEntryVelBrakeIntegral returns maximum entry speed for braking to speed
// vExit within distance brakeDist. For any wind. For small Δvel ~ 0.5, this is
// the "exact" function for the purpose. And quite slow compred to other MaxEntryVel function.
func (c *BikeCalc) MaxEntryVelBrakeIntegral(brakeDist, Δvel, vExit float64) (vEntry float64) {
	if Δvel > 0 {
		Δvel = -Δvel
	}
	if vExit < -Δvel {
		vExit = -Δvel
	}
	dist := 0.0
	for {
		vExit -= Δvel
		Δdist, _, _ := c.DeltaVelBrake(Δvel, vExit)
		if Δdist <= 0 { // brakig not possible
			return vExit // or NaN ?
		}
		dist += Δdist
		if dist >= brakeDist { // interpolate speed at brakeDist
			vExit += Δvel * (dist - brakeDist) / Δdist
			break
		}
	}
	return vExit
}

// MaxBrakeStopVelNoWind returns approx. maximum speed for braking to stop
// within distance dist. For wind = 0. Up to 50 meters distances,
// MaxBrakeStopVelNoWind gives speeds, for which the exact braking distance
// differs less than 1% from original parameter dist.
// MaxBrakeStopVelNoWind(dist) = MaxEntryVel(dist, 0).
// Solves vEntry from equation:
// 0.5 * mass * vEntry^2 = dist*(fGR + fBrake + 0.5*cDrag*vEntry^2)
func (c *BikeCalc) MaxBrakeStopVelNoWind(dist float64) float64 {
	var (
		d = 2 / (c.massKin - dist*c.cDrag)
		f = c.fGR + c.fBrake
	)
	if f <= 0 || dist <= 0 {
		return 0
	}
	if d <= 0 {
		return maxVel
	}
	return math.Sqrt(dist * f * d)
}

// MaxEntryVelWind solves v from equation:
// 0.5 * mass * v^2 = dist*(fGR + fBrake + 0.5*cDrag*(v+wind)^2
func (c *BikeCalc) MaxBrakeStopVelWind(dist float64) float64 {
	var (
		s = -dist * c.cDrag
		b = s * c.wind
		z = -dist*(c.fGR+c.fBrake) + b*c.wind
		a = c.massKin + s
	)
	if z >= 0 || dist <= 0 {
		return 0
	}
	d := b*b - 2*a*z
	if d <= 0 {
		return maxVel
	}
	return (-b + math.Sqrt(d)) / a
}

// Only either one of the conditons c.wind == 0 or c.wind != 0 is always true.
// In some applications calculating a bicycle route, at least. So, branch prediction
// should work perfectly. The two functions wrapped this way are 2.5 x faster
// than than the Wind - functons alone. In no wind conditions
func (c *BikeCalc) MaxBrakeStopVel(dist float64) float64 {
	if c.wind == 0 {
		return c.MaxBrakeStopVelNoWind(dist)
	}
	return c.MaxBrakeStopVelWind(dist)
}

// BrakeForce returns braking force for braking from speed v0 to speed v1
// using exactly distance dist.
func (c *BikeCalc) BrakeForce(v0, v1, dist float64) (fBrake float64) {

	vm := (v0 + v1) * 0.5
	fBrake = (v0-v1)*c.massKin*vm/dist - c.fGR - c.cDrag*c.signSq(vm)
	return
}
