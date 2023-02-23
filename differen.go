package motion

import "math"

// Acceleration returns acceleration at speed v and power power.
// 1. derivative of speed v with respect to time t: dv/dt = acceleration a = F/m.
// where m is mass and F is the sum of resisting/assisting forces.
func (c *BikeCalc) Acceleration(v, power float64) float64 {
	power /= v
	v += c.wind
	return (power - c.fGR - c.cDrag*math.Abs(v)*v) * c.omassKin
}

/*
 MeanFdrag calculates
    integral(cDrag * v^2 dv from v0 to v0 + Δvel) / Δvel
 Adjusted to right sign, this is
    v0 += wind
    v1 = v0 + Δvel
    cDrag * 1/3 * (abs(v1^3) - abs(v0^3)) / Δvel
This works for small abs(Δvel) until v0 == v1, and gives zero after that
and NaN for Δvel=0. In practice for abs(Δvel) < ~0.2 the mean speed
air drag is as good and faster.
Without wind we can use
    1/3 * cDrag * (3 * v0 * (v0 + Δvel) + Δvel * Δvel
whicht correctly reduces to cDrag * v0^2, when Δvel = 0.
*/

// MeanFdrag returns exact mean air drag force for speed interval (v0, v0 + Δvel).
func (c *BikeCalc) MeanFdrag(Δvel, v0 float64) (fDrag float64) {

	if c.wind == 0 {
		return (1.0 / 3) * c.cDrag * (3*v0*(v0+Δvel) + Δvel*Δvel)
	}
	x := (1.0 / 3) / Δvel
	v0 += c.wind
	v1 := v0 + Δvel
	if v0 < 0 {
		v0 = -v0
	}
	if v1 < 0 {
		v1 = -v1
	}
	return c.cDrag * (v1*v1*v1 - v0*v0*v0) * x
	// return (1.0 / 3) * c.cDrag * (math.Abs(v1*v1*v1) - math.Abs(v0*v0*v0)) / Δvel
}

/*
MeanFrider return integral (power * 1/v dv from v0 to v0 + Δvel) divided by Δvel.
It should be mean rider force in the interval. This does not increase accuracy
and the log functions cost too much anyway.
*/

func MeanFrider(Δvel, v0, power float64) (fRider float64) {
	fRider = power * (math.Log(v0+Δvel) - math.Log(v0)) / Δvel
	return
}

// accelerationVelocityDerivatives returns acceleration and
// its 1. and 2. derivative with respect to speed.
func (c *BikeCalc) accelerationVelocityDerivatives(v, power float64) (a, da_dv, d2a_dv2 float64) {
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
}

// DeltaTime returns change in speed (m/s) Δvel and distance (m) Δdist for
// driving Δtime seconds with power (w) power and starting with initial 
// speed v0 (m/s).
func (c *BikeCalc) DeltaTime(Δtime, v0, power float64) (Δvel, Δdist float64) {

	a, da_dv, d2a_dv2 := c.accelerationVelocityDerivatives(v0, power)

	// 2. and 3. derivative of velocity with respect to time
	var (
		da_dt   = a * da_dv
		d2a_dt2 = a*(a*d2a_dv2) + da_dt*da_dv
	)
	Δvel = Δtime * (a + Δtime*(0.5*da_dt+(1.0/6)*Δtime*d2a_dt2))
	Δdist = Δtime * (v0 + 0.5*Δvel)
	return
}

// DeltaDist returns change in speed (m/s) Δvel and time (s) Δtime for driving
// Δdist meters with power (w) power and starting with initial speed v0 (m/s).
func (c *BikeCalc) DeltaDist(Δdist, v0, power float64) (Δvel, Δtime float64) {

	a, da_dv, d2a_dv2 := c.accelerationVelocityDerivatives(v0, power)

	// 1. derivative of speed v with respect to distance x: dv/dx = a/v
	// 2. and 3. derivative are calculated below
	var (
		ov      = 1.0 / v0
		dv_dx   = a * ov                     //  dv/dx = f
		f1      = (da_dv - dv_dx) * ov       //  df/dv
		f2      = (d2a_dv2 - 2*f1) * ov      //  d2f/dv2
		d2v_dx2 = dv_dx * f1                 //  df/dx		= d2v/dx2
		d3v_dx3 = dv_dx * (dv_dx*f2 + f1*f1) //  d2f/dx2 	= d3v/dx3
	)
	Δvel = Δdist * (dv_dx + Δdist*(0.5*d2v_dx2+(1.0/6)*Δdist*d3v_dx3))
	Δtime = Δdist / (v0 + 0.5*Δvel)
	return
}

/*
 Solving Δtime and Δdist from v0, v1 = v0 + Δvel and rider power.
 The point: we know the speeds and can estimate speed dependant mean air drag
 force and rider power force in the interval (v0, v1). And then assume linear
 acceleration by the mean acceleration calculated from the estimates.
 Forces below are resisting, negative forces push forward.

 vMean  = (v0 + v1)/2                    - mean speed
 fGR    = fGrav + fRoll                  - gravity + rolling resistance force
 fRider	= -(power/v0 + power/v1) / 2     - mean rider pedaling force.
 vAir   = vMean+wind                     - mean air speed
 fDrag 	= cDrag * Abs(vAir) * vAir       - mean air drag force
 ΔKE    = 0.5 * mass * (v1^2 - v0^2)     - change in kinetic energy
 forces	= fGR + fRider + fDrag           - resisting/assisting forces
 acceleration = -forces / mass           - mean acceleration

 Solve Δdist and Δtime by Δdist from
 ΔKE          = Δdist * -forces or
 v1^2 - v0^2  = 2 * Δdist * acceleration
 Δtime        = Δdist / vMean

 Or solve Δtime and Δdist by Δtime from simpler equivalent formulas
 Δvel  = Δtime * acceleration
 Δdist = Δtime * vMean
*/

// DeltaVel returns distance and time needed to ac/decelerate from initial
// speed (m/s) v0 to speed v0 + Δvel using power (w) power. DeltaVel also returns
// air resistance force fDrag and the sum all resisting/assisting forces fSum.
// fSum can be or must to be used to check if the indented ac/deceleration
// was possible with the given parameters. For acceleration (Δvel > 0) fSum
// must have nonzero negative value and for deceleration (Δvel < 0) it must
// be nonzero positive. When fSum = 0 also accelaration is 0 and no change
// in speed is possible. When v0 = 0 acceleration is goes to infinity and
// DeltaVel returns zero distance and time for any Δvel.
func (c *BikeCalc) DeltaVel(Δvel, v0, power float64) (Δdist, Δtime, fDrag, fSum float64) {

	vMean := v0 + 0.5*Δvel
	fDrag = c.cDrag * signedSquare(vMean+c.wind)
	fSum = c.fGR - power*vMean/(v0*(v0+Δvel)) + fDrag
	Δtime = Δvel * c.massKin / -fSum
	Δdist = Δtime * vMean
	return
}

// DeltaVelMD is same as DeltaVel, but it uses accurate MeanFdrag function
// for air drag force. DeltaVelMD perfoms slightly better than DeltaVel with
// bigger Δvel's but for Δvel < 1 the difference is practically none.
func (c *BikeCalc) DeltaVelMD(Δvel, v0, power float64) (Δdist, Δtime, fDrag, fSum float64) {

	vMean := v0 + 0.5*Δvel
	fDrag = c.MeanFdrag(Δvel, v0)
	fSum = c.fGR - power*vMean/(v0*(v0+Δvel)) + fDrag
	Δtime = Δvel * c.massKin / -fSum
	Δdist = Δtime * vMean
	return
}

// DeltaVelBrake functions are otherwise same as DeltaVel functions,
// but the speed dependant rider force fRider is replaced by constant
// brake force fBrake.

// DeltaVelBrake returns braking distance and time from speed (m/s) v0
// to v0 + Δvel. DeltaVelBrake also returns air resistance force fDrag and
//
//	sum of all resisting/assisting forces fSum. If fSum is near zero or
//
// negative, braking friction coefficient Cbf is too small to brake/slow down.
// The air resistance is calculated from the average speed + wind.
func (c *BikeCalc) DeltaVelBrake(Δvel, v0 float64) (Δdist, Δtime, fDrag, fSum float64) {

	vMean := v0 + 0.5*Δvel
	fDrag = c.cDrag * signedSquare(vMean+c.wind)
	fSum = c.fGR + fDrag + c.fBrake
	Δtime = Δvel * c.massKin / -fSum
	Δdist = Δtime * vMean
	return
}

// DeltaVelBrakeMD returns braking distance and time from speed (m/s) v0
// to v0 + Δvel. The air resistance is calculated exactly by MeanFdrag function.
// Accuracy is slightly better than by DeltaVelBrake but for Δvel < 2 m/s
// there seems to be no practical difference between the methods.
// By a single step DeltaVelBrakeMD can quite accurately solve speed
// differences (Δvel) up to 5+ m/s (18+ km/h).
func (c *BikeCalc) DeltaVelBrakeMD(Δvel, v0 float64) (Δdist, Δtime, fDrag, fSum float64) {

	vMean := v0 + 0.5*Δvel
	fDrag = c.MeanFdrag(Δvel, v0)
	fSum = c.fGR + c.fBrake + fDrag
	Δtime = Δvel * c.massKin / -fSum
	Δdist = Δtime * vMean
	return
}

// Brake returns distance (m), time (s) and wind resistance energy (J) for
// braking from speed v0 to v1. Braking is calculated by steps steps.
// If braking friction coefficientr is too small enough, math.NaN is returned
// for dist and time.
func (c *BikeCalc) Brake(v0, v1 float64, steps int) (dist, time, jouleDrag float64) {
	if v0 <= v1 {
		return
	}
	Δvel := (v1 - v0) / float64(steps)
	for steps > 0 {
		Δdist, Δtime, fDrag, fSum := c.DeltaVelBrakeMD(Δvel, v0)
		if fSum <= 0 {
			dist = math.NaN()
			time = dist
			return
		}
		dist += Δdist
		time += Δtime
		jouleDrag += Δdist * fDrag
		v0 += Δvel
		steps--
	}
	return
}

// Accelerate returns distance (m), time (s) and wind resistance energy (J) for
// accelerating from speed v0 to v1 with power power. Acceleration is
// calculated by steps steps. (v1 - v0) / steps should be < 1.
// If power is not enough, math.NaN is returned for dist and time.
func (c *BikeCalc) Accelerate(v0, v1, power float64, steps int) (dist, time, jouleDrag float64) {
	if v0 >= v1 {
		return
	}
	if c.Acceleration(v1, power) <= 0 {
		dist = math.NaN()
		time = dist
		return
	}
	Δvel := (v1 - v0) / float64(steps)
	for steps > 0 {
		Δdist, Δtime, fDrag, _ := c.DeltaVelMD(Δvel, v0, power)
		dist += Δdist
		time += Δtime
		jouleDrag += Δdist * fDrag
		v0 += Δvel
		steps--
	}
	return
}

// MaxEntryVel returns maximum entry speed for braking to speed
// vExit within distance dist.
func (c *BikeCalc) MaxEntryVel(dist, vExit float64) float64 {
	// Solve vEntry from equation:
	// 0.5*mass*(vEntry^2-vExit^2) = dist*(fGR+fBrake+0.5*cDrag*(vEntry^2+vExit^2))
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
// within distance dist. MaxBrakeStopVel(dist) == MaxEntryVel(dist, 0)
func (c *BikeCalc) MaxBrakeStopVel(dist float64) float64 {
	// Solve vEntry from equation:
	// 0.5 * mass  * vEntry^2 = dist*(fGR + fBrake + 0.5*cDrag*vEntry^2)
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
