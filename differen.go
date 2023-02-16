package motion

import "math"

// 1. derivative of speed v with respect to time t: dv/dt = acceleration a = F/m.
// where m is mass and F is the sum of resisting/assisting forces.
func (c *BikeCalc) Acceleration(v, power float64) float64 {
	power /= v
	v += c.wind
	return (power - c.fGR - c.cDrag*math.Abs(v)*v) * c.omassKin
}

// MeanFdrag returns exact mean air drag force between speeds v0 and v0 + Δvel.
// It is integral (cDrag * v^2 dv from v0 to v0 + Δvel) divided by Δvel.
func (c *BikeCalc) MeanFdrag(Δvel, v0 float64) (fDrag float64) {
	ov := 1 / Δvel
	v0 += c.wind
	v1 := v0 + Δvel
	// return (1.0 / 3) * (math.Abs(v1*v1*v1) - math.Abs(v0*v0*v0)) * c.cDrag / Δvel
	if v0 < 0 {
		v0 = -v0
	}
	if v1 < 0 {
		v1 = -v1
	}
	return (1.0 / 3) * ov * c.cDrag * (v1*v1*v1 - v0*v0*v0)
}

// MeanFride return integral (power * 1/v dv from v0 to v0 + Δvel) divided by Δvel.
// It should be mean rider force in the interval. This does not increase accuracy
// and the log functions cost too much anyway.
func MeanFrider(Δvel, v0, power float64) (fRider float64) {
	fRider = power * (math.Log(v0+Δvel) - math.Log(v0)) / Δvel
	return
}

// accelerationVelocityDerivatives returns acceleration and its 1. and 2. derivative
// with respect to speed.
func (c *BikeCalc) accelerationVelocityDerivatives(v, power float64) (a, da_dv, d2a_dv2 float64) {
	ov := 1.0 / v
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

// DeltaTime returns change in speed Δvel and distance Δdist when driving
// Δtime seconds with power power and starting with initial speed v0.
func (c *BikeCalc) DeltaTime(Δtime, v0, power float64) (Δvel, Δdist float64) {

	a, da_dv, d2a_dv2 := c.accelerationVelocityDerivatives(v0, power)

	// 2. and 3. derivative of velocity with respect to time
	da_dt := a * da_dv
	d2a_dt2 := a*(a*d2a_dv2) + da_dt*da_dv

	Δvel = Δtime * (a + Δtime*(0.5*da_dt+(1.0/6)*Δtime*d2a_dt2))
	Δdist = Δtime * (v0 + 0.5*Δvel)
	return
}

// DeltaDist returns change in speed Δvel and time Δtime when driving
// Δdist meters with power power and starting with initial speed v0.
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

// Solving Δtime and Δdist from v0, v1 = v0 + Δvel and rider power.
// The point: we know the speeds and can estimate speed dependant mean air drag
// force and rider power force in the interval (v0, v1).
// And then assume linear acceleration by the mean acceleration calculated
// from the estimates. Forces are resisting, negative forces push forward.
//
// vMean 	= (v0 + v1)/2
// fGR 		= fGrav + fRoll					- gravity + rolling resistance force
// fRider	= -(power/v0 + power/v1) / 2  	- mean rider pedaling force.
// fDrag 	= cDrag *Abs(vMean+wind) *
// 					 	(vMean+wind) 		- mean air drag force
// ΔKE 		= 0.5 * mass * (v1^2 - v0^2)  	- change in kinetic energy
// forces	= fGR + fRider + fDrag
// acceleration = -forces / mass			- mean acceleration
//
// Solve Δdist and Δtime by Δdist from
// ΔKE			= Δdist * -forces or
// v1^2 - v0^2	= 2 * Δdist * acceleration
// Δtime		= Δdist / vMean
//
// Or solve Δtime and Δdist by Δtime from simpler equivalent formulas
// Δvel		= Δtime * acceleration
// Δdist 	= Δtime * vMean

// DeltaVel returns distance and time needed to ac/decelerate from initial
// speed v0 to speed v0 + Δvel using power power. DeltaVel also returns
// air resistance force fDrag and the sum all resisting/assisting forces fSum.
// fSum can be or must to be used to check if the indented ac/deceleration
// was possible with the given parameters. For acceleration fSum must have
// nonzero negative value and for deceleration it must be nonzero positive.
// When fSum = 0 also accelaration is 0 and no change in the speed is possible.
// fSum is denominator in Δdist and Δtime formulas, so near and around zero fSum
// produces very unstable and meaningless Δdist and Δtime, eg. huge negative
// time and distance. When v0 = 0 acceleration is infinity and DeltaVel
// returns zero distance and time for any Δvel.
func (c *BikeCalc) DeltaVel(Δvel, v0, power float64) (Δdist, Δtime, fDrag, fSum float64) {

	vMean := v0 + 0.5*Δvel
	fDrag = c.cDrag * signedSquare(vMean+c.wind)
	fSum = c.fGR - power*vMean/(v0*(v0+Δvel)) + fDrag
	Δtime = Δvel * c.massKin / -fSum
	Δdist = Δtime * vMean
	return
}

// DeltaVelMF is same as DeltaVel, but uses accurate MeanFdrag function
// for air drag force. DeltaVelMF perfoms better than DeltaVel with bigger
// Δvel's but for Δvel < 1 the difference is practically none.
func (c *BikeCalc) DeltaVelMF(Δvel, v0, power float64) (Δdist, Δtime, fDrag, fSum float64) {

	vMean := v0 + 0.5*Δvel
	fDrag = c.MeanFdrag(Δvel, v0)
	fSum = c.fGR - power*vMean/(v0*(v0+Δvel)) + fDrag
	Δtime = Δvel * c.massKin / -fSum
	Δdist = Δtime * vMean
	return
}

// DeltaVelBrake functions are otherwise same as DeltaVel functions,
// but speed dependant rider force fRider is replaced by constant
// brake force fBrake.

// DeltaVelBrake returns braking distance from speed v0 to v0 + Δvel.
// DeltaVelBrake also returns air resistance force fDrag and sum all
// resisting/assisting forces fSum. If fSum is near zero or negative,
// braking friction coefficient Cbf is too small to brake/slow down.
// By a single step DeltaVelBrake can quite accurately solve velocity
// differences (Δvel) up to 5+ m/s (18+ km/h).
func (c *BikeCalc) DeltaVelBrake(Δvel, v0 float64) (Δdist, Δtime, fDrag, fSum float64) {

	vMean := v0 + 0.5*Δvel
	fDrag = c.MeanFdrag(Δvel, v0)
	fSum = c.fGR + c.fBrake + fDrag
	Δtime = Δvel * c.massKin / -fSum
	Δdist = Δtime * vMean
	return
}

// DeltaVelBrakeFast is slightly (2 x) faster version of DeltaVelBrake.
// The air resistance is calculated from the average speed + wind.
// Accuracy is slightly worse but for Δvel < 2 m/s there seems to be no practical
// difference between the methods.
func (c *BikeCalc) DeltaVelBrakeFast(Δvel, v0 float64) (Δdist, Δtime, fDrag, fSum float64) {

	vMean := v0 + 0.5*Δvel
	fDrag = c.cDrag * signedSquare(vMean+c.wind)
	fSum = c.fGR + fDrag + c.fBrake
	Δtime = Δvel * c.massKin / -fSum
	Δdist = Δtime * vMean
	return
}

// Brake returns distance (m), time (s) and wind resistance energy (J) for braking
// from speed v0 to v1. Braking is calculated by steps steps.
func (c *BikeCalc) Brake(v0, v1 float64, steps int) (dist, time, jouleDrag float64) {
	if v0 <= v1 {
		return
	}
	Δvel := (v1 - v0) / float64(steps)
	for steps > 0 {
		Δdist, Δtime, fDrag, fSum := c.DeltaVelBrake(Δvel, v0)
		if fSum < 1e-15 {
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

// MaxEntryVel returns maximum speed for braking to speed vExit within distance dist.
// Solve vEntry from equation:
// 0.5 * mass  * (vEntry^2 - vExit^2) = dist*(fGR + fBrake + 0.5*cDrag*(vEntry^2 + vExit^2))
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

// MaxBrakeStopVel returns maximum speed for braking to stop within distance dist.
// MaxBrakeStopVel(dist) =  MaxEntryVel(dist, 0)
// Solve vEntry from equation:
// 0.5 * mass  * vEntry^2 = dist*(fGR + fBrake + 0.5*cDrag*vEntry^2)
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
