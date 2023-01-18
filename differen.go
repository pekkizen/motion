package motion

import "math"

// 1. derivative of speed v with respect to time t: dv/dt = acceleration a = F/m.
// where m is mass and F is sum of resisting/assisting forces.
func (c *BikeCalc) Acceleration(v, power float64) float64 {
	power /= v
	v += c.wind
	return (power - c.fGR - c.cDrag*math.Abs(v)*v) * c.omassKin
}

func (c *BikeCalc) AccelerationPower(v, acceleration float64) (power float64) {
	v += c.wind
	power = v * (c.massKin*acceleration + c.fGR + c.cDrag*math.Abs(v)*v)
	return
}

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
	d2a_dv2 = c.omassKin * (power*ov*ov - cDrag) * 2
	return
}

func (c *BikeCalc) accelerationTimeDerivatives(v, power float64) (a, da_dt, d2a_dt2 float64) {
	var da_dv, d2a_dv2 float64

	a, da_dv, d2a_dv2 = c.accelerationVelocityDerivatives(v, power)
	da_dt = a * da_dv
	d2a_dt2 = a * (a*d2a_dv2 + da_dv*da_dv)
	return
}

func (c *BikeCalc) DeltaTime(Δtime, v0, power float64) (Δvel, Δdist float64) {

	a, da_dt, d2a_dt2 := c.accelerationTimeDerivatives(v0, power)

	Δvel = Δtime * (a + Δtime*(0.5*da_dt+(1.0/6)*Δtime*d2a_dt2))
	Δdist = Δtime * (v0 + 0.5*Δvel)
	return
}

func (c *BikeCalc) velocityDistanceDerivatives(v, power float64) (dv_dx, d2v_dx2, d3v_dx3 float64) {
	ov := 1.0 / v
	a, da_dv, d2a_dv2 := c.accelerationVelocityDerivatives(v, power)

	// 1. derivative of speed v with respect to distance x:  dv/dx = a/v
	dv_dx = a * ov                       //  dv/dx = f
	f1 := (da_dv - dv_dx) * ov           //  df/dv
	f2 := (d2a_dv2 - 2*f1) * ov          //  d2f/dv2
	d2v_dx2 = dv_dx * f1                 //  df/dx		= d2v/dx2
	d3v_dx3 = dv_dx * (dv_dx*f2 + f1*f1) //  d2f/dx2 	= d3v/dx3
	return
}

func (c *BikeCalc) DeltaDist(Δdist, v0, power float64) (Δvel, Δtime float64) {

	dv_dx, d2v_dx2, d3v_dx3 := c.velocityDistanceDerivatives(v0, power)

	Δvel = Δdist * (dv_dx + Δdist*(0.5*d2v_dx2+(1.0/6)*Δdist*d3v_dx3))
	Δtime = Δdist / (v0 + 0.5*Δvel)
	return
}

// Change in kinetic energy ΔKE = distance x force
// fGR 		= gravity and rolling resistance force
// fRider 	= -(power/v0 + power/v1) / 2
// fDrag 	= cDrag * ((v0 + v1)/2 + wind)
// ΔKE 		= 0.5 * mass * (v0^2 - v1^2)
// Δdist	= ΔKE / (fGR + fRider + fDrag)
// Δtime	= Δvel / avg. acceleration

func (c *BikeCalc) DeltaVel(Δvel, v0, power float64) (Δdist, Δtime, fDrag, fSum float64) {
	//this is inlineable. Go 1.19.3 windows/amd64

	vm := v0 + 0.5*Δvel
	v1 := v0 + Δvel
	fSum = -power * vm / (v0 * v1)
	vm += c.wind
	fDrag = c.cDrag * vm * vm
	if vm < 0 {
		fDrag = -fDrag
	}
	fSum += c.fGR + fDrag
	Δdist = 0.5 * c.massKin * (v0*v0 - v1*v1) / fSum
	Δtime = -Δvel * c.massKin / fSum
	return
}

// Solve dist from
// fDrag = cDrag * ((v0 + v1)/2 + wind)
// 0.5 * mass * (v0^2 - v1^2) = dist*(fGR + fBrake + fDrag)

func (c *BikeCalc) DeltaVelBrake(Δvel, v0 float64) (Δdist, fDrag, fSum float64) {

	v1 := v0 + Δvel
	vm := v0 + 0.5*Δvel + c.wind
	fDrag = c.cDrag * vm * vm
	if vm < 0 {
		fDrag = -fDrag
	}
	fSum = c.fGR + fDrag + c.fBrake
	Δdist = 0.5 * c.massKin * (v0*v0 - v1*v1) / fSum
	return
}

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
