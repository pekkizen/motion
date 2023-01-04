package motion

import "math"

// 1. derivative of speed v with respect to time t: dv/dt = acceleration a = F/m.
// where m is mass and F is sum of resisting/assisting forces.
func (c *BikeCalc) Acceleration(v, power float64) float64 {
	power /= v
	v += c.wind
	return (power - c.fGR - c.cDrag*math.Abs(v)*v) * c.oMassKin
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
	a = c.oMassKin * (power - c.fGR - cDrag*v*v)
	da_dv = c.oMassKin * (-power*ov - 2*cDrag*v)
	d2a_dv2 = c.oMassKin * (power*ov*ov - cDrag) * 2
	return
}

func (c *BikeCalc) VelocityTimeDerivatives(v, power float64) (a, da_dt, d2a_dt2 float64) {
	var da_dv, d2a_dv2 float64

	a, da_dv, d2a_dv2 = c.accelerationVelocityDerivatives(v, power)
	da_dt = a * da_dv
	d2a_dt2 = a * (a*d2a_dv2 + da_dv*da_dv)
	return
}

func (c *BikeCalc) DeltaTime(Δtime, v0, power float64) (Δvel, Δdist float64) {

	a, da_dt, d2a_dt2 := c.VelocityTimeDerivatives(v0, power)

	Δvel = Δtime * (a + Δtime*(0.5*da_dt+(1.0/6)*Δtime*d2a_dt2))
	Δdist = Δtime * (v0 + 0.5*Δvel)
	return
}

func (c *BikeCalc) VelocityDistanceDerivatives(v, power float64) (dv_dx, d2v_dx2, d3v_dx3 float64) {
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

	dv_dx, d2v_dx2, d3v_dx3 := c.VelocityDistanceDerivatives(v0, power)

	// Δvel = Δdist * (dv_dx + Δdist*(0.5*d2v_dx2+(1.0/6)*Δdist*d3v_dx3))
	Δvel = Δdist * (dv_dx + Δdist*(0.5*d2v_dx2+(1.0/6)*Δdist*d3v_dx3))
	Δtime = Δdist / (v0 + 0.5*Δvel)
	return
}

// Change in kinetic energy ΔKE = distance x force
// fGR 		= gravity and rolling resistance force
// fRider 	= -(power/v0 + power/v1) / 2
// fDrag 	= Fdrag((v0 + v1)/2 + wind)
// ΔKE 		= 0.5 * mass * (v0^2 - v1^2)
// Δdist	= ΔKE / (fGR + fRider + fDrag)
// Δtime	= Δvel / avg. acceleration
func (c *BikeCalc) DeltaVel(Δvel, v0, power float64) (Δdist, Δtime, fDrag, fSum float64) {

	vm := v0 + 0.5*Δvel
	v1 := v0 + Δvel
	vAir := vm + c.wind
	Δdist = 0.5 * c.massKin * (v0*v0 - v1*v1)
	Δtime = -Δvel * c.massKin
	fRider := -power * vm / (v0 * v1)
	fDrag = c.cDrag * vAir * math.Abs(vAir)
	fSum = c.fGR + fRider + fDrag
	Δdist /= fSum
	Δtime /= fSum
	return
}

// DistBrake(MaxEntryVel(dist, vExit), vExit) == dist
// DistBrake(VelFromBrakeDist(dist), 0) == dist
// Funktioilla VelFromBrakeDist ja MaxEntryVel lasketaan nopeusrajoituksia tyynessä.
//
// fDrag = 0.5*cDrag*((vEntry+wind)^2 + (vExit+wind)^2)
// 0.5 * mass * (vEntry^2 - vExit^2) = dist*(fGR + fBrake + fDrag)
func (c *BikeCalc) DistBrake(v0, v1 float64) (Δdist, fDrag, fSum float64) {

	ΔKE := 0.5 * c.massKin * (v0*v0 - v1*v1)
	v0 += c.wind
	v1 += c.wind
	fDrag = c.cDrag * 0.5 * (v0*math.Abs(v0) + v1*math.Abs(v1))
	fSum = c.fGR + fDrag + c.fBrake
	Δdist = ΔKE / fSum
	return
}

// MaxEntryVel laskee maksiminopeuden, josta etäisyydellä dist on
// mahdollista jarruttaa nopeuteen vExit (tuuli=0).
// Solve vEntry from equation:
// 0.5 * mass  * (vEntry^2 - vExit^2) = dist*(fGR + fBrake + 0.5*cDrag*(vEntry^2 + vExit^2))
func (c *BikeCalc) MaxEntryVel(dist, vExit float64) float64 {

	vv := vExit * vExit
	fSum := c.fGR + c.fBrake + 0.5*c.cDrag*vv
	if fSum <= 0 {
		return vExit
	}
	s := 0.5 * (c.massKin - dist*c.cDrag)
	if s <= 0 {
		return maxVEL
	}
	return math.Sqrt((dist*fSum + 0.5*c.massKin*vv) / s)
}

// VelFromBrakeDist laskee vauhdin, josta etäisyydellä dist voidaan jarruttamalla pysähtyä (tuuli=0).
// VelFromBrakeDist(dist) == MaxEntryVel(dist, 0)
//
// Solve vEntry from equation:
// 0.5 * mass  * vEntry^2 = dist*(fGR + fBrake + 0.5*cDrag*vEntry^2)
func (c *BikeCalc) VelFromBrakeDist(dist float64) float64 {

	fSum := c.fGR + c.fBrake
	if fSum <= 0 {
		return 0
	}
	s := 0.5 * (c.massKin - dist*c.cDrag)
	if s <= 0 {
		return maxVEL
	}
	return math.Sqrt(dist * fSum / s)
}
