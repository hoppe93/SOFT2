/**
 * Methods for initializing momentum-space
 * coordinates in the 'Particle' object.
 */

#include <softlib/Vector.h>
#include "config.h"
#include "PhaseSpace/Particle.h"

using namespace std;

/**
 * Returns the name corresponding to a
 * given coordinate type as a string.
 *
 * t: Coordinate type ID (Particle::COORDINATE_???).
 */
const char* Particle::GetCoordinateName(const int t) {
	switch (t) {
		case Particle::COORDINATE_GAMMA: return "gamma";
		case Particle::COORDINATE_P: return "p";
		case Particle::COORDINATE_PPAR: return "ppar";
		case Particle::COORDINATE_PPERP: return "pperp";
		case Particle::COORDINATE_THETAP: return "thetap";
        case Particle::COORDINATE_ITHETAP: return "ithetap";
		case Particle::COORDINATE_XI: return "xi";
		default: return "n/a";
	}
}

/**
 * Initialize the momentum properties of this particle
 * from arbitrary coordinates. This function initializes
 * the cartesian momentum, the magnitude of the momentum
 * vector, the parallel and perpendicular momentum, the
 * pitch angle, the cosine of the pitch angle, and the
 * momentum differential elements (dparam1, dparam2).
 * 
 * t1, t2:   Coordinate types (Particle::COORDINATE_???).
 * p1, p2:   Coordinate values.
 * zeta:     Gyro phase.
 * dp1, dp2: Change in p1/p2 direction relative to previous particle.
 */
void Particle::InitializeMomentum(const int t1, const int t2, const slibreal_t p1, const slibreal_t p2, const slibreal_t zeta, const slibreal_t dp1, const slibreal_t dp2, const slibreal_t dzeta) {
	this->momentum1 = t1;
	this->momentum2 = t2;
    this->zeta = zeta;
    this->SetDZeta(dzeta);

	switch (t1) {
		case Particle::COORDINATE_GAMMA:
			if (t2 == Particle::COORDINATE_PPAR) Init_GammaPPar(p1, p2, &(this->dparam1), &(this->dparam2));
			else if (t2 == Particle::COORDINATE_THETAP) Init_GammaThetap(p1, p2, &(this->dparam1), &(this->dparam2));
            else if (t2 == Particle::COORDINATE_ITHETAP) Init_GammaThetap(p1, M_PI-p2, &(this->dparam1), &(this->dparam2));
			else if (t2 == Particle::COORDINATE_XI) Init_GammaXi(p1, p2, &(this->dparam1), &(this->dparam2));
			else
				throw ParticleException("Invalid combination of coordinates! Coordinate 'gamma' can be combined with 'ppar', 'thetap' and 'xi'.");
			break;
		case Particle::COORDINATE_P:
			if (t2 == Particle::COORDINATE_PPAR) Init_PPPar(p1, p2, &(this->dparam1), &(this->dparam2));
			else if (t2 == Particle::COORDINATE_THETAP) Init_PThetap(p1, p2, &(this->dparam1), &(this->dparam2));
			else if (t2 == Particle::COORDINATE_ITHETAP) Init_PThetap(p1, M_PI-p2, &(this->dparam1), &(this->dparam2));
			else if (t2 == Particle::COORDINATE_XI) Init_PXi(p1, p2, &(this->dparam1), &(this->dparam2));
			else
				throw ParticleException("Invalid combination of coordinates! Coordinate 'p' can be combined with 'ppar', 'thetap' and 'xi'.");
			break;
		case Particle::COORDINATE_PPAR:
			if (t2 == Particle::COORDINATE_GAMMA) Init_GammaPPar(p2, p1, &(this->dparam2), &(this->dparam1));
			else if (t2 == Particle::COORDINATE_P) Init_PPPar(p2, p1, &(this->dparam2), &(this->dparam1));
			else if (t2 == Particle::COORDINATE_PPERP) Init_PParPPerp(p1, p2, &(this->dparam1), &(this->dparam2));
			else if (t2 == Particle::COORDINATE_THETAP) Init_PParThetap(p1, p2, &(this->dparam1), &(this->dparam2));
			else if (t2 == Particle::COORDINATE_ITHETAP) Init_PParThetap(p1, M_PI-p2, &(this->dparam1), &(this->dparam2));
			else if (t2 == Particle::COORDINATE_XI) Init_PParXi(p1, p2, &(this->dparam1), &(this->dparam2));
			else
				throw ParticleException("Invalid combination of coordinates! Coordinate 'ppar' can be combined with 'gamma', 'p', 'pperp', 'thetap' and 'xi'.");
			break;
		case Particle::COORDINATE_PPERP:
			if (t2 == Particle::COORDINATE_PPAR) Init_PParPPerp(p2, p1, &(this->dparam2), &(this->dparam1));
			else if (t2 == Particle::COORDINATE_THETAP) Init_PPerpThetap(p1, p2, &(this->dparam1), &(this->dparam2));
			else if (t2 == Particle::COORDINATE_ITHETAP) Init_PPerpThetap(p1, M_PI-p2, &(this->dparam1), &(this->dparam2));
			else if (t2 == Particle::COORDINATE_XI) Init_PPerpXi(p1, p2, &(this->dparam1), &(this->dparam2));
			else
				throw ParticleException("Invalid combination of coordinates! Coordinate 'pperp' can be combined with 'ppar', 'thetap' and 'xi'.");
			break;
		case Particle::COORDINATE_THETAP:
			if (t2 == Particle::COORDINATE_GAMMA) Init_GammaThetap(p2, p1, &(this->dparam2), &(this->dparam1));
			else if (t2 == Particle::COORDINATE_P) Init_PThetap(p2, p1, &(this->dparam2), &(this->dparam1));
			else if (t2 == Particle::COORDINATE_PPAR) Init_PParThetap(p2, p1, &(this->dparam2), &(this->dparam1));
			else if (t2 == Particle::COORDINATE_PPERP) Init_PPerpThetap(p2, p1, &(this->dparam2), &(this->dparam1));
			else
				throw ParticleException("Invalid combination of coordinates! Coordinate 'thetap' can be combined with 'gamma', 'p', 'ppar' and 'pperp'.");
			break;
		case Particle::COORDINATE_ITHETAP:
			if (t2 == Particle::COORDINATE_GAMMA) Init_GammaThetap(p2, M_PI-p1, &(this->dparam2), &(this->dparam1));
			else if (t2 == Particle::COORDINATE_P) Init_PThetap(p2, M_PI-p1, &(this->dparam2), &(this->dparam1));
			else if (t2 == Particle::COORDINATE_PPAR) Init_PParThetap(p2, M_PI-p1, &(this->dparam2), &(this->dparam1));
			else if (t2 == Particle::COORDINATE_PPERP) Init_PPerpThetap(p2, M_PI-p1, &(this->dparam2), &(this->dparam1));
			else
				throw ParticleException("Invalid combination of coordinates! Coordinate 'ithetap' can be combined with 'gamma', 'p', 'ppar' and 'pperp'.");
			break;
		case Particle::COORDINATE_XI:
			if (t2 == Particle::COORDINATE_GAMMA) Init_GammaXi(p2, p1, &(this->dparam2), &(this->dparam1));
			else if (t2 == Particle::COORDINATE_P) Init_PXi(p2, p1, &(this->dparam2), &(this->dparam1));
			else if (t2 == Particle::COORDINATE_PPAR) Init_PParXi(p2, p1, &(this->dparam2), &(this->dparam1));
			else if (t2 == Particle::COORDINATE_PPERP) Init_PPerpXi(p2, p1, &(this->dparam2), &(this->dparam1));
			else
				throw ParticleException("Invalid combination of coordinates! Coordinate 'xi' can be combined with 'gamma', 'p', 'ppar' and 'pperp'.");
			break;
		default:
			throw ParticleException("Unrecognized type of first momentum coordinate: %d", t1);
	}

    this->dparam1 *= fabs(dp1);
    this->dparam2 *= fabs(dp2);

	this->momentum_initialized = true;
}

/**
 * Convert the given pair of momentum coordinates into a
 * pair of ppar/pperp coordinates.
 * 
 * p1, p2:      Parameter values.
 * t1, t2:      Parameter types.
 * ppar, pperp: On return, contains the corresponding
 *              ppar/pperp values.
 */
void Particle::ToPP(
    const slibreal_t p1, const slibreal_t p2,
    const int t1, const int t2,
    slibreal_t *ppar, slibreal_t *pperp
) {
	switch (t1) {
		case Particle::COORDINATE_GAMMA:
			if (t2 == Particle::COORDINATE_PPAR) ToPP_GammaPPar(p1, p2, ppar, pperp);
			else if (t2 == Particle::COORDINATE_THETAP) ToPP_GammaThetap(p1, p2, ppar, pperp);
			else if (t2 == Particle::COORDINATE_ITHETAP) ToPP_GammaThetap(p1, M_PI-p2, ppar, pperp);
			else if (t2 == Particle::COORDINATE_XI) ToPP_GammaXi(p1, p2, ppar, pperp);
			else
				throw ParticleException("Invalid combination of coordinates! Coordinate 'gamma' can be combined with 'ppar', 'thetap' and 'xi'.");
			break;
		case Particle::COORDINATE_P:
			if (t2 == Particle::COORDINATE_PPAR) ToPP_PPPar(p1, p2, ppar, pperp);
			else if (t2 == Particle::COORDINATE_THETAP) ToPP_PThetap(p1, p2, ppar, pperp);
			else if (t2 == Particle::COORDINATE_ITHETAP) ToPP_PThetap(p1, M_PI-p2, ppar, pperp);
			else if (t2 == Particle::COORDINATE_XI) ToPP_PXi(p1, p2, ppar, pperp);
			else
				throw ParticleException("Invalid combination of coordinates! Coordinate 'p' can be combined with 'ppar', 'thetap' and 'xi'.");
			break;
		case Particle::COORDINATE_PPAR:
			if (t2 == Particle::COORDINATE_GAMMA) ToPP_GammaPPar(p2, p1, ppar, pperp);
			else if (t2 == Particle::COORDINATE_P) ToPP_PPPar(p2, p1, ppar, pperp);
			else if (t2 == Particle::COORDINATE_PPERP) ToPP_PParPPerp(p1, p2, ppar, pperp);
			else if (t2 == Particle::COORDINATE_THETAP) ToPP_PParThetap(p1, p2, ppar, pperp);
			else if (t2 == Particle::COORDINATE_ITHETAP) ToPP_PParThetap(p1, M_PI-p2, ppar, pperp);
			else if (t2 == Particle::COORDINATE_XI) ToPP_PParXi(p1, p2, ppar, pperp);
			else
				throw ParticleException("Invalid combination of coordinates! Coordinate 'ppar' can be combined with 'gamma', 'p', 'pperp', 'thetap' and 'xi'.");
			break;
		case Particle::COORDINATE_PPERP:
			if (t2 == Particle::COORDINATE_PPAR) ToPP_PParPPerp(p2, p1, ppar, pperp);
			else if (t2 == Particle::COORDINATE_THETAP) ToPP_PPerpThetap(p1, p2, ppar, pperp);
			else if (t2 == Particle::COORDINATE_ITHETAP) ToPP_PPerpThetap(p1, M_PI-p2, ppar, pperp);
			else if (t2 == Particle::COORDINATE_XI) ToPP_PPerpXi(p1, p2, ppar, pperp);
			else
				throw ParticleException("Invalid combination of coordinates! Coordinate 'pperp' can be combined with 'ppar', 'thetap' and 'xi'.");
			break;
		case Particle::COORDINATE_THETAP:
			if (t2 == Particle::COORDINATE_GAMMA) ToPP_GammaThetap(p2, p1, ppar, pperp);
			else if (t2 == Particle::COORDINATE_P) ToPP_PThetap(p2, p1, ppar, pperp);
			else if (t2 == Particle::COORDINATE_PPAR) ToPP_PParThetap(p2, p1, ppar, pperp);
			else if (t2 == Particle::COORDINATE_PPERP) ToPP_PPerpThetap(p2, p1, ppar, pperp);
			else
				throw ParticleException("Invalid combination of coordinates! Coordinate 'thetap' can be combined with 'gamma', 'p', 'ppar' and 'pperp'.");
			break;
		case Particle::COORDINATE_ITHETAP:
			if (t2 == Particle::COORDINATE_GAMMA) ToPP_GammaThetap(p2, M_PI-p1, ppar, pperp);
			else if (t2 == Particle::COORDINATE_P) ToPP_PThetap(p2, M_PI-p1, ppar, pperp);
			else if (t2 == Particle::COORDINATE_PPAR) ToPP_PParThetap(p2, M_PI-p1, ppar, pperp);
			else if (t2 == Particle::COORDINATE_PPERP) ToPP_PPerpThetap(p2, M_PI-p1, ppar, pperp);
			else
				throw ParticleException("Invalid combination of coordinates! Coordinate 'thetap' can be combined with 'gamma', 'p', 'ppar' and 'pperp'.");
			break;
		case Particle::COORDINATE_XI:
			if (t2 == Particle::COORDINATE_GAMMA) ToPP_GammaXi(p2, p1, ppar, pperp);
			else if (t2 == Particle::COORDINATE_P) ToPP_PXi(p2, p1, ppar, pperp);
			else if (t2 == Particle::COORDINATE_PPAR) ToPP_PParXi(p2, p1, ppar, pperp);
			else if (t2 == Particle::COORDINATE_PPERP) ToPP_PPerpXi(p2, p1, ppar, pperp);
			else
				throw ParticleException("Invalid combination of coordinates! Coordinate 'xi' can be combined with 'gamma', 'p', 'ppar' and 'pperp'.");
			break;
		default:
			throw ParticleException("Unrecognized type of first momentum coordinate: %d", t1);
	}
}

/**
 * Verify velocity coordinate specifications.
 */
void Particle::VerifyCoordinateSpecification(
    slibreal_t v0, slibreal_t v1, unsigned int n, int type
) {
    switch (type) {
        case Particle::COORDINATE_GAMMA:
            if (v0 <= 1.0 || v1 <= 1.0)
                throw ParticleException("Phasespace coordinate 'gamma' must be greater than 1.");
            break;
        case Particle::COORDINATE_P:
            if (v0 <= 0.0 || v1 <= 0.0)
                throw ParticleException("Phasespace coordinate 'p' must be greater than 0.");
            break;
        case Particle::COORDINATE_PPAR: break;
        case Particle::COORDINATE_PPERP:
            if (v0 <= 0.0 || v1 <= 0.0)
                throw ParticleException("Phasespace coordinate 'pperp' must be greater than 0.");
            break;
        case Particle::COORDINATE_THETAP:
            if ((v0 < 0.0 || v0 > M_PI) ||
                (v1 < 0.0 || v1 > M_PI))
                throw ParticleException("Phasespace coordinate 'thetap' must be given on the interval [0,PI].");
            break;
        case Particle::COORDINATE_ITHETAP:
            if ((v0 < 0.0 || v0 > M_PI) ||
                (v1 < 0.0 || v1 > M_PI))
                throw ParticleException("Phasespace coordinate 'ithetap' must be given on the interval [0,PI].");
            break;
        case Particle::COORDINATE_XI:
            if ((v0 < -1.0 || v0 > 1.0) ||
                (v1 < -1.0 || v1 > 1.0))
                throw ParticleException("Phasespace coordinate 'xi' must be given on the interval [-1,1].");
            break;
        default:
            throw ParticleException("Unrecognized coordinate type specified: %d.", type);
    }

    if (n == 0)
        throw ParticleException("Invalid number of points for coordinate '%s'.", Particle::GetCoordinateName(type));
    else if (n == 1 && v0 != v1)
        throw ParticleException(
            "Number of points in %s is 1, but %s0 != %s1.",
            Particle::GetCoordinateName(type), Particle::GetCoordinateName(type),
            Particle::GetCoordinateName(type)
        );
}

/***********************************
 * METHODS FOR SETTING COORDINATES *
 *                                 *
 * The methods below set the       *
 * following properties of this:   *
 *                                 *
 *   gamma                         *
 *   pmag                          *
 *   ppar                          *
 *   pperp                         *
 *   thetap                        *
 *   xi                            *
 *                                 *
 * Pointers to 'this->dparam1' and *
 * 'this->dparam2' should also be  *
 * passed in the appropriate order *
 * to each of these routines.      *
 ***********************************/
/* \gamma & p_\parallel */
void Particle::ToPP_GammaPPar(
    const slibreal_t gamma, const slibreal_t ppar,
    slibreal_t *pparo, slibreal_t *pperpo
) {
    *pparo = ppar;
    *pperpo = sqrt(gamma*gamma - 1.0 - ppar*ppar);
}
void Particle::Init_GammaPPar(const slibreal_t gamma, const slibreal_t ppar, slibreal_t *dgamma, slibreal_t *dppar) {
	slibreal_t p2 = gamma*gamma-1.0,
			   ppar2 = ppar*ppar;

	this->gamma = gamma;
	this->ppar = ppar;

	this->pmag = sqrt(p2);
	this->xi = ppar / this->pmag;
	this->thetap = acos(this->xi);
	this->pperp = sqrt(p2 - ppar2);

	if (fabs(ppar) > this->pmag)
		throw ParticleException("Invalid value set for coordinate 'ppar'. |ppar| > p.");

	*dgamma = gamma;
	*dppar  = 1.0;
}

/********************************
 * \gamma & \theta_{\mathrm{p}} *
 ********************************/
void Particle::ToPP_GammaThetap(
    const slibreal_t gamma, const slibreal_t thetap,
    slibreal_t *pparo, slibreal_t *pperpo
) {
    slibreal_t p = sqrt(gamma*gamma - 1.0);
    *pparo = p*cos(thetap);
    *pperpo = p*sin(thetap);
}
void Particle::Init_GammaThetap(const slibreal_t gamma, const slibreal_t thetap, slibreal_t *dgamma, slibreal_t *dthetap) {
	slibreal_t sinThetap;

	this->gamma = gamma;
	this->thetap = thetap;

	sinThetap = sin(thetap);
	this->xi = cos(thetap);
	this->pmag = sqrt(gamma*gamma-1.0);
	this->ppar = this->pmag*this->xi;
	this->pperp = this->pmag * sinThetap;
	
	*dgamma  = gamma * this->pmag;
	*dthetap = sinThetap;
}

/****************
 * \gamma & \xi *
 ****************/
void Particle::ToPP_GammaXi(
    const slibreal_t gamma, const slibreal_t xi,
    slibreal_t *pparo, slibreal_t *pperpo
) {
    slibreal_t p = sqrt(gamma*gamma - 1.0);
    *pparo = p * xi;
    *pperpo = p * sqrt(1.0 - xi*xi);
}
void Particle::Init_GammaXi(const slibreal_t gamma, const slibreal_t xi, slibreal_t *dgamma, slibreal_t *dxi) {
	this->gamma = gamma;
	this->xi = xi;

	this->thetap = acos(xi);
	this->pmag = sqrt(gamma*gamma-1.0);
	this->ppar = this->pmag * this->xi;
	this->pperp = this->pmag * sin(this->thetap);

	*dgamma = gamma * this->pmag;
	*dxi = 1.0;
}

/*******************
 * p & p_\parallel *
 *******************/
void Particle::ToPP_PPPar(
    const slibreal_t p, const slibreal_t ppar,
    slibreal_t *pparo, slibreal_t *pperpo
) {
    *pparo = ppar;
    *pperpo = sqrt(p*p - ppar*ppar);
}
void Particle::Init_PPPar(const slibreal_t p, const slibreal_t ppar, slibreal_t *dp, slibreal_t *dppar) {
	slibreal_t p2 = p*p;

	if (fabs(ppar) > p)
		throw ParticleException("Invalid value set for coordinate 'ppar'. |ppar| > p.");

	this->pmag = p;
	this->ppar = ppar;
	this->xi = ppar / p;
	this->thetap = acos(this->xi);
	this->pperp = sqrt(p2 - ppar*ppar);
	this->gamma = sqrt(p2 + 1.0);

	*dp = p;
	*dppar = 1.0;
}

/***************************
 * p & \theta_{\mathrm{p}} *
 ***************************/
void Particle::ToPP_PThetap(
    const slibreal_t p, const slibreal_t thetap,
    slibreal_t *pparo, slibreal_t *pperpo
) {
    *pparo = p*cos(thetap);
    *pperpo = p*sin(thetap);
}
void Particle::Init_PThetap(const slibreal_t p, const slibreal_t thetap, slibreal_t *dp, slibreal_t *dthetap) {
	slibreal_t sinThetap, p2 = p*p;

	this->pmag = p;
	this->thetap = thetap;

	sinThetap = sin(thetap);
	this->xi = cos(thetap);
	this->ppar = p * this->xi;
	this->pperp = p * sinThetap;
	this->gamma = sqrt(p2 + 1.0);

	*dp = p2;
	*dthetap = sinThetap;
}

/***********
 * p & \xi *
 ***********/
void Particle::ToPP_PXi(
    const slibreal_t p, const slibreal_t xi,
    slibreal_t *pparo, slibreal_t *pperpo
) {
    *pparo  = p * xi;
    *pperpo = p * sqrt(1.0 - xi*xi);
}
void Particle::Init_PXi(const slibreal_t p, const slibreal_t xi, slibreal_t *dp, slibreal_t *dxi) {
	slibreal_t p2 = p*p;

	this->pmag = p;
	this->xi = xi;

	this->thetap = acos(xi);
	this->ppar = p * xi;
	this->pperp = p * sqrt(1.0 - xi*xi);
	this->gamma = sqrt(p2 + 1.0);

	*dp = p2;
	*dxi = 1.0;
}

/*************************
 * p_\parallel & p_\perp *
 *************************/
void Particle::ToPP_PParPPerp(
    const slibreal_t ppar, const slibreal_t pperp,
    slibreal_t *pparo, slibreal_t *pperpo
) {
    *pparo = ppar;
    *pperpo = pperp;
}
void Particle::Init_PParPPerp(const slibreal_t ppar, const slibreal_t pperp, slibreal_t *dppar, slibreal_t *dpperp) {
	slibreal_t p2 = ppar*ppar + pperp*pperp;

	this->ppar = ppar;
	this->pperp = pperp;

	this->pmag = sqrt(p2);
	this->gamma = sqrt(p2 + 1.0);
	this->xi = ppar / this->pmag;
	this->thetap = acos(this->xi);

	*dppar = 1.0;
	*dpperp = pperp;
}

/*************************************
 * p_\parallel & \theta_{\mathrm{p}} *
 *************************************/
void Particle::ToPP_PParThetap(
    const slibreal_t ppar, const slibreal_t thetap,
    slibreal_t *pparo, slibreal_t *pperpo
) {
    *pparo = ppar;
    *pperpo = ppar * tan(thetap);
}
void Particle::Init_PParThetap(const slibreal_t ppar, const slibreal_t thetap, slibreal_t *dppar, slibreal_t *dthetap) {
	slibreal_t sinThetap = sin(thetap);

	this->ppar = ppar;
	this->thetap = thetap;

	this->xi = cos(thetap);
	this->pmag = ppar / this->xi;
	this->pperp = this->pmag * sinThetap;
	this->gamma = sqrt(this->pmag*this->pmag + 1.0);

	if (signbit(ppar) != signbit(this->xi))
		throw ParticleException("Invalid value set for coordinate pair ppar / thetap. sign(ppar) != sign(xi).");

	*dppar = ppar*ppar;
	*dthetap = sinThetap / fabs(this->xi*this->xi*this->xi);
}

/*********************
 * p_\parallel & \xi *
 *********************/
void Particle::ToPP_PParXi(
    const slibreal_t ppar, const slibreal_t xi,
    slibreal_t *pparo, slibreal_t *pperpo
) {
    *pparo = ppar;
    *pperpo = ppar / xi * sqrt(1.0 - xi*xi);
}
void Particle::Init_PParXi(const slibreal_t ppar, const slibreal_t xi, slibreal_t *dppar, slibreal_t *dxi) {
	slibreal_t xi2 = xi*xi;

	if (signbit(ppar) != signbit(xi))
		throw ParticleException("Invalid value set for coordinate pair ppar / xi. sign(ppar) != sign(xi).");

	this->ppar = ppar;
	this->xi = xi;

	this->thetap = acos(xi);
	this->pmag = ppar / xi;
	this->pperp = this->pmag * sqrt(1.0 - xi2);
	this->gamma = sqrt(this->pmag*this->pmag + 1.0);

	*dppar = ppar*ppar;
	*dxi   = 1.0 / fabs(xi2*xi);
}

/*********************************
 * p_\perp & \theta_{\mathrm{p}} *
 *********************************/
void Particle::ToPP_PPerpThetap(
    const slibreal_t pperp, const slibreal_t thetap,
    slibreal_t *pparo, slibreal_t *pperpo
) {
    *pperpo = pperp;
    *pparo = pperp / tan(thetap);
}
void Particle::Init_PPerpThetap(const slibreal_t pperp, const slibreal_t thetap, slibreal_t *dpperp, slibreal_t *dthetap) {
	slibreal_t sinThetap = sin(thetap);
	this->pperp = pperp;
	this->thetap = thetap;

	this->xi = cos(thetap);
	this->pmag = pperp / sinThetap;
	this->ppar = this->pmag * this->xi;
	this->gamma = sqrt(this->pmag*this->pmag + 1.0);

	*dpperp = pperp*pperp;
	*dthetap = 1.0 / (sinThetap*sinThetap);
}

/*****************
 * p_\perp & \xi *
 *****************/
void Particle::ToPP_PPerpXi(
    const slibreal_t pperp, const slibreal_t xi,
    slibreal_t *pparo, slibreal_t *pperpo
) {
    *pparo = pperp * xi / sqrt(1.0 - xi*xi);
    *pperpo = pperp;
}
void Particle::Init_PPerpXi(const slibreal_t pperp, const slibreal_t xi, slibreal_t *dpperp, slibreal_t *dxi) {
	slibreal_t sinThetap2 = 1.0-xi*xi,
			   sinThetap  = sqrt(sinThetap2);
	
	this->pperp = pperp;
	this->xi = xi;

	this->pmag = pperp / sinThetap;
	this->ppar = this->pmag * xi;
	this->thetap = acos(xi);
	this->gamma = sqrt(this->pmag*this->pmag + 1.0);

	*dpperp = pperp*pperp;
	*dxi    = 1.0 / fabs(sinThetap2 * sinThetap);
}

