//------------------------------------------------------------------------------
// AUSM+-up flux function (Liou, 2006)
//
// All-speed AUSM+ scheme with pressure/velocity diffusion terms for low
// Mach number flows.  Function signature matches KurganovFlux / TadmorFlux
// so it can be registered in the same fluxMap.
//------------------------------------------------------------------------------
void AUSMFlux
(
    const bool& downgradeSingleT,
    const bool& downgradeSingleTv,
    const bool& meshMoving,
    const scalar& mesh_phi,
    const vector& Sf,
    const scalar& cSf_pos,
    const scalar& cSf_neg,
    const scalar& rho_pos,
    const scalar& rho_neg,
    const PtrList<scalar>& pD_pos,
    const PtrList<scalar>& pD_neg,
    const vector& U_pos,
    const vector& U_neg,
    const vector& rhoU_pos,
    const vector& rhoU_neg,
    const scalar& p_pos,
    const scalar& p_neg,
    const PtrList<scalar>& evk_pos,
    const PtrList<scalar>& evk_neg,
    const scalar& e_pos,
    const scalar& e_neg,
    scalar& amaxSf,
    vector& aByU,
    scalar& rhoFlux,
    vector& rhoUFlux,
    PtrList<scalar>& rhoEvkFlux,
    scalar& rhoEFlux
)
{
    // --- Mach number splitting functions ---
    // M+(4) : Eq.(22) in Liou 2006
    auto Mplus4 = [](const scalar M) -> scalar
    {
        if (mag(M) >= 1.0)
        {
            return 0.5*(M + mag(M));
        }
        else
        {
            return 0.25*sqr(M + 1.0) + 0.125*sqr(sqr(M) - 1.0);
        }
    };

    // M-(4) : Eq.(22) in Liou 2006
    auto Mminus4 = [](const scalar M) -> scalar
    {
        if (mag(M) >= 1.0)
        {
            return 0.5*(M - mag(M));
        }
        else
        {
            return -0.25*sqr(M - 1.0) - 0.125*sqr(sqr(M) - 1.0);
        }
    };

    // P+(5) : Eq.(24) with alpha = 3/16
    auto Pplus5 = [](const scalar M, const scalar alpha = 3.0/16.0) -> scalar
    {
        if (mag(M) >= 1.0)
        {
            return 0.5*(1.0 + sign(M));
        }
        else
        {
            return 0.25*sqr(M + 1.0)*(2.0 - M)
                 + alpha*M*sqr(sqr(M) - 1.0);
        }
    };

    // P-(5) : Eq.(24) with alpha = 3/16
    auto Pminus5 = [](const scalar M, const scalar alpha = 3.0/16.0) -> scalar
    {
        if (mag(M) >= 1.0)
        {
            return 0.5*(1.0 - sign(M));
        }
        else
        {
            return 0.25*sqr(M - 1.0)*(2.0 + M)
                 - alpha*M*sqr(sqr(M) - 1.0);
        }
    };

    // --- Face normal unit vector and face area ---
    const scalar magSf = mag(Sf);
    const vector n = Sf / magSf;

    // --- Face-projected velocities ---
    scalar un_pos = U_pos & n;
    scalar un_neg = U_neg & n;

    // Make velocities relative to mesh motion
    scalar mesh_un = 0.0;
    if (meshMoving)
    {
        mesh_un = mesh_phi / magSf; // mesh flux → velocity normal to face
        un_pos -= mesh_un;
        un_neg -= mesh_un;
    }

    // --- Interface speed of sound (geometric mean for robustness) ---
    // c_pos = cSf_pos / |Sf|, c_neg = cSf_neg / |Sf|
    const scalar c_pos = cSf_pos / magSf;
    const scalar c_neg = cSf_neg / magSf;
    const scalar c_half = sqrt(c_pos * c_neg);

    // --- Left/Right Mach numbers ---
    const scalar M_L = un_pos / c_half;
    const scalar M_R = un_neg / c_half;

    // --- Basic AUSM+ splits ---
    const scalar Mp = Mplus4(M_L);
    const scalar Mm = Mminus4(M_R);

    // --- All-speed pressure diffusion term (Kp) ---
    // Kp = 0.25, sigma = 1.0
    const scalar Kp = 0.25;
    const scalar M_bar_sq = 0.5*(sqr(un_pos) + sqr(un_neg)) / sqr(c_half);
    const scalar f_a = max(1.0 - 1.0*M_bar_sq, 0.0); // sigma = 1.0

    // --- Interface Mach number with pressure diffusion ---
    const scalar rho_bar = 0.5*(rho_pos + rho_neg);
    const scalar Mp_corr = -Kp * f_a * (p_neg - p_pos)
                         / (rho_bar * sqr(c_half) + SMALL);

    const scalar M_half = Mp + Mm + Mp_corr;

    // --- Interface pressure with velocity diffusion ---
    // Ku = 0.75
    const scalar Ku = 0.75;
    const scalar Pu_corr = -Ku * Pplus5(M_L)*Pminus5(M_R)
                          * (rho_pos + rho_neg)
                          * c_half * (un_neg - un_pos);

    const scalar P_half = Pplus5(M_L)*p_pos + Pminus5(M_R)*p_neg + Pu_corr;

    // --- Mass flux ---
    // ṁ = c_half * M_half * rho_upwind  [per unit area]
    scalar mDot;       // mass flux per unit area
    scalar rho_upwind; // upwind density
    vector U_upwind;   // upwind velocity
    scalar e_upwind;   // upwind total specific energy
    scalar p_upwind;   // upwind pressure (for enthalpy)

    if (M_half > 0.0)
    {
        // Flow from left (pos) to right (neg)
        rho_upwind = rho_pos;
        U_upwind   = U_pos;
        e_upwind   = e_pos;
        p_upwind   = p_pos;
        mDot       = c_half * M_half * rho_pos;
    }
    else
    {
        // Flow from right (neg) to left (pos)
        rho_upwind = rho_neg;
        U_upwind   = U_neg;
        e_upwind   = e_neg;
        p_upwind   = p_neg;
        mDot       = c_half * M_half * rho_neg;
    }

    // Total enthalpy H = e + p/rho
    const scalar H_upwind = e_upwind + p_upwind/rho_upwind;

    // --- Face fluxes ---
    // Relative velocity at face (for kinetic energy in mesh-moving case)
    const scalar un_upwind = U_upwind & n - mesh_un;

    // Density flux: ṁ * |Sf|
    rhoFlux = mDot * magSf;

    // Momentum flux: ṁ * U_upwind * |Sf| + P_half * Sf
    rhoUFlux = mDot * U_upwind * magSf + P_half * Sf;

    // Energy flux: ṁ * H_upwind * |Sf|
    // In mesh-moving case, add mesh_phi * P_half for pressure work
    rhoEFlux = mDot * H_upwind * magSf;
    if (meshMoving)
    {
        rhoEFlux += mesh_phi * P_half;
    }

    // Vibrational energy fluxes
    if (downgradeSingleT)
    {
        rhoEvkFlux.clear();
    }
    else if (downgradeSingleTv)
    {
        // Single vib. mode: ṁ * ev_upwind * |Sf|
        if (M_half > 0.0)
        {
            rhoEvkFlux[0] = mDot * evk_pos[0] * magSf;
        }
        else
        {
            rhoEvkFlux[0] = mDot * evk_neg[0] * magSf;
        }
    }
    else
    {
        // Per-species vib. energy: ṁ * Y_s * ev_s * |Sf|
        // mDot = c*M*rho_upwind, Y_s = pD_s/rho, so the rho factors cancel:
        // mDot * Y * evk * |Sf| = c*M * pD * evk * |Sf|  [J/s]
        forAll(rhoEvkFlux, speciei)
        {
            if (M_half > 0.0)
            {
                rhoEvkFlux[speciei] = mDot * (pD_pos[speciei]/rho_pos) * evk_pos[speciei] * magSf;
            }
            else
            {
                rhoEvkFlux[speciei] = mDot * (pD_neg[speciei]/rho_neg) * evk_neg[speciei] * magSf;
            }
        }
    }

    // --- amaxSf for CFL: (c + max(|un_L|, |un_R|)) * |Sf| ---
    amaxSf = (c_half + max(mag(un_pos), mag(un_neg))) * magSf;

    // --- aByU for viscous stress (σ·U term) ---
    aByU = U_upwind;
}
