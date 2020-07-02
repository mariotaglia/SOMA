
#include "serv_util.h"


bool need_to_do(unsigned int obs_delta, unsigned int time){
    return (obs_delta != 0 && time % obs_delta != 0);
}

bool no_observables(const Ana_Info * ai, unsigned int t){


    if (need_to_do(ai->delta_mc_Re, t))
        return false;
    if (need_to_do(ai->delta_mc_Rg, t))
        return false;
    if (need_to_do(ai->delta_mc_b_anisotropy, t))
        return false;
    if (need_to_do(ai->delta_mc_density_field, t))
        return false;
    if (need_to_do(ai->delta_mc_acc_ratio, t))
        return false;
    if (need_to_do(ai->delta_mc_MSD, t))
        return false;
    if (need_to_do(ai->delta_mc_dump, t))
        return false;
    if (need_to_do(ai->delta_mc_density_var, t))
        return false;
    if (need_to_do(ai->delta_mc_non_bonded_energy, t))
        return false;
    if (need_to_do(ai->delta_mc_bonded_energy, t))
        return false;
    if (need_to_do(ai->delta_mc_umbrella_field, t))
        return false;
    if (need_to_do(ai->delta_mc_dynamical_structure, t))
        return false;
    if (need_to_do(ai->delta_mc_static_structure, t))
        return false;

    return true;

}

bool has_poly_obs(const Ana_Info * ai, unsigned int t)
{

    if (need_to_do(ai->delta_mc_MSD, t))
        return true;
    if (need_to_do(ai->delta_mc_Re, t))
        return true;
    if (need_to_do(ai->delta_mc_b_anisotropy, t))
        return true;
    if (need_to_do(ai->delta_mc_Rg, t))
        return true;
    if (need_to_do(ai->delta_mc_dynamical_structure, t))
        return true;
    if (need_to_do(ai->delta_mc_static_structure, t))
        return true;
    if (need_to_do(ai->delta_mc_bonded_energy, t))
        return true;
    return false;
}

bool has_field_obs(const Ana_Info * ai, unsigned int t)
{
    if (need_to_do(ai->delta_mc_density_var, t))
        return true;
    if (need_to_do(ai->delta_mc_non_bonded_energy, t))
        return true;

    return false;
}

bool has_omega_field_obs(const Ana_Info * ai, unsigned int t)
{
    return need_to_do(ai->delta_mc_non_bonded_energy, t);
}