#include "soma_config.h"
#include "phase.h"

#ifndef ELECRO_HELPER_H
#define ELECTRO_HELPER_H

struct Phase;

int call_PB(struct Phase *const p);
int call_J(struct Phase *const p);
//! Calls kinsol

int call_EN(const struct Phase *const p);
//! Calls electroneutrality efield solver

extern void calc_ions(struct Phase *const p);
extern void calc_born_S(struct Phase *const p);
void calc_invbls(struct Phase *const p);

void update_invblav(const struct Phase *const p);
void update_d_invblav(const struct Phase *const p);
void update_exp_born(const struct Phase *const p);
void update_rhoF(const struct Phase *const p);
void update_NB(const struct Phase *const p);

void update_electric_field(const struct Phase *const p);
void calc_J_umbrella(const struct Phase *const p);


#endif
