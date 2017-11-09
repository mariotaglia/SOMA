/* Copyright (C) 2016 Ludwig Schneider
   Copyright (C) 2016 Ulrich Welling
   Copyright (C) 2016 Marcel Langenberg
   Copyright (C) 2016 Fabien Leonforte
   Copyright (C) 2016 Juan Orozco
   Copyright (C) 2016 Yongzhi Ren

 This file is part of SOMA.

 SOMA is free software: you can redistribute it and/or modify
 it under the terms of the GNU Lesser General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 SOMA is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU Lesser General Public License for more details.

 You should have received a copy of the GNU Lesser General Public License
 along with SOMA.  If not, see <http://www.gnu.org/licenses/>.
*/

//! \file autotuner.c
//! \brief Implementation of autotuner.h

#include "autotuner.h"
#include <assert.h>
#include <stdio.h>
#include <mpi.h>
#include <stdbool.h>
#include <time.h>

//! Step size between different elements.
static const unsigned int AUTO_TUNER_STEP = 16;
//! Number of iterations an autotuner runs, before an optimal value is determined.
static const unsigned int AUTO_TUNER_ITERATIONS = 2000;

int init_autotuner(Autotuner*a)
    {
    assert(a);
    for(unsigned int i=0; i < AUTO_TUNER_N_ELEMENTS; i++)
	{
	a->trials[i] = (i+1)*AUTO_TUNER_STEP;
	a->trial_times[i] = 0.;
	}
    a->last_start = 0.;
    a->started = false;
    a->iteration = 0;
    a->next_trial = 0;
    a->equilibrated = false;
    //No need to set a value.
    return 0;
    }

int start_autotuner(Autotuner*a)
    {
    assert(a);
    if(a->equilibrated)
	return 0;
    if( a->started )
	{
	fprintf(stderr,"ERROR: %s:%d Already started Autotuner start.\n",__FILE__,__LINE__);
	return -1;
	}
    assert(a->next_trial < AUTO_TUNER_N_ELEMENTS);
    a->value = a->trials[ a->next_trial ];

    a->last_start = MPI_Wtime();
    a->started = true;
    return 0;
    }

int end_autotuner(Autotuner*a)
    {
    assert(a);
    if(a->equilibrated)
	return 0;
    const double end = MPI_Wtime();
    if( ! a->started)
	{
	fprintf(stderr,"ERROR: %s:%d No started Autotuner end.\n",__FILE__,__LINE__);
	return -1;
	}
    a->started = false;

    assert(a->next_trial < AUTO_TUNER_N_ELEMENTS);
    a->trial_times[ a->next_trial ] += end - a->last_start;

    a->next_trial = (a->next_trial+1) % AUTO_TUNER_N_ELEMENTS;
    if(a->next_trial == 0) a->iteration++;

    if( a->iteration > AUTO_TUNER_ITERATIONS)
	return evaluate_autotuner(a);
    return 0;
    }

int evaluate_autotuner(Autotuner*a)
    {
    assert( a );
    if( a->iteration <= AUTO_TUNER_ITERATIONS)
	{
	fprintf(stderr,"ERROR: %s:%d Evaluation of a not finished Autotuner.\n",__FILE__,__LINE__);
	return -1;
	}
    if( a->equilibrated)
	{
	fprintf(stderr,"ERROR: %s:%d Evaluation of an already finished Autotuner.\n",__FILE__,__LINE__);
	return -1;
	}

    assert( AUTO_TUNER_N_ELEMENTS > 0);
    time_t min = a->trial_times[0];
    unsigned int argmin = 0;
    for(unsigned int i=1 ; i < AUTO_TUNER_N_ELEMENTS; i++)
	if( a->trial_times[i] < min)
	    {
	    min = a->trial_times[i];
	    argmin = i;
	    }
    a->next_trial = argmin;
    a->value = a->trials[argmin];

    a->equilibrated = true;
    return 0;
    }

int restart_autotuner(Autotuner*a)
    {
    return init_autotuner(a);
    }
