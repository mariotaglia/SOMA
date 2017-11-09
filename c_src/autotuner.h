/* Copyright (C) 2016-2017 Ludwig Schneider
   Copyright (C) 2016 Ulrich Welling
   Copyright (C) 2016-2017 Marcel Langenberg
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

//! \file autotuner.h Autotuner related structures and functions.
#ifndef SOMA_AUTO_TUNER_H
#define SOMA_AUTO_TUNER_H

//! \file autotuner.h
//! \brief Code for the Autotuner.

#include <stdbool.h>

//! Number of elements that are tried for the autotuner.
//! \warning There is a problem with one kernel, if you increase this.
//! \todo Ask PGI about this.
#define AUTO_TUNER_N_ELEMENTS  12U

//! \brief Autotuner to struct to hold data, which is necessary to determine
//! an optimal parameter for ACC kernels.
typedef struct
    {
        //! Array holding all tuning parameter that are tried. \private
        unsigned int trials[AUTO_TUNER_N_ELEMENTS];
        //! Summed timings for each of the values to try. \private
        double trial_times[AUTO_TUNER_N_ELEMENTS];
        //! Temporary variable store the last start of timing. \private
        double last_start;
        //! Variable to ensure proper user usage. \private
        bool started;
        //! Store the actual iteration that is performed. \private
        unsigned int iteration;
        //! Index which trial is tried out next.\private
        unsigned int next_trial;
        //! Indication, whethter the optimal value has been found.
        bool equilibrated;
        //! Value to use for next iteration. It changes during equilibration, but is constant afterwards.
        unsigned int value;
    }Autotuner;

//! Init of an Autotuner.
//!
//! (C++ constructor)
//! \param a Autotuner to init.
//! \return Errorcode.
int init_autotuner(Autotuner*a);

//! Start the timing of an Autotuner before the function to time.
//!
//! Prepares the next value for iteration.
//! Returns if Autotuner is equilibrated.
//! \param a Autotuner to start.
//! \return Errorcode
int start_autotuner(Autotuner*a);

//! End the timing of an Autotuner after the function to time.
//!
//! Finishes the timing of a function and does some book keeping.
//! If equilibration is done, evaluation is performed.
//! Returns if equilibrated.
//! \param a Autotuner to end.
//! \return Errorcode
int end_autotuner(Autotuner*a);

//! Evaluate a fully equilibrated Autotuner and set the parameter to final configuration.
//!
//! \param a Autotuner to evaluate.
//! \return Errorcode
int evaluate_autotuner(Autotuner*a);

//! Restart an Autotuner.
//!
//! Trigger reevaluation of optimal parameter.
//! \param a Autotuner to restart.
//! \return Errorcode.
int restart_autotuner(Autotuner*a);

#endif//SOMA_AUTO_TUNER_H
