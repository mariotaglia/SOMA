/* Copyright (C) 2016-2018 Ludwig Schneider

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

//! \file signal.c
//! \brief Implementation of signal.h


#include "signal.h"
#include <stdlib.h>
#include <stdio.h>

//! previous sigint handler
void (*prev_sigint_handler)(int) = NULL;
//! pervious sigterm handler
void (*prev_sigterm_handler)(int) = NULL;

//! Store if a SIGINT has been received
volatile sig_atomic_t sigint_recvd = 0;
//! Store if a SIGTERM has been received
volatile sig_atomic_t sigterm_recvd = 0;

//! Handle Function for SIGINT
//! \param sig Signal
void sigint_handler(int sig)
    {
    if(sig != SIGINT)
	return;

    fflush(stdout);
    if(prev_sigint_handler &&
       prev_sigint_handler != SIG_ERR &&
       prev_sigint_handler != SIG_DFL &&
       prev_sigint_handler != SIG_IGN)
	prev_sigint_handler(sig);

    sigint_recvd = 1;
    }

//! Handle Function for SIGTERM
//! \param sig Signal
void sigterm_handler(int sig)
    {
    if(sig != SIGTERM)
	return;

    if(prev_sigterm_handler &&
       prev_sigterm_handler != SIG_ERR &&
       prev_sigterm_handler != SIG_DFL &&
       prev_sigterm_handler != SIG_IGN)
	prev_sigterm_handler(sig);

    sigterm_recvd = 1;
    }

int init_soma_signal(void)
    {
    int ret =0;
    //Set up SIGINT
    void (*retval)(int);
    retval = signal(SIGINT,sigint_handler);
    if(retval == SIG_ERR)
	{
	fprintf(stderr,"ERROR: cannot set signal handler for SIGINT.\n");
	ret += 1;
	}
    if(retval != sigint_handler )
	prev_sigint_handler = retval;
    else
	prev_sigint_handler = NULL;

    //Set up SIGTERM
    retval = signal(SIGTERM,sigterm_handler);
    if(retval == SIG_ERR)
	{
	fprintf(stderr,"ERROR: cannot set signal handler for SIGTERM.\n");
	ret += 2;
	}
    if(retval != sigterm_handler )
	prev_sigterm_handler = retval;
    else
	prev_sigterm_handler = NULL;

    return ret;
    }

int check_signal_stop(void)
    {
    return (sigint_recvd == 1) || (sigterm_recvd == 1);
    }
