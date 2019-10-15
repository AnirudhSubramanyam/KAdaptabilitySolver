/***************************************************************************************/
/*                                                                                     */
/*  Copyright 2018 by Anirudh Subramanyam, Chrysanthos Gounaris and Wolfram Wiesemann  */
/*                                                                                     */
/*  Licensed under the FreeBSD License (the "License").                                */
/*  You may not use this file except in compliance with the License.                   */
/*  You may obtain a copy of the License at                                            */
/*                                                                                     */
/*  https://www.freebsd.org/copyright/freebsd-license.html                             */
/*                                                                                     */
/***************************************************************************************/

#ifndef CONSTANTS_H_
#define CONSTANTS_H_

// Exception numbers
#define EXCEPTION_CPXINIT 1    // CPLEX could not be initialized
#define EXCEPTION_CPXEXIT 2    // CPLEX could not be closed
#define EXCEPTION_CPXNEWCOLS 3 // Error while trying to create new columns
#define EXCEPTION_CPXNEWROWS 4 // Error while trying to create new rows
#define EXCEPTION_X 5          // Invalid size of solution vector
#define EXCEPTION_K 6          // Invalid policy index
#define EXCEPTION_Q 7          // Invalid size of uncertain parameters

// Defaults
#define TIME_LIMIT 7200
#define MEMORY_LIMIT 3072
#define NUM_THREADS 1

#endif
