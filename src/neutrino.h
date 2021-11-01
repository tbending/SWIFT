/*******************************************************************************
 * This file is part of SWIFT.
 * Coypright (c) 2021 Willem Elbers (whe@willemelbers.com)
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 ******************************************************************************/
#ifndef SWIFT_NEUTRINO_H
#define SWIFT_NEUTRINO_H

/* Config parameters. */
#include "../config.h"

/* Select the correct neutrino model */
#if defined(NEUTRINO_NONE)
#include "./neutrino/none/neutrino.h"
#elif defined(NEUTRINO_DEFAULT)
#include "./neutrino/Default/neutrino.h"
#else
#error "Invalid choice of neutrino model"
#endif

#endif
