/**
 * @file seisdata_load.h
 * @author Benjamin Yang @ National Taiwan University (b98204032@gmail.com)
 * @brief
 * @version 1.0.0
 * @date 2024-04-07
 *
 * @copyright Copyright (c) 2024
 *
 */

#pragma once

#include <postmajor.h>
/* */
int seisdata_load_sac( SNL_INFO *, const char * );
int seisdata_load_ms( SNL_INFO *, const char * );
int seisdata_load_tank( SNL_INFO *, const char * );
