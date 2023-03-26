//-----------------------------------------------------------------------------------
// d-SEAMS - Deferred Structural Elucidation Analysis for Molecular Simulations
//
// Copyright (c) 2018--present d-SEAMS core team
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the MIT License as published by
// the Open Source Initiative.
//
// A copy of the MIT License is included in the LICENSE file of this repository.
// You should have received a copy of the MIT License along with this program.
// If not, see <https://opensource.org/licenses/MIT>.
//-----------------------------------------------------------------------------------

#ifndef _OPT_PARSER_H_
#define _OPT_PARSER_H_

// This is bundled
#include <cxxopts.hpp>

cxxopts::ParseResult parse(int argc, char *argv[]);

#endif // __OPTIONS_PARSER_H_
