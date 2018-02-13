/*
  Copyright (c) 2016 Statoil ASA
  Copyright (c) 2013-2015 Andreas Lauser
  Copyright (c) 2013 SINTEF ICT, Applied Mathematics.
  Copyright (c) 2013 Uni Research AS
  Copyright (c) 2015 IRIS AS

  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.
*/
#ifndef RESTART_IO_HPP
#define RESTART_IO_HPP

#include <vector>
#include <map>
#include <cstddef>
#include <time.h>
#include <ctime>
#include <cstring>
#include <type_traits>
#include <sys/types.h> 
#include <sys/stat.h>

//#ifndef ERT_VECTOR_H
//#define ERT_VECTOR_H

//#include <ert/util/node_data.h>
//#include <ert/util/type_macros.h>
//#include <ert/util/int_vector.h>

//jal - comment on the following include to avoid duplication
#include <opm/parser/eclipse/Units/UnitSystem.hpp>
#include <opm/parser/eclipse/EclipseState/Runspec.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/Well.hpp>

#include <opm/output/data/Cells.hpp>
#include <opm/output/data/Solution.hpp>
#include <opm/output/data/Wells.hpp>
#include <opm/output/eclipse/libECLRestart.hpp>
#include <opm/output/eclipse/RestartValue.hpp>

//#include <ert/ecl/EclKW.hpp>
#include <ert/ecl/FortIO.hpp>
//#include <ert/ecl/ecl_rsthead.h>
//#include <ert/ecl/ecl_rst_file.h>
//#include <ert/util/util.h>
#include <ert/ecl/fortio.h>

namespace Opm {

class EclipseGrid;
class EclipseState;
class Phases;
class Schedule;

namespace RestartIO {

    
void save(const std::string& filename,
          int report_step,
          double seconds_elapsed,
          data::Solution cells,
          data::Wells wells,
          const EclipseState& es,
          const EclipseGrid& grid,
	  const Schedule& schedule,
          std::map<std::string, std::vector<double>> extra_data,
	  bool write_double);
// orig version jal  	  bool write_double = false);
// orig version jal          std::map<std::string, std::vector<double>> extra_data = {},

RestartValue load( const std::string& filename,
                   int report_step,
                   const std::map<std::string, RestartKey>& keys,
                   const EclipseState& es,
                   const EclipseGrid& grid,
		           const Schedule& schedule,
                   const std::map<std::string, bool>& extra_keys);
// orig version jal            const std::map<std::string, bool>& extra_keys = {});
 
}
}
#endif // RESTART_IO_HPP