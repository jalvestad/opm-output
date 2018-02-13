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
#include <opm/output/eclipse/RestartIO.hpp>
#include <string>
#include <vector>
#include <cstddef>
#include <time.h>
#include <ctime>
#include <cstring>
#include <stdlib.h> 
#include <iostream>
#include <stdio.h>
#include <sys/stat.h>
#include <fstream> 
#include <fcntl.h>
#include <cstdint>
#include <cstdlib>
#include <cstdarg>
#include <cstdio>
#include <stdarg.h>
#include <assert.h>
#include <signal.h>
#include <pthread.h>
//#include <ert/util/@TYPE@_vector.h>



//  This function is purely a helper function for util_abort().
// start define for util_abort
//#define __USE_GNU       // Must be defined to get access to the dladdr() function; Man page says the symbol should be: _GNU_SOURCE but that does not seem to work?
//#define _GNU_SOURCE     // Must be defined _before_ #include <errno.h> to get the symbol 'program_invocation_name'.

//#include <errno.h>
//#include <stdlib.h>
//#include <string.h>

//#include <ert/util/util.h>
//#include <ert/util/test_util.h>

//#include <stdbool.h>

//#include <dlfcn.h>
//#include <execinfo.h>
//#include <pthread.h>
//#include <pwd.h>
//#include <signal.h>
//#include <unistd.h>

//#if !defined(__GLIBC__)         /* note: not same as __GNUC__ */
//#  if defined (__APPLE__)
//#    include <stdlib.h>         /* alloca   */
//#    include <sys/syslimits.h>  /* PATH_MAX */
//#    include <mach-o/dyld.h>    /* _NSGetExecutablePath */
//#  elif defined (__LINUX__)
//#    include <stdlib.h>         /* alloca   */
//#    include <limits.h>         /* PATH_MAX */
//#    include <unistd.h>         /* readlink */
//#  else
//#    error No known program_invocation_name in runtime library
//#  endif
//#endif  

#include <opm/parser/eclipse/EclipseState/Grid/EclipseGrid.hpp>
#include <opm/parser/eclipse/EclipseState/EclipseState.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/Schedule.hpp>

#include <opm/output/eclipse/libECLRestart.hpp>

//#include <ert/ecl/EclKW.hpp>
#include <ert/ecl/FortIO.hpp>
#include <ert/ecl/EclFilename.hpp>
//#include <ert/ecl/ecl_kw_magic.h>
//#include <ert/ecl/ecl_init_file.h>
//#include <ert/ecl/ecl_file.h>
//#include <ert/ecl/ecl_kw.h>
//#include <ert/ecl/ecl_type.h>
//#include <ert/ecl/ecl_grid.h>
//#include <ert/ecl/ecl_util.h>
//#include <ert/ecl/ecl_rft_file.h>
//#include <ert/ecl/ecl_file_view.h>
//#include <ert/ecl_well/well_const.h>
//#include <ert/ecl/ecl_rsthead.h>
//#include <ert/util/util.h>
#include <ert/ecl/fortio.h>
#define OPM_XWEL      "OPM_XWEL"
#define OPM_IWEL      "OPM_IWEL"
#define IWEL_KW      "IWEL"
#define ZWEL_KW      "ZWEL"
#define ICON_KW      "ICON"


namespace {
    Opm::UnitSystem::UnitType
    getUnitSystem(const Opm::RestartIO::ecl_kw_type* intehead)
    {
        using UType = Opm::UnitSystem::UnitType;

        //switch (Opm::RestartIO::ecl_kw_iget_int(intehead, INTEHEAD_UNIT_INDEX)) {
	switch (Opm::RestartIO::ecl_kw_iget_type<int>(intehead, Opm::RestartIO::ECL_INT_TYPE, INTEHEAD_UNIT_INDEX)) {
        case 1: return UType::UNIT_TYPE_METRIC;
        case 2: return UType::UNIT_TYPE_FIELD;
        case 3: return UType::UNIT_TYPE_LAB;
        case 4: return UType::UNIT_TYPE_PVT_M;
        }

        return UType::UNIT_TYPE_METRIC; // questionable…
    }
}

namespace Opm {
namespace RestartIO  {

//namespace {

    static const int NIWELZ = 11; //Number of data elements per well in IWEL array in restart file
    static const int NZWELZ = 3;  //Number of 8-character words per well in ZWEL array restart file
    static const int NICONZ = 15; //Number of data elements per completion in ICON array restart file

    /**
     * The constants NIWELZ and NZWELZ referes to the number of
     * elements per well that we write to the IWEL and ZWEL eclipse
     * restart file data arrays. The constant NICONZ refers to the
     * number of elements per completion in the eclipse restart file
     * ICON data array.These numbers are written to the INTEHEAD
     * header.
     *
     * Observe that all of these values are our "current-best-guess"
     * for how many numbers are needed; there might very well be third
     * party applications out there which have a hard expectation for
     * these values.
     */


    inline int to_ert_welltype( const Well& well, size_t timestep ) {
        if( well.isProducer( timestep ) ) return IWEL_PRODUCER;

        switch( well.getInjectionProperties( timestep ).injectorType ) {
        case WellInjector::WATER:
            return IWEL_WATER_INJECTOR;
        case WellInjector::GAS:
            return IWEL_GAS_INJECTOR;
        case WellInjector::OIL:
            return IWEL_OIL_INJECTOR;
        default:
            return IWEL_UNDOCUMENTED_ZERO;
        }
    }

  std::vector<double> double_vector( const ::Opm::RestartIO::ecl_kw_type * ecl_kw ) {
        size_t size = ::Opm::RestartIO::ecl_kw_get_size( ecl_kw );

        if (::Opm::RestartIO::ecl_type_get_type( ::Opm::RestartIO::ecl_kw_get_data_type( ecl_kw ) ) == ECL_DOUBLE_TYPE ) {
            //const double * ecl_data = ::Opm::RestartIO::ecl_kw_get_double_ptr( ecl_kw );
            const double * ecl_data = ::Opm::RestartIO::ecl_kw_get_type_ptr<double>( ecl_kw, ECL_DOUBLE_TYPE );
            return { ecl_data , ecl_data + size };
        } else {
            //const float * ecl_data = ::Opm::RestartIO::ecl_kw_get_float_ptr( ecl_kw );
	    const float * ecl_data = ::Opm::RestartIO::ecl_kw_get_type_ptr<float>( ecl_kw, ECL_FLOAT_TYPE );
            return { ecl_data , ecl_data + size };
        }

    }


  inline data::Solution restoreSOLUTION( ::Opm::RestartIO::ecl_file_view_type* file_view,
                                           const std::map<std::string, RestartKey>& keys,
                                           const UnitSystem& units,
                                           int numcells) {

        data::Solution sol;
        for (const auto& pair : keys) {
            const std::string& key = pair.first;
            UnitSystem::measure dim = pair.second.dim;
            bool required = pair.second.required;

            if( !::Opm::RestartIO::ecl_file_view_has_kw( file_view, key.c_str() ) ) {
                if (required)
                    throw std::runtime_error("Read of restart file: "
                                             "File does not contain "
                                             + key
                                             + " data" );
                else
                    continue;
            }

            const Opm::RestartIO::ecl_kw_type * ecl_kw = ::Opm::RestartIO::ecl_file_view_iget_named_kw( file_view , key.c_str() , 0 );
            if( Opm::RestartIO::ecl_kw_get_size(ecl_kw) != numcells)
                throw std::runtime_error("Restart file: Could not restore "
                                         + std::string( Opm::RestartIO::ecl_kw_get_header( ecl_kw ) )
                                         + ", mismatched number of cells" );

            std::vector<double> data = ::Opm::RestartIO::double_vector( ecl_kw );
            units.to_si( dim , data );

            sol.insert( key, dim, data , data::TargetType::RESTART_SOLUTION );
        }

        return sol;
    }


using rt = data::Rates::opt;
data::Wells restore_wells( const ::Opm::RestartIO::ecl_kw_type * opm_xwel,
                           const ::Opm::RestartIO::ecl_kw_type * opm_iwel,
                           int restart_step,
                           const EclipseState& es,
                           const EclipseGrid& grid,
			   const Schedule& schedule) {

    //const auto& sched_wells = es.getSchedule().getWells( restart_step );
    const auto& sched_wells = schedule.getWells( restart_step );
    std::vector< rt > phases;
    {
        const auto& phase = es.runspec().phases();
        if( phase.active( Phase::WATER ) ) phases.push_back( rt::wat );
        if( phase.active( Phase::OIL ) )   phases.push_back( rt::oil );
        if( phase.active( Phase::GAS ) )   phases.push_back( rt::gas );
    }

    const auto well_size = [&]( size_t acc, const Well* w ) {
        return acc
            + 2 + phases.size()
            + (  w->getCompletions( restart_step ).size()
              * (phases.size() + data::Completion::restart_size) );
    };

    const auto expected_xwel_size = std::accumulate( sched_wells.begin(),
                                                     sched_wells.end(),
                                                     0,
                                                     well_size );

    if( ::Opm::RestartIO::ecl_kw_get_size( opm_xwel ) != expected_xwel_size ) {
        throw std::runtime_error(
                "Mismatch between OPM_XWEL and deck; "
                "OPM_XWEL size was " + std::to_string( ::Opm::RestartIO::ecl_kw_get_size( opm_xwel ) ) +
                ", expected " + std::to_string( expected_xwel_size ) );
    }

    if( ::Opm::RestartIO::ecl_kw_get_size( opm_iwel ) != int(sched_wells.size()) )
        throw std::runtime_error(
                "Mismatch between OPM_IWEL and deck; "
                "OPM_IWEL size was " + std::to_string( ::Opm::RestartIO::ecl_kw_get_size( opm_iwel ) ) +
                ", expected " + std::to_string( sched_wells.size() ) );

    data::Wells wells;
    //const double * opm_xwel_data = ::Opm::RestartIO::ecl_kw_get_double_ptr( opm_xwel );
    const double * opm_xwel_data = ::Opm::RestartIO::ecl_kw_get_type_ptr<double>( opm_xwel, ECL_DOUBLE_TYPE );
    //const int * opm_iwel_data = ::Opm::RestartIO::ecl_kw_get_int_ptr( opm_iwel );
    const int * opm_iwel_data = ::Opm::RestartIO::ecl_kw_get_type_ptr<int>( opm_iwel, ECL_INT_TYPE );
    for( const auto* sched_well : sched_wells ) {
        data::Well& well = wells[ sched_well->name() ];

        well.bhp = *opm_xwel_data++;
        well.temperature = *opm_xwel_data++;
        well.control = *opm_iwel_data++;

        for( auto phase : phases )
            well.rates.set( phase, *opm_xwel_data++ );

        for( const auto& sc : sched_well->getCompletions( restart_step ) ) {
            const auto i = sc.getI(), j = sc.getJ(), k = sc.getK();
            if( !grid.cellActive( i, j, k ) || sc.getState() == WellCompletion::SHUT ) {
                opm_xwel_data += data::Completion::restart_size + phases.size();
                continue;
            }

            const auto active_index = grid.activeIndex( i, j, k );

            well.completions.emplace_back();
            auto& completion = well.completions.back();
            completion.index = active_index;
            completion.pressure = *opm_xwel_data++;
            completion.reservoir_rate = *opm_xwel_data++;
            for( auto phase : phases )
                completion.rates.set( phase, *opm_xwel_data++ );
        }
    }

    return wells;
}
//}


//* should take grid as argument because it may be modified from the simulator */
RestartValue load( const std::string& filename,
                   int report_step,
                   const std::map<std::string, RestartKey>& keys,
                   const EclipseState& es,
                   const EclipseGrid& grid,
		   const Schedule& schedule,
                   const std::map<std::string, bool>& extra_keys) {

    const bool unified                   = ( ::Opm::RestartIO::EclFiletype( filename ) == ::Opm::RestartIO::ECL_UNIFIED_RESTART_FILE );
    ::Opm::RestartIO::ert_unique_ptr< ::Opm::RestartIO::ecl_file_type, ::Opm::RestartIO::ecl_file_close > file(::Opm::RestartIO::ecl_file_open( filename.c_str(), 0 ));
    ::Opm::RestartIO::ecl_file_view_type * file_view;

    if( !file )
        throw std::runtime_error( "Restart file " + filename + " not found!" );

    if( unified ) {
        file_view = ::Opm::RestartIO::ecl_file_get_restart_view( file.get() , -1 , report_step , -1 , -1 );
        if (!file_view)
            throw std::runtime_error( "Restart file " + filename
                                      + " does not contain data for report step "
                                      + std::to_string( report_step ) + "!" );
    } else
        file_view = ::Opm::RestartIO::ecl_file_get_global_view( file.get() );

    const ::Opm::RestartIO::ecl_kw_type * intehead = ::Opm::RestartIO::ecl_file_view_iget_named_kw( file_view , "INTEHEAD", 0 );
    const ::Opm::RestartIO::ecl_kw_type * opm_xwel = ::Opm::RestartIO::ecl_file_view_iget_named_kw( file_view , "OPM_XWEL", 0 );
    const ::Opm::RestartIO::ecl_kw_type * opm_iwel = ::Opm::RestartIO::ecl_file_view_iget_named_kw( file_view, "OPM_IWEL", 0 );
    //UnitSystem units( static_cast<::Opm::RestartIO::ert_ecl_unit_enum>(::Opm::RestartIO::ecl_kw_iget_int( intehead , INTEHEAD_UNIT_INDEX )));
    //UnitSystem units( static_cast<::Opm::RestartIO::ert_ecl_unit_enum>(::Opm::RestartIO::ecl_kw_iget_type<int>( intehead , ECL_INT_TYPE, INTEHEAD_UNIT_INDEX )));
    UnitSystem units(getUnitSystem(intehead));
    RestartValue rst_value( ::Opm::RestartIO::restoreSOLUTION( file_view, keys, units , grid.getNumActive( )),
                            ::Opm::RestartIO::restore_wells( opm_xwel, opm_iwel, report_step , es, grid, schedule));

    for (const auto& pair : extra_keys) {
        const std::string& key = pair.first;
        bool required = pair.second;

        if (ecl_file_view_has_kw( file_view , key.c_str())) {
            const ::Opm::RestartIO::ecl_kw_type * ecl_kw = ::Opm::RestartIO::ecl_file_view_iget_named_kw( file_view , key.c_str() , 0 );
            //const double * data_ptr = ::Opm::RestartIO::ecl_kw_get_double_ptr( ecl_kw );
            const double * data_ptr = ::Opm::RestartIO::ecl_kw_get_type_ptr<double>( ecl_kw, ECL_DOUBLE_TYPE );
            const double * end_ptr  = data_ptr + ::Opm::RestartIO::ecl_kw_get_size( ecl_kw );
            rst_value.extra[ key ] = { data_ptr, end_ptr };
        } else if (required)
            throw std::runtime_error("No such key in file: " + key);

    }

    return rst_value;
}

//namespace {


std::vector<int> serialize_ICON( int report_step,
                                 int ncwmax,
                                 const std::vector<const Well*>& sched_wells) {

    size_t well_offset = 0;
    std::vector<int> data( sched_wells.size() * ncwmax * NICONZ , 0 );
    for (const Well* well : sched_wells) {
        const auto& completions = well->getCompletions( report_step );
        size_t completion_offset = 0;
        for( const auto& completion : completions) {
            size_t offset = well_offset + completion_offset;
            data[ offset + ICON_IC_INDEX ] = 1;

            data[ offset + ICON_I_INDEX ] = completion.getI() + 1;
            data[ offset + ICON_J_INDEX ] = completion.getJ() + 1;
            data[ offset + ICON_K_INDEX ] = completion.getK() + 1;
            data[ offset + ICON_DIRECTION_INDEX ] = completion.getDirection();
            {
                const auto open = WellCompletion::StateEnum::OPEN;
                data[ offset + ICON_STATUS_INDEX ] = completion.getState() == open
                    ? 1
                    : 0;
            }

            completion_offset += NICONZ;
        }
        well_offset += ncwmax * NICONZ;
    }
    return data;
}

std::vector<int> serialize_IWEL( size_t step,
                                 const std::vector<const Well *>& wells) {

    std::vector<int> data( wells.size() * NIWELZ , 0 );
    size_t offset = 0;
    for (const auto well : wells) {
        const auto& completions = well->getCompletions( step );

        data[ offset + IWEL_HEADI_INDEX ] = well->getHeadI( step ) + 1;
        data[ offset + IWEL_HEADJ_INDEX ] = well->getHeadJ( step ) + 1;
        data[ offset + IWEL_CONNECTIONS_INDEX ] = completions.size();
        data[ offset + IWEL_GROUP_INDEX ] = 1;

        data[ offset + IWEL_TYPE_INDEX ] = to_ert_welltype( *well, step );
        data[ offset + IWEL_STATUS_INDEX ] =
            well->getStatus( step ) == WellCommon::OPEN ? 1 : 0;

        offset += NIWELZ;
    }
    return data;
}










std::vector< int > serialize_OPM_IWEL( const data::Wells& wells,
                                       const std::vector< const Well* >& sched_wells ) {

    const auto getctrl = [&]( const Well* w ) {
        const auto itr = wells.find( w->name() );
        return itr == wells.end() ? 0 : itr->second.control;
    };

    std::vector< int > iwel( sched_wells.size(), 0.0 );
    std::transform( sched_wells.begin(), sched_wells.end(), iwel.begin(), getctrl );
    return iwel;
}

std::vector< double > serialize_OPM_XWEL( const data::Wells& wells,
                                          int report_step,
                                          const std::vector< const Well* >& sched_wells,
                                          const Phases& phase_spec,
                                          const EclipseGrid& grid ) {

    using rt = data::Rates::opt;

    std::vector< rt > phases;
    if( phase_spec.active( Phase::WATER ) ) phases.push_back( rt::wat );
    if( phase_spec.active( Phase::OIL ) )   phases.push_back( rt::oil );
    if( phase_spec.active( Phase::GAS ) )   phases.push_back( rt::gas );

    std::vector< double > xwel;
    for( const auto* sched_well : sched_wells ) {

        if( wells.count( sched_well->name() ) == 0 ) {
            const auto elems = (sched_well->getCompletions( report_step ).size()
                               * (phases.size() + data::Completion::restart_size))
                + 2 /* bhp, temperature */
                + phases.size();

            // write zeros if no well data is provided
            xwel.insert( xwel.end(), elems, 0.0 );
            continue;
        }

        const auto& well = wells.at( sched_well->name() );

        xwel.push_back( well.bhp );
        xwel.push_back( well.temperature );
        for( auto phase : phases )
            xwel.push_back( well.rates.get( phase ) );

        for( const auto& sc : sched_well->getCompletions( report_step ) ) {
            const auto i = sc.getI(), j = sc.getJ(), k = sc.getK();

            const auto rs_size = phases.size() + data::Completion::restart_size;
            if( !grid.cellActive( i, j, k ) || sc.getState() == WellCompletion::SHUT ) {
                xwel.insert( xwel.end(), rs_size, 0.0 );
                continue;
            }

            const auto active_index = grid.activeIndex( i, j, k );
            const auto at_index = [=]( const data::Completion& c ) {
                return c.index == active_index;
            };
            const auto& completion = std::find_if( well.completions.begin(),
                                                   well.completions.end(),
                                                   at_index );

            if( completion == well.completions.end() ) {
                xwel.insert( xwel.end(), rs_size, 0.0 );
                continue;
            }

            xwel.push_back( completion->pressure );
            xwel.push_back( completion->reservoir_rate );
            for( auto phase : phases )
                xwel.push_back( completion->rates.get( phase ) );
        }
    }

    return xwel;
};


std::vector<const char*> serialize_ZWEL( const std::vector<const Well *>& wells) {
    std::vector<const char*> data( wells.size( ) * NZWELZ , "");
    size_t offset = 0;

    for (const auto& well : wells) {
        data[ offset ] = well->name().c_str();
        offset += NZWELZ;
    }
    return data;
}

template< typename T >
void write_kw(::Opm::RestartIO::ecl_rst_file_type * rst_file , Opm::RestartIO::EclKW< T >&& kw) {
  ::Opm::RestartIO::ecl_rst_file_add_kw( rst_file, kw.get() );
}

void writeHeader(::Opm::RestartIO::ecl_rst_file_type * rst_file,
                 int report_step,
                 time_t posix_time,
                 double sim_days,
                 int ert_phase_mask,
                 const UnitSystem& units,
                 const Schedule& schedule,
                 const EclipseGrid& grid) {

    Opm::RestartIO::ecl_rsthead_type rsthead_data = {};

    rsthead_data.sim_time    = posix_time;
    rsthead_data.nactive     = grid.getNumActive();
    rsthead_data.nx          = grid.getNX();
    rsthead_data.ny          = grid.getNY();
    rsthead_data.nz          = grid.getNZ();
    rsthead_data.nwells      = schedule.numWells(report_step);
    rsthead_data.niwelz      = NIWELZ;
    rsthead_data.nzwelz      = NZWELZ;
    rsthead_data.niconz      = NICONZ;
    rsthead_data.ncwmax      = schedule.getMaxNumCompletionsForWells(report_step);
    rsthead_data.phase_sum   = ert_phase_mask;
    rsthead_data.sim_days    = sim_days;
    rsthead_data.unit_system = units.getEclType( );
    
    // this function should be moved to opm (does not need to be changed)
    ::Opm::RestartIO::ecl_util_set_date_values( rsthead_data.sim_time,
                              &rsthead_data.day,
                              &rsthead_data.month,
                              &rsthead_data.year );

    ::Opm::RestartIO::ecl_rst_file_fwrite_header( rst_file, report_step , &rsthead_data );
}


  ::Opm::RestartIO::ert_unique_ptr< ::Opm::RestartIO::ecl_kw_type, ecl_kw_free > ecl_kw( const std::string& kw, const std::vector<double>& data, bool write_double) {
    Opm::RestartIO::ert_unique_ptr< ::Opm::RestartIO::ecl_kw_type, ecl_kw_free > kw_ptr;

      if (write_double) {
	::Opm::RestartIO::ecl_kw_type * ecl_kw = ::Opm::RestartIO::ecl_kw_alloc( kw.c_str() , data.size() , ECL_DOUBLE );
	  ::Opm::RestartIO::ecl_kw_set_memcpy_data( ecl_kw , data.data() );
	  kw_ptr.reset( ecl_kw );
      } else {
	::Opm::RestartIO::ecl_kw_type * ecl_kw = ::Opm::RestartIO::ecl_kw_alloc( kw.c_str() , data.size() , ECL_FLOAT );
	  //float * float_data = ecl_kw_get_float_ptr( ecl_kw );
	  float * float_data = ecl_kw_get_type_ptr<float>( ecl_kw, ECL_FLOAT_TYPE );
	  for (size_t i=0; i < data.size(); i++)
	      float_data[i] = data[i];
	  kw_ptr.reset( ecl_kw );
      }

      return kw_ptr;
  }


  void writeSolution(::Opm::RestartIO::ecl_rst_file_type* rst_file, const data::Solution& solution, bool write_double) {
    ::Opm::RestartIO::ecl_rst_file_start_solution( rst_file );
    for (const auto& elm: solution) {
        if (elm.second.target == data::TargetType::RESTART_SOLUTION)
            ::Opm::RestartIO::ecl_rst_file_add_kw( rst_file , ecl_kw(elm.first, elm.second.data, write_double).get());
     }
     ::Opm::RestartIO::ecl_rst_file_end_solution( rst_file );

     for (const auto& elm: solution) {
        if (elm.second.target == data::TargetType::RESTART_AUXILIARY)
            ::Opm::RestartIO::ecl_rst_file_add_kw( rst_file , ecl_kw(elm.first, elm.second.data, write_double).get());
     }
  }


void writeExtraData(::Opm::RestartIO::ecl_rst_file_type* rst_file, const std::map<std::string,std::vector<double>>& extra_data) {
    for (const auto& pair : extra_data) {
        const std::string& key = pair.first;
        const std::vector<double>& data = pair.second;
        {
            ::Opm::RestartIO::ecl_kw_type * ecl_kw = ::Opm::RestartIO::ecl_kw_alloc_new_shared( key.c_str() , data.size() , ECL_DOUBLE , const_cast<double *>(data.data()));
            ::Opm::RestartIO::ecl_rst_file_add_kw( rst_file , ecl_kw);
            ::Opm::RestartIO::ecl_kw_free( ecl_kw );
        }
    }
}


void writeWell(::Opm::RestartIO::ecl_rst_file_type* rst_file, int report_step, const EclipseState& es , const EclipseGrid& grid, const Schedule& schedule, const data::Wells& wells) {

    //const auto& schedule = es.getSchedule();
    const auto sched_wells  = schedule.getWells(report_step);
    const auto& phases = es.runspec().phases();
    const size_t ncwmax = schedule.getMaxNumCompletionsForWells(report_step);

    const auto opm_xwel  = serialize_OPM_XWEL( wells, report_step, sched_wells, phases, grid );
    const auto opm_iwel  = serialize_OPM_IWEL( wells, sched_wells );
    const auto iwel_data = serialize_IWEL(report_step, sched_wells);
    const auto icon_data = serialize_ICON(report_step , ncwmax, sched_wells);
    const auto zwel_data = serialize_ZWEL( sched_wells );

    ::Opm::RestartIO::write_kw( rst_file, ::Opm::RestartIO::EclKW< int >( IWEL_KW, iwel_data) );
    ::Opm::RestartIO::write_kw( rst_file, ::Opm::RestartIO::EclKW< const char* >(ZWEL_KW, zwel_data ) );
    ::Opm::RestartIO::write_kw( rst_file, ::Opm::RestartIO::EclKW< double >( OPM_XWEL, opm_xwel ) );
    ::Opm::RestartIO::write_kw( rst_file, ::Opm::RestartIO::EclKW< int >( OPM_IWEL, opm_iwel ) );
    ::Opm::RestartIO::write_kw( rst_file, ::Opm::RestartIO::EclKW< int >( ICON_KW, icon_data ) );
}

void checkSaveArguments(const data::Solution& cells,
                        const EclipseGrid& grid,
                        const std::map<std::string, std::vector<double>>& extra_data) {

    const std::set<std::string> reserved_keys = {"LOGIHEAD", "INTEHEAD" ,"DOUBHEAD", "IWEL", "XWEL","ICON", "XCON" , "OPM_IWEL" , "OPM_XWEL", "ZWEL"};

    for (const auto& pair : extra_data) {
        const std::string& key = pair.first;
        if (key.size() > 8)
            throw std::runtime_error("The keys in extra data must have maximum eight characaters");

        if (cells.has( key ))
            throw std::runtime_error("The keys used must unique across Solution and extra_data");

        if (reserved_keys.find( key ) != reserved_keys.end())
            throw std::runtime_error("The extra_data uses a reserved key");
    }

    for (const auto& elm: cells)
        if (elm.second.data.size() != grid.getNumActive())
            throw std::runtime_error("Wrong size on solution vector: " + elm.first);
}


//}


void save(const std::string& filename,
          int report_step,
          double seconds_elapsed,
          data::Solution cells,
          data::Wells wells,
          const EclipseState& es,
          const EclipseGrid& grid,
	  const Schedule& schedule,
          std::map<std::string, std::vector<double>> extra_data,
	  bool write_double)
{
    ::Opm::RestartIO::checkSaveArguments( cells, grid, extra_data );
    {
        int ert_phase_mask = es.runspec().eclPhaseMask( );
        //const Schedule& schedule = es.getSchedule();
        const auto& units = es.getUnits();
        time_t posix_time = schedule.posixStartTime() + seconds_elapsed;
        const auto sim_time = units.from_si( UnitSystem::measure::time, seconds_elapsed );
        ::Opm::RestartIO::ert_unique_ptr< ::Opm::RestartIO::ecl_rst_file_type, ecl_rst_file_close > rst_file;

        if (::Opm::RestartIO::EclFiletype( filename ) == ECL_UNIFIED_RESTART_FILE)
	  rst_file.reset( ::Opm::RestartIO::ecl_rst_file_open_write_seek( filename.c_str(), report_step ) );
        else
	  rst_file.reset( ::Opm::RestartIO::ecl_rst_file_open_write( filename.c_str() ) );


        cells.convertFromSI( units );
        ::Opm::RestartIO::writeHeader(rst_file.get() , report_step, posix_time , sim_time, ert_phase_mask, units, schedule , grid );
        ::Opm::RestartIO::writeWell( rst_file.get() , report_step, es , grid, schedule, wells);
        ::Opm::RestartIO::writeSolution( rst_file.get() , cells , write_double );
        ::Opm::RestartIO::writeExtraData( rst_file.get() , extra_data );
    }
}

}
}