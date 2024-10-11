/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Bo Zhang (University of Utah)
------------------------------------------------------------------------- */

#ifdef DUMP_CLASS
// clang-format off
DumpStyle(dspaces,DumpDSpaces);
// clang-format on
#else

#ifndef LMP_DUMP_DSPACES_H
#define LMP_DUMP_DSPACES_H

#include "dump_custom.h"
#include "dspaces.h"

namespace LAMMPS_NS {

class DumpDSpaces : public DumpCustom {
  public:
    DumpDSpaces(class LAMMPS *, int, char **);
    ~DumpDSpaces() {}
    static void finalize();

  private:
    //pack all meta data into a struct
    struct dspaces_lmp_meta { 
      double boxxlo;
      double boxxhi;
      double boxylo;
      double boxyhi;
      double boxzlo;
      double boxzhi;
      
      int triclinic;
      double boxxy;
      double boxxz;
      double boxyz;

      int boundary[3][2];
    };

    static dspaces_client_t client;

    // Use 2D array to store atom data in dspaces
    // (C array has to use the reverse order)
    // gdim[0] = # of field listed by user
    // gdim[1] = total # of atmos in a snapshot
    // Use user-defined dump id as the var name
    uint64_t dspaces_gdim[2];
    uint64_t dspaces_lb[2];
    uint64_t dspaces_ub[2];

    struct dspaces_lmp_meta meta;

    void write() override;
    void init_style() override;
};
}   // namespace LAMMPS_NS

#endif // LMP_DUMP_DSPACES_H
#endif // DUMP_CLASS