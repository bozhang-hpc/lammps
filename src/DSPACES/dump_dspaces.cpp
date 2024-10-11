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

#include "dump_dspaces.h"

#include "domain.h"
#include "error.h"
#include "fix.h"
#include "input.h"
#include "memory.h"
#include "modify.h"
#include "update.h"
#include "variable.h"

using namespace LAMMPS_NS;

dspaces_client_t DumpDSpaces::client = dspaces_CLIENT_NULL;

/* ---------------------------------------------------------------------- */

DumpDSpaces::DumpDSpaces(LAMMPS *lmp, int narg, char **arg) : DumpCustom(lmp, narg, arg)
{
  int dspaces_ret;

  try {
#if defined(MPI_STUBS)
    dspaces_ret = dspaces_init(0, &client);
#else
    dspaces_ret = dspaces_init_mpi(world, &client);
#endif
    if (dspaces_ret != dspaces_SUCCESS) throw dspaces_ret;
  } catch (int err) {
    error->all(FLERR, "Error: dspaces_init(), Error Code = {}", err);
  }
}

/* ---------------------------------------------------------------------- */

void DumpDSpaces::finalize()
{
  if(client != dspaces_CLIENT_NULL) {
    dspaces_fini(client);
  }
}

/* ---------------------------------------------------------------------- */

void DumpDSpaces::write()
{
  int dspaces_ret;

  // nme = # of dump lines this proc contributes to dump
  nme = count();

  // ntotal = total # of atoms in snapshot
  // atomOffset = sum of # of atoms up to this proc (exclusive prefix sum)
  bigint bnme = nme;
  MPI_Allreduce(&bnme, &ntotal, 1, MPI_LMP_BIGINT, MPI_SUM, world);

  bigint atomOffset;    // sum of all atoms on processes 0..me-1
  MPI_Scan(&bnme, &atomOffset, 1, MPI_LMP_BIGINT, MPI_SUM, world);
  atomOffset -= nme;    // exclusive prefix sum needed

  dspaces_gdim[0] = static_cast<uint64_t>(size_one);
  dspaces_gdim[1] = static_cast<uint64_t>(ntotal);

  dspaces_define_gdim(client, id, 2, dspaces_gdim);

  // ensure filewriter proc can receive everyone's info
  // limit nmax*size_one to int since used as arg in MPI_Rsend() below
  // pack my data into buf
  // if sorting on IDs also request ID list from pack()
  // sort buf as needed

  if (nme > maxbuf) {
    if ((bigint) nme * size_one > MAXSMALLINT) error->all(FLERR, "Too much per-proc info for dump");
    maxbuf = nme;
    memory->destroy(buf);
    memory->create(buf, (maxbuf * size_one), "dump:buf");
  }
  if (sort_flag && sortcol == 0 && nme > maxids) {
    maxids = nme;
    memory->destroy(ids);
    memory->create(ids, maxids, "dump:ids");
  }

  if (sort_flag && sortcol == 0)
    pack(ids);
  else
    pack(nullptr);
  if (sort_flag) sort();

  // dspaces meta-data = {
  //   .var_name = user-defined filename
  //   .version = current timestep
  //   .elem_size = all cells are in double
  //   .ndim = 2 (rows are atoms, columns are user-quired fields)
  //   .lb = lb[0] starts from the 1st user-quired field
  //         lb[1] is the atom offest where my rank starts to write
  //   .ub = ub[0] is the end index of the user-quired field
  //         ub[1] is the end index of the local atoms
  //   (C array has to use the reverse order for fast->slow dimensions)
  // }

  // We also need to put the # of total atoms as a metadata
  dspaces_ret = dspaces_put_meta(client, "natoms", update->ntimestep,
                                          (const void*)(&ntotal), sizeof(int64_t));
  if (dspaces_ret != dspaces_SUCCESS) {
    error->one(FLERR, "Error: dspaces_put_tag(), Error Code = {}", dspaces_ret);
  }                          

  // TODO: Use dspaces_put_tag() and set tag to the original lammps rank 
  dspaces_lb[0] = 0;
  dspaces_lb[1] = static_cast<size_t>(atomOffset);
  dspaces_ub[0] = static_cast<uint64_t>(size_one) - 1;
  dspaces_ub[1] = dspaces_lb[1] + static_cast<size_t>(nme) - 1;

  dspaces_ret = dspaces_put(client, filename, update->ntimestep, sizeof(double), 2, dspaces_lb, dspaces_ub, buf);
  if (dspaces_ret != dspaces_SUCCESS) {
    error->one(FLERR, "Error: dspaces_put_tag(), Error Code = {}", dspaces_ret);
  }
}

/* ---------------------------------------------------------------------- */

void DumpDSpaces::init_style()
{
  int dspaces_ret;

  // assemble column string from defaults and user values
  delete[] columns;
  std::string combined;
  int icol = 0;
  for (const auto &item : utils::split_words(columns_default)) {
    if (combined.size()) combined += " ";
    if (keyword_user[icol].size())
      combined += keyword_user[icol];
    else
      combined += item;
    ++icol;
  }
  columns = utils::strdup(combined);

  // remove % from filename since ADIOS always writes a global file with
  // data/metadata
  char *ptr = strchr(filename, '%');
  if (ptr) {
    while (*ptr) {
      ptr[0] = ptr[1];
      ++ptr;
    }
  }

  /* The next four loops are copied from dump_custom_mpiio, but nothing is
   * done with them.
   * It is unclear why we need them here.
   * For metadata, variable[] will be written out as an ADIOS attribute if
   * nvariable>0
   */
  // find current ptr for each compute,fix,variable
  // check that fix frequency is acceptable
  for (int i = 0; i < ncompute; i++) {
    compute[i] = modify->get_compute_by_id(id_compute[i]);
    if (!compute[i])
      error->all(FLERR, "Could not find dump dspaces compute ID {}", id_compute[i]);
  }

  for (int i = 0; i < nfix; i++) {
    fix[i] = modify->get_fix_by_id(id_fix[i]);
    if (!fix[i]) error->all(FLERR, "Could not find dump dspaces fix ID {}", id_fix[i]);
    if (nevery % fix[i]->peratom_freq)
      error->all(FLERR, "dump dspaces and fix {} with ID {} not computed at compatible times",
                 fix[i]->style, id_fix[i]);
  }

  int ivariable;
  for (int i = 0; i < nvariable; i++) {
    ivariable = input->variable->find(id_variable[i]);
    if (ivariable < 0) error->all(FLERR, "Could not find dump dspaces variable name");
    variable[i] = ivariable;
  }

  // set index and check validity of region
  if (idregion && !domain->get_region_by_id(idregion))
    error->all(FLERR, "Region {} for dump dspaces does not exist", idregion);

  // dspaces meta data preparation & put
  if (domain->triclinic == 0) {
    meta.triclinic = 0;
    meta.boxxlo = domain->boxlo[0];
    meta.boxxhi = domain->boxhi[0];
    meta.boxylo = domain->boxlo[1];
    meta.boxyhi = domain->boxhi[1];
    meta.boxzlo = domain->boxlo[2];
    meta.boxzhi = domain->boxhi[2];
  } else {
    meta.triclinic = 1;
    meta.boxxlo = domain->boxlo_bound[0];
    meta.boxxhi = domain->boxhi_bound[0];
    meta.boxylo = domain->boxlo_bound[1];
    meta.boxyhi = domain->boxhi_bound[1];
    meta.boxzlo = domain->boxlo_bound[2];
    meta.boxzhi = domain->boxhi_bound[2];
    meta.boxxy = domain->xy;
    meta.boxxz = domain->xz;
    meta.boxyz = domain->yz;
  }
  memcpy(meta.boundary, domain->boundary, 6*sizeof(int));

  std::string column_names;
  int column_names_str_size;

  /* Put the column names as the metadata*/
  // column_names.clear();
  // for(int i=0; i<nfield; i++) {
  //   column_names.append(earg[i]);
  //   if(i != nfield) {
  //     column_names.append(1, ' ');
  //   }
  // }
  column_names_str_size = strlen(columns);

  // Use user-defined dump id as the var name
  if (me==0) {
    dspaces_ret = dspaces_put_meta(client, "column_names_str_size", 0,
                                    &column_names_str_size, sizeof(int));
    if(dspaces_ret != dspaces_SUCCESS) {
      error->one(FLERR, "Error: dspaces_put_meta(column_names_size) failed, "
                        "Error Code = {}", dspaces_ret);
    }

    dspaces_ret = dspaces_put_meta(client, "column_names", 0,
                                    columns, column_names_str_size);
    if(dspaces_ret != dspaces_SUCCESS) {
      error->one(FLERR, "Error: dspaces_put_meta(column_names) failed, "
                        "Error Code = {}", dspaces_ret);
    }

    dspaces_ret = dspaces_put_meta(client, "ncolumns", 0, &nfield, sizeof(int));
    if(dspaces_ret != dspaces_SUCCESS) {
      error->one(FLERR, "Error: dspaces_put_meta(ncolumns) failed, "
                        "Error Code = {}", dspaces_ret);
    }

    dspaces_ret = dspaces_put_meta(client, id, 0, &meta, sizeof(struct dspaces_lmp_meta));
    if (dspaces_ret != dspaces_SUCCESS) {
      error->one(FLERR, "Error: dspaces_put_meta(dspaces_lmp_meta), "
                        "Error Code = {}", dspaces_ret);
    }
  }
}