// SPDX-FileCopyrightText: Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
// SPDX-License-Identifier: BSD-3-Clause
#ifndef vtk_exodusII_h
#define vtk_exodusII_h

/* Use the libproj library configured for VTK.  */
#define VTK_MODULE_USE_EXTERNAL_vtkexodusII 1

#if VTK_MODULE_USE_EXTERNAL_vtkexodusII
# include <exodusII.h>
#else
# include <vtkexodusII/include/exodusII.h>
#endif

#endif
