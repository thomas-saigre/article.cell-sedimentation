# ##  CMakeLists.txt; coding: utf-8 ---

# Author(s): Thomas Saigre <thomas.saigre@cemosis.fr>
# Date: 26 Nov 2024
#
# Copyright (C) 2024 Feel++ Consortium
#
# Distributed under the GPL(GNU Public License):
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
#
#

feelpp_add_application(heatfluid_wss SRCS main.cpp
  LINK_LIBRARIES
  feelpp_toolbox_heatfluid_lib_2dP1G1_2dP2P1G1
  feelpp_toolbox_heatfluid_lib_2dP2G1_2dP2P1G1
  feelpp_toolbox_heatfluid_lib_2dP1G1_2dP1P1G1

  feelpp_toolbox_heatfluid_lib_3dP1G1_3dP2P1G1
  feelpp_toolbox_heatfluid_lib_3dP2G1_3dP2P1G1
  feelpp_toolbox_heatfluid_lib_3dP1G1_3dP1P1G1
)
