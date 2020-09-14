#!/usr/bin/env bash
#   Copyright (C) 2016-2019 Ludwig Schneider
#
# This file is part of SOMA.
#
# SOMA is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# SOMA is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with SOMA.  If not, see <http://www.gnu.org/licenses/>.


# This file is outside of the building protocal, yet. Because we don't want to force
# 'gengetopt' as a dependency. If you as a developer change something at the cmdline
# interface run this script in the c source directory after you modified the "soma.ggo"
gengetopt -i soma.ggo -a som_args --set-version="hash-unknown-version-hash"
#Do some tricks, to get our git version
sed -i '1i#include <soma_config.h>' cmdline.h
sed -i s/"\"hash-unknown-version-hash\""/"get_soma_version()"/g cmdline.h
sed -i s/"\"hash-unknown-version-hash\""/"get_soma_version()"/g cmdline.c
#Don't exit, but return with result.
sed -i s/"exit ("/"return ("/g cmdline.h
sed -i s/"exit ("/"return ("/g cmdline.c
sed -i s/"EXIT_SUCCESS"/" 1 "/g cmdline.h
sed -i s/"EXIT_SUCCESS"/" 1 "/g cmdline.c
sed -i s/"EXIT_FAILURE"/" -1 "/g cmdline.h
sed -i s/"EXIT_FAILURE"/" -1 "/g cmdline.c
#Refromat the arg dumping to include a tab in front of each line
sed -i s/'%s=\\\"%s\\"'/'\\t%s=\\"%s\\"'/g cmdline.c
