# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" Contains the results classes for the nrps_pks module """

import logging
from typing import Any, Dict, List, Optional, Union

from antismash.common.module_results import ModuleResults

class ClusterCompareResults(ModuleResults):
    """ The results of cluster comparison """
    _schema_version = 1

