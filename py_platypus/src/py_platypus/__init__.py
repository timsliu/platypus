import py_platypus.models.pic_1d as pic_1d
import py_platypus.models.pic_2d as pic_2d
import py_platypus.models.pic_2d_em as pic_2d_em
import py_platypus.models.pic_3d as pic_3d
import py_platypus.utils.params as params
import py_platypus.utils.run_sim as run_sim
import py_platypus.utils.io_utils  as io_utils
import py_platypus.utils.math_utils  as math_utils
import py_platypus.utils.charge_step  as charge_step
import py_platypus.vis.plotter as plotter 
import py_platypus.vis.vis_util  as vis_util

import os
PLATYPUS_HOME = os.getenv("PLATYPUS_HOME")
DIM_MAP = {"x": 1, "y": 0, "z": 2}
