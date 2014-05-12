ProbeDeformer plugins for Maya
/**
 * @brief Cage Deformer plugins for Maya
 * @section LICENSE The MIT License
 * @section requirements:  Eigen 3:  http://eigen.tuxfamily.org/
 * @section Autodesk Maya: http://www.autodesk.com/products/autodesk-maya/overview
 * @section AffineLib: https://github.com/shizuo-kaji/AffineLib
 * @version 0.10
 * @date  3/Nov/2013
 * @author Shizuo KAJI
 */

There are two versions of deformers;
one is simple, and the other adds ARAP modifiation.
For the detail of the algorithm, refer to the paper 
"in preparation."

How to compile:
Look at the included Xcode project file.
For Windows users, please refer to Autodesk's web page.

How to use:
put the plugin in "MAYA_PLUG_IN_PATH"
put the python script in "MAYA_SCRIPT_PATH"
open script editor in Maya and type in the following Python command:
#
import plugin_deformer.ui_cageDeformer as ui
ui.UI_CageDeformer()
#




