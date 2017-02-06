ProbeDeformer plugins for Maya
/**
 * @brief Cage Deformer plugins for Maya
 * @section LICENSE The MIT License
 * @section requirements:  Eigen 3:  http://eigen.tuxfamily.org/
 * @section Autodesk Maya: http://www.autodesk.com/products/autodesk-maya/overview
 * @section (included) AffineLib: https://github.com/shizuo-kaji/AffineLib
 * @author Shizuo KAJI
 */

## The plugins
There are two versions of deformers;
one is simple, and the other adds ARAP modifiation.

## How to compile:
- Mac OS X: Look at the included Xcode project file
- Windows: Look at the included VisualStudio project file. __DO NOT__ turn on AVX instructions.
- Other: Look at the included Makefile
- on some systems, specifying the compiler option -DEIGEN_DONT_VECTORIZE may be necessary to avoid compilation errors (thank giordi91 for this information)

## How to use:
- put the plugins in "MAYA_PLUG_IN_PATH"
- put the UI python script in "MAYA_SCRIPT_PATH"
- open script editor in Maya and type in the following Python command:

```python
import ui_cageDeformer as ui
ui.UI_CageDeformer()
```

## LIMITATION:
The ARAP version works only on "clean" meshes.
First apply "Cleanup" from "Mesh" menu
to remove zero faces and edges, non-manifold geometry, etc.
