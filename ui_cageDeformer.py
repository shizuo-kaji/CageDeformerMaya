# -*- coding: utf-8 -*-
# @file ui_cageDeformer.py
# @brief GUI for cageDeformer maya plugins
# @section LICENSE The MIT License
# @section requirements:  CageDeformer and CageDeformerARAP maya plugins
# @section Load from Maya's script editor:
# @section import plugin_deformer.ui_cageDeformer as ui; ui.UI_CageDeformer()
# @version 0.25
#  @author      Shizuo KAJI
#  @date        13/May/2013


#import debugmaya
#debugmaya.startDebug()

# Maya modules
import maya.cmds as cmds
import pymel.core as pm

deformerTypes = ["cageDeformer","cageDeformerARAP"]

for type in deformerTypes:
    try:
        cmds.loadPlugin(type)
    except:
        print("Plugin %s already loaded" %(type))


class UI_CageDeformer:
    uiID = "CageDeformer"
    title = "CageDeformerPlugin"

    def __init__(self):
        if pm.window(self.uiID, exists=True):
            pm.deleteUI(self.uiID)
        win = pm.window(self.uiID, title=self.title, menuBar=True)
        with win:
            pm.menu( label='CageDeformer', tearOff=True )
            for type in deformerTypes:
                pm.menuItem( label=type, c=pm.Callback( self.initPlugin, type) )
            self._parentLayout = pm.columnLayout( adj=True )
            with self._parentLayout:
                self.createUISet()

    def createUISet(self):
        self._childLayout = pm.columnLayout( adj=True )
        with self._childLayout:
            pm.text(l="Click cage mesh first, then shift+click target mesh, and click the menu item.")
            # cageDeformer specific
            deformers = pm.ls(type=deformerTypes[0])
            for node in deformers:
                frameLayout = pm.frameLayout( label=node.name(), collapsable = True)
                with frameLayout:
                    self.createCommonAttr(node)

            # cageDeformerARAP specific
            deformers = pm.ls(type=deformerTypes[1])
            for node in deformers:
                frameLayout = pm.frameLayout( label=node.name(), collapsable = True)
                with frameLayout:
                    self.createCommonAttr(node)
                    with pm.rowLayout(numberOfColumns=3) :
                        pm.attrControlGrp( label="constraint mode", attribute= node.ctm)
                        pm.attrFieldSliderGrp( label="constraint weight", min=1e-10, max=1000, attribute=node.cw)
                        pm.attrFieldSliderGrp(label="constraint radius", min=0.001, max=10.0, attribute=node.cr)
                    with pm.rowLayout(numberOfColumns=3) :
                        pm.attrFieldSliderGrp( label="iteration", min=1, max=20, attribute=node.it)
                        pm.attrControlGrp( label="tet mode", attribute= node.tm)
                        pm.attrFieldSliderGrp( label="translation weight", min=0.0, max=1.0, attribute=node.tw)

    # initialise deformer
    def initPlugin(self, deformerType):
        meshes = pm.selected( type="transform" )
        if not meshes:
            return
        pm.select( meshes[-1])
        deformer = cmds.deformer(type=deformerType)[0]
        if len(meshes)>1:
            shape=meshes[-2].getShapes()[0]
            cmds.connectAttr(shape+".worldMesh[0]", deformer+".cageMesh")
        # Make deformer weights paintable
        cmds.makePaintable(deformerType, 'weights', attrType='multiFloat', shapeMode='deformer')
        self.updateUI()

    def setCage(self, node):
        meshes = pm.selected( type="transform" )
        if not meshes:
            return
        shape=meshes[-1].getShapes()[0]
        cmds.connectAttr(shape+".worldMesh[0]", node.name()+".cageMesh")
        self.updateUI()

    # draw common attributes
    def createCommonAttr(self,node):
        with pm.rowLayout(numberOfColumns=4):
            pm.button( l="Delete deformer", c=pm.Callback( self.deleteNode, node ))
            pm.button( l="Set selected as cage", c=pm.Callback( self.setCage, node))
            pm.attrControlGrp( label="cage mode", attribute= node.cgm)
            pm.attrControlGrp( label="blend mode", attribute= node.bm)
        with pm.rowLayout(numberOfColumns=4) :
            pm.attrControlGrp( label="normalise cage tet", attribute= node.nr)
            pm.attrControlGrp( label="area weight", attribute= node.aw)
            pm.attrControlGrp( label="positive weight", attribute= node.posw)
            pm.attrControlGrp( label="normalise weight", attribute= node.nw)
        with pm.rowLayout(numberOfColumns=3) :
            pm.attrControlGrp( label="rotation consistency", attribute= node.rc)
            pm.attrControlGrp( label="Frechet sum", attribute= node.fs)
            pm.attrControlGrp( label="symmetrise face", attribute= node.sf)
        with pm.rowLayout(numberOfColumns=4) :
            pm.attrControlGrp( label="Weight mode", attribute= node.wtm)
            pm.attrControlGrp( label="normExponent", attribute=node.ne)
            pm.attrFieldSliderGrp(label="effect radius", min=0.001, max=20.0, attribute=node.er)


    # delete deformer node
    def deleteNode(self,node):
        cmds.delete(node.name())
        self.updateUI()

    # update UI
    def updateUI(self):
        pm.deleteUI( self._childLayout )
        pm.setParent(self._parentLayout)
        self.createUISet()
