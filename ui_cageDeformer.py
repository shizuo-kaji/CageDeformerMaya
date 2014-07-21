# -*- coding: utf-8 -*-
# @file ui_cageDeformer.py
# @brief GUI for cageDeformer maya plugins
# @section LICENSE The MIT License
# @section requirements:  CageDeformer and CageDeformerARAP maya plugins
# @section Load from Maya's script editor:
# @section import plugin_deformer.ui_cageDeformer as ui; ui.UI_CageDeformer()
# @version 0.25
#  @author      Shizuo KAJI
#  @date        2013/5/13


#import debugmaya
#debugmaya.startDebug()

# Maya modules
import maya.cmds as cmds
import pymel.core as pm

try:
    cmds.loadPlugin("cageDeformer")
    cmds.loadPlugin("cageDeformerARAP")

except:
    print("Plugin already loaded")
    
class UI_CageDeformer:
    uiID = "CageDeformer"
    title = "CageDeformerPlugin"

    deformers = []
    
    def __init__(self):
        if pm.window(self.uiID, exists=True):
            pm.deleteUI(self.uiID)
        win = pm.window(self.uiID, title=self.title, menuBar=True)
        with win:
            pm.menu( label='CageDeformer', tearOff=True )
            pm.menuItem( label="CageDeformer", c=pm.Callback( self.initPlugin, "cage") )
            pm.menuItem( label="CageDeformerARAP", c=pm.Callback( self.initPlugin, "cageARAP") )
            self._parentLayout = pm.columnLayout( adj=True )
            with self._parentLayout:
                self.createUISet()

    def createUISet(self):
        self._childLayout = pm.columnLayout( adj=True )
        with self._childLayout:
            pm.text(l="Click cage mesh first, then shift+click target mesh, and click the menu item.")
            self.deformers = pm.ls(type="cage")
            for i in range(len(self.deformers)):
                frameLayout = pm.frameLayout( label=self.deformers[i].name(), collapsable = True)
                with frameLayout:
                    with pm.rowLayout(numberOfColumns=3) :
                        pm.button( l="Del", c=pm.Callback( self.deleteNode, self.deformers[i].name()))
                        pm.attrControlGrp( label="cage mode", attribute= self.deformers[i].cgm)
                        pm.attrControlGrp( label="blend mode", attribute= self.deformers[i].bm)
                    with pm.rowLayout(numberOfColumns=3) :
                        pm.attrControlGrp( label="rotation consistency", attribute= self.deformers[i].rc)
                        pm.attrControlGrp( label="Frechet sum", attribute= self.deformers[i].fs)                        
            self.deformers = pm.ls(type="cageARAP")
            for i in range(len(self.deformers)):
                frameLayout = pm.frameLayout( label=self.deformers[i].name(), collapsable = True)
                with frameLayout:
                    with pm.rowLayout(numberOfColumns=4) :
                        pm.button( l="Del", c=pm.Callback( self.deleteNode, self.deformers[i].name()))
                        pm.attrControlGrp( label="cage mode", attribute= self.deformers[i].cgm)
                        pm.attrControlGrp( label="blend mode", attribute= self.deformers[i].bm)
                        pm.attrControlGrp( label="constraint mode", attribute= self.deformers[i].constraintMode)
                    with pm.rowLayout(numberOfColumns=3) :
                        pm.attrControlGrp( label="rotation consistency", attribute= self.deformers[i].rc)
                        pm.attrControlGrp( label="Frechet sum", attribute= self.deformers[i].fs)                        
                        pm.attrControlGrp( label="constraint weight", attribute= self.deformers[i].constraintWeight)
                    
    def initPlugin(self, deformerType):
        meshes = pm.selected( type="transform" )
        if len(meshes)<2:
            return
        pm.select( meshes[-1])
        deformer = cmds.deformer(type=deformerType)[0]
        shape=meshes[-2].getShapes()[0]
        cmds.connectAttr(shape+".worldMesh[0]", deformer+".cageMesh")
        # Make deformer weights paintable
        cmds.makePaintable(deformerType, 'weights', attrType='multiFloat', shapeMode='deformer')
        self.updateUI()

    # delete deformer node
    def deleteNode(self,name):
        cmds.delete(name)
        self.updateUI()

    def updateUI(self):
        # update UI
        pm.deleteUI( self._childLayout )
        pm.setParent(self._parentLayout)
        self.createUISet()
                             