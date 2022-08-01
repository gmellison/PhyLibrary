function plotFigtree(treestring::Vector{String}, outfile="figtree.tre")  
    treemessage="#NEXUS
    begin trees;
        tree tree_1
    end;

    begin figtree;
        set appearance.backgroundColorAttribute=\"Default\";
        set appearance.backgroundColour=#ffffff;
        set appearance.branchColorAttribute=\"User selection\";
        set appearance.branchColorGradient=false;
        set appearance.branchLineWidth=3.0;
        set appearance.branchMinLineWidth=0.0;
        set appearance.branchWidthAttribute=\"Fixed\";
        set appearance.foregroundColour=#000000;
        set appearance.hilightingGradient=false;
        set appearance.selectionColour=#2d3680;
        set branchLabels.colorAttribute=\"User selection\";
        set branchLabels.displayAttribute=\"label\";
        set branchLabels.fontName=\"sansserif\";
        set branchLabels.fontSize=12;
        set branchLabels.fontStyle=0;
        set branchLabels.isShown=true;
        set branchLabels.significantDigits=4;
        set layout.expansion=100;
        set layout.layoutType=\"RECTILINEAR\";
        set layout.zoom=0;
        set legend.attribute=\"label\";
        set legend.fontSize=10.0;
        set legend.isShown=false;
        set legend.significantDigits=4;
        set nodeBars.barWidth=4.0;
        set nodeBars.displayAttribute=null;
        set nodeBars.isShown=false;
        set nodeLabels.colorAttribute=\"User selection\";
        set nodeLabels.displayAttribute=\"nodenumber\";
        set nodeLabels.fontName=\"sansserif\";
        set nodeLabels.fontSize=12;
        set nodeLabels.fontStyle=0;
        set nodeLabels.isShown=true;
        set nodeLabels.significantDigits=4;
        set nodeShapeExternal.colourAttribute=\"User selection\";
        set nodeShapeExternal.isShown=false;
        set nodeShapeExternal.minSize=10.0;
        set nodeShapeExternal.scaleType=Width;
        set nodeShapeExternal.shapeType=Circle;
        set nodeShapeExternal.size=4.0;
        set nodeShapeExternal.sizeAttribute=\"Fixed\";
        set nodeShapeInternal.colourAttribute=\"User selection\";
        set nodeShapeInternal.isShown=false;
        set nodeShapeInternal.minSize=10.0;
        set nodeShapeInternal.scaleType=Width;
        set nodeShapeInternal.shapeType=Circle;
        set nodeShapeInternal.size=4.0;
        set nodeShapeInternal.sizeAttribute=\"Fixed\";
        set polarLayout.alignTipLabels=false;
        set polarLayout.angularRange=0;
        set polarLayout.rootAngle=0;
        set polarLayout.rootLength=100;
        set polarLayout.showRoot=true;
        set radialLayout.spread=0.0;
        set rectilinearLayout.alignTipLabels=false;
        set rectilinearLayout.curvature=0;
        set rectilinearLayout.rootLength=100;
        set scale.offsetAge=0.0;
        set scale.rootAge=1.0;
        set scale.scaleFactor=1.0;
        set scale.scaleRoot=false;
        set scaleAxis.automaticScale=true;
        set scaleAxis.fontSize=8.0;
        set scaleAxis.isShown=false;
        set scaleAxis.lineWidth=1.0;
        set scaleAxis.majorTicks=1.0;
        set scaleAxis.minorTicks=0.5;
        set scaleAxis.origin=0.0;
        set scaleAxis.reverseAxis=false;
        set scaleAxis.showGrid=true;
        set scaleBar.automaticScale=true;
        set scaleBar.fontSize=10.0;
        set scaleBar.isShown=true;
        set scaleBar.lineWidth=1.0;
        set scaleBar.scaleRange=0.0;
        set tipLabels.colorAttribute=\"User selection\";
        set tipLabels.displayAttribute=\"Names\";
        set tipLabels.fontName=\"sansserif\";
        set tipLabels.fontSize=14;
        set tipLabels.fontStyle=1;
        set tipLabels.isShown=true;
        set tipLabels.significantDigits=4;
        set trees.order=true;
        set trees.orderType=\"increasing\";
        set trees.rooting=false;
        set trees.rootingType=\"User Selection\";
        set trees.transform=true;
        set trees.transformType=\"cladogram\";
    end;"    

    str = "tree tree_" .* string.(1:length(treestring)) .* " = " .* treestring
    str = join(str,"\n\t")
    treemessage1 = replace(treemessage, "tree tree_1" => str)

    out = open(outfile, "w")
    write(out, treemessage1)
    close(out)

    run(`figtree $outfile`)
end 
