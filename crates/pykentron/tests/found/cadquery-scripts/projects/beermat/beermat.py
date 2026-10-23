from cadquery import *

def make_beermat(mD: float, mTh: float, mThTol: float, mNumber: float) -> Workplane:
    offsetR = 3 
    excentricityD = 3
    nD = 2
    
    a = 5
    holder = (Workplane("YZ")
              .circle(mD/2+offsetR*2)
              .revolve(a, (0, nD * mD, 0), (1, nD * mD, 0)))
    holder += (Workplane("YZ")
               .move(0, - mD /2 - 2*offsetR)
               .hLine(mD /4 + offsetR / 2)
               .line(mD /4 + offsetR / 2, mD /2 + offsetR)
               .hLine(-mD /2 - offsetR)
               .mirrorY()
               .consolidateWires()
               .revolve(a, (0, nD * mD, 0), (1, nD * mD, 0)))

    
    holder -= (Workplane("YZ")
               .move(0, offsetR*2)
               .circle(mD/2 - offsetR)
               .revolve(360, (0, nD * mD, 0), (1, nD * mD, 0)))
    
    holder -= (Workplane("YZ")
               .move(0, offsetR*2)
               .circle(mD/2)
               .extrude(mNumber * (mTh + mThTol)/ 2.0, both=True))
    """
    holder -= (Workplane("XY", origin = (0, 0, 2 * offsetR))
               .rect(mNumber * (mTh + mThTol), mD)
               .workplane(offset=mD)
               .rect((mNumber)* (mTh + mThTol), mD)
               .loft())
    """
    
    holder -= (Workplane("XY", origin = (0, 0, 2 * offsetR))
               .rect(mNumber * (mTh + mThTol), mD)
               .extrude(mD))
    
    holder -= (Workplane("YZ")
               .line(mD, mD)
               .hLine(-mD*2)
               .close()
               .revolve(360, (0, nD * mD, 0), (1, nD * mD, 0))
               )

    holder -= Workplane("XY", origin= (0, 0, -mD/2)).rect(100, 100).extrude(-10)
    
    
    #holder = holder.faces(">Z").edges("not %LINE").fillet(1)
    
    #show_object(holder.edges("not (<Z or >X)"))
    
    holder = holder.edges("not (<Z or >X)").fillet(1.5)
    
    holder = holder.union(holder.mirror("YZ"), clean=False)
   
    
    return holder

if __name__ == "__cq_main__":
    mat_diameter = 108
    mat_thickness = 2.5
    mat_tol_th = 0.25
    n_mats = 8
    beermat = make_beermat(mat_diameter, mat_thickness, mat_tol_th, n_mats)
    show_object(beermat)
    exporters.export(beermat, "beermat.3mf")
    exporters.export(beermat, "beermat.stl")
    

if __name__ == "__main__":
    mat_diameter = 108
    mat_thickness = 2.5
    mat_tol_th = 0.25
    n_mats = 8
    beermat = make_beermat(mat_diameter, mat_thickness, mat_tol_th, n_mats)
    exporters.export(beermat, "beermat.3mf")
    exporters.export(beermat, "beermat.stl")