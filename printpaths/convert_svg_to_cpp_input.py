from svgpathtools import svg2paths
from svgpathtools import Path, Line, QuadraticBezier, CubicBezier, Arc
import numpy as np
import sys,glob
import matplotlib.pyplot as plt


flipY = True # True for orchid and folding flower

def writeComplexPoint(fout, point, first=False):
    if(not first): fout.write("\t")
    if(flipY): fout.write("%10.10e %10.10e" % (np.real(point), -np.imag(point)))
    else: fout.write("%10.10e %10.10e" % (np.real(point), np.imag(point)))

fnames_i = sys.argv[1:]
for fname_i in fnames_i:
    fname_o = fname_i.replace('.svg', '.dat')

    print 'converting ',fname_i,' into ',fname_o

    paths, attributes = svg2paths(fname_i)
    nPaths = len(paths)

    with open(fname_o, 'w') as fout:
        fout.write("%d\n" % nPaths)
        for i,path in enumerate(paths):
            nSegs = len(path)
            if(nSegs==0): continue
            try:
                path_length = path.length()
            except: continue

            print 'path ',i,' has length ',path.length(),' and ',nSegs,' number of segments'
            fout.write("%d\n" % nSegs)
            for s,seg in enumerate(path):
                
                tStart = path.t2T(s, 0)
                tEnd = path.t2T(s, 1)
                fout.write("%10.10e \t %10.10e\n" % (tStart, tEnd))
                if(type(seg) == Line):
                    fout.write("1\n")
                    writeComplexPoint(fout, seg.start, True)
                    writeComplexPoint(fout, seg.end)
            #                print 'linear       : ',seg.start, seg.end, tStart, tEnd
                elif(type(seg) == QuadraticBezier):
                    fout.write("2\n")
                    writeComplexPoint(fout, seg.start, True)
                    writeComplexPoint(fout, seg.control1)
                    writeComplexPoint(fout, seg.end)
            #                print 'quad bezier  : ',seg.start, seg.control1, seg.end, tStart, tEnd
                elif(type(seg) == CubicBezier):
                    fout.write("3\n")
                    writeComplexPoint(fout, seg.start, True)
                    writeComplexPoint(fout, seg.control1)
                    writeComplexPoint(fout, seg.control2)
                    writeComplexPoint(fout, seg.end)
            #                print 'cubic bezier : ',seg.start, seg.control1, seg.control2, seg.end, tStart, tEnd
                elif(type(seg) == Arc):
                    fout.write("4\n")
                    writeComplexPoint(fout, seg.start, True)
                    writeComplexPoint(fout, seg.end)
                    writeComplexPoint(fout, seg.radius)
                    fout.write("\t %10.10e %d %d\n" % (seg.rotation, seg.large_arc, seg.sweep))
            #                print 'arc          : ',seg.start, seg.end, seg.radius, seg.rotation, seg.large_arc, seg.sweep, tStart, tEnd
                else:
                    print 'unknown segment type found      : ',seg,tStart, tEnd
                fout.write("\n")

            
